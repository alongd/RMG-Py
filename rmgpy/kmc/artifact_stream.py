"""Read a compiled event-set artifact one record at a time (I-295).

An emission shard used to parse the whole artifact (``json.loads`` of its
text) and certify every record, so its memory and time grew with the whole
catalog and more shards did not help: a 5.3 GB artifact failed every shard
inside the JSON parse.  This module lets a shard hold only what its own
records need.

``build_index`` streams the artifact once (memory: the index, one record and
one read block): every top-level member other than ``records`` (the *head*,
which binding and the certification rules read), and per record its event
id, family, byte span and *links* (the other records' event ids it carries as
string values: the strings ``compiler._replace_event_ids`` rewrites when
certification republishes a record).  The SHA-256 of the bytes is taken in
the same pass.

``read_records`` materialises the records of chosen entries by their byte
spans.  ``emit_prepare`` builds on it.
"""

from __future__ import annotations

import codecs
import hashlib
import json
import os
import tempfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Optional

INDEX_SCHEMA = "kmc-artifact-index/1"
BLOCK = 1 << 24
# a value is accepted once this many characters follow it (or at end of
# file): a number cut by the block boundary would otherwise decode short
LOOKAHEAD = 64
WHITESPACE = " \t\n\r"


class ArtifactFormatError(ValueError):
    """The artifact is not a JSON object with a ``records`` array of records."""


class IndexMismatch(ValueError):
    """An index does not belong to the artifact it is used with."""


@dataclass(frozen=True)
class Entry:
    event_id: str
    family: Optional[str]
    offset: int
    length: int
    links: tuple

    def as_list(self) -> list:
        return [self.event_id, self.family, self.offset, self.length, list(self.links)]

    @classmethod
    def from_list(cls, item) -> "Entry":
        event_id, family, offset, length, links = item
        return cls(event_id, family, int(offset), int(length), tuple(links))


@dataclass
class ArtifactIndex:
    sha256: str
    size: int
    head: dict
    entries: list

    def as_mapping(self) -> dict:
        return {"schema": INDEX_SCHEMA, "artifact_sha256": self.sha256,
                "artifact_bytes": self.size, "head": self.head,
                "entries": [e.as_list() for e in self.entries]}


def _byte_length(text: str) -> int:
    return len(text) if text.isascii() else len(text.encode("utf-8"))


class _Reader:
    """A JSON text consumed left to right from a binary file, with the byte
    offset of the read position and the digest of every byte read."""

    def __init__(self, handle, block: int):
        self.handle, self.block = handle, block
        self.decoder = codecs.getincrementaldecoder("utf-8")()
        self.json = json.JSONDecoder()
        self.digest = hashlib.sha256()
        self.size = 0
        self.text, self.pos, self.byte, self.eof = "", 0, 0, False

    def fill(self, at_least: int = 0) -> bool:
        """Append the next block (at least ``at_least`` bytes); False at end of file."""
        if self.eof:
            return False
        if self.pos:
            self.text, self.pos = self.text[self.pos:], 0
        raw = self.handle.read(max(self.block, at_least))
        if not raw:
            self.text += self.decoder.decode(b"", final=True)
            self.eof = True
            return False
        self.size += len(raw)
        self.digest.update(raw)
        self.text += self.decoder.decode(raw)
        return True

    def advance(self, end: int) -> None:
        self.byte += _byte_length(self.text[self.pos:end])
        self.pos = end

    def peek(self) -> str:
        """The next non-whitespace character ('' at end of file), not consumed."""
        while True:
            text, pos = self.text, self.pos
            while pos < len(text) and text[pos] in WHITESPACE:
                pos += 1
            self.advance(pos)
            if pos < len(text):
                return text[pos]
            if not self.fill():
                return ""

    def expect(self, chars: str) -> str:
        found = self.peek()
        if not found or found not in chars:
            raise ArtifactFormatError(f"expected one of {chars!r} at byte {self.byte}, "
                                      f"found {found!r}")
        self.advance(self.pos + 1)
        return found

    def value(self):
        """(value, first byte, end byte) of the JSON value at the position."""
        self.peek()
        while True:
            try:
                obj, end = self.json.raw_decode(self.text, self.pos)
            except json.JSONDecodeError as error:
                if self.fill(len(self.text) - self.pos):
                    continue
                raise ArtifactFormatError(f"unreadable value at byte {self.byte}: {error}")
            if self.eof or len(self.text) - end >= LOOKAHEAD:
                start = self.byte
                self.advance(end)
                return obj, start, self.byte
            self.fill(len(self.text) - self.pos)


def stream(path, block: int = BLOCK) -> Iterator[tuple]:
    """("member", key, value) for every top-level member but ``records``, and
    ("record", value, first byte, end byte) for every record, in file order;
    last ("end", sha256 hex digest, size in bytes), after every byte is read."""
    with open(path, "rb") as handle:
        reader = _Reader(handle, block)
        reader.expect("{")
        seen_records = False
        if reader.peek() == "}":
            reader.advance(reader.pos + 1)
        else:
            while True:
                key, start, _end = reader.value()
                if not isinstance(key, str):
                    raise ArtifactFormatError(f"a member name at byte {start} is not a string")
                reader.expect(":")
                if key == "records":
                    if seen_records:
                        raise ArtifactFormatError("two records members")
                    seen_records = True
                    reader.expect("[")
                    if reader.peek() == "]":
                        reader.advance(reader.pos + 1)
                    else:
                        while True:
                            record, first, end = reader.value()
                            yield "record", record, first, end
                            if reader.expect(",]") == "]":
                                break
                else:
                    value, _first, _end = reader.value()
                    yield "member", key, value
                if reader.expect(",}") == "}":
                    break
        if reader.peek():
            raise ArtifactFormatError(f"data after the artifact object at byte {reader.byte}")
        if not seen_records:
            raise ArtifactFormatError("the artifact has no records member")
        yield "end", reader.digest.hexdigest(), reader.size


def _strings(value) -> Iterator[str]:
    """Every string value (not member names) in ``value``: the strings
    ``compiler._replace_event_ids`` visits."""
    stack = [value]
    while stack:
        item = stack.pop()
        if isinstance(item, str):
            yield item
        elif isinstance(item, dict):
            stack.extend(item.values())
        elif isinstance(item, (list, tuple)):
            stack.extend(item)


def record_links(record: dict) -> set:
    """The event-id-shaped strings ``record`` carries other than its own id."""
    from rmgpy.kmc.event_record import _is_full_event_id
    own = record.get("event_id")
    return {s for s in _strings(record) if s != own and _is_full_event_id(s)}


def build_index(path, block: int = BLOCK) -> ArtifactIndex:
    """The index of the artifact at ``path``, from one streaming pass."""
    head, entries, raw_links = {}, [], []
    sha = size = None
    for item in stream(path, block):
        if item[0] == "member":
            head[item[1]] = item[2]
        elif item[0] == "record":
            record, start, end = item[1:]
            if not isinstance(record, dict) or not isinstance(record.get("event_id"), str):
                raise ArtifactFormatError(f"the record at byte {start} has no event id")
            entries.append((record["event_id"], record.get("family"), start, end - start))
            raw_links.append(record_links(record))
        else:
            sha, size = item[1], item[2]
    counts = Counter(e[0] for e in entries)
    duplicated = sorted(i for i, n in counts.items() if n > 1)
    if duplicated:
        raise ArtifactFormatError(f"refused:partition: two records of event id {duplicated[:3]}")
    return ArtifactIndex(sha, size, head, [
        Entry(event_id, family, start, length, tuple(sorted(links & counts.keys())))
        for (event_id, family, start, length), links in zip(entries, raw_links)])


def write_index(index: ArtifactIndex, path) -> Path:
    """Write ``index`` atomically (a failed write leaves no file)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, prefix=path.name + ".", suffix=".tmp")
    try:
        with os.fdopen(fd, "w") as handle:
            json.dump(index.as_mapping(), handle, sort_keys=True, separators=(",", ":"))
            handle.write("\n")
        os.replace(tmp, path)
    except BaseException:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise
    return path


def artifact_digest(path, block: int = BLOCK) -> tuple:
    """(sha256 hex digest, size) of the file, read in blocks."""
    digest, size = hashlib.sha256(), 0
    with open(path, "rb") as handle:
        while True:
            raw = handle.read(block)
            if not raw:
                return digest.hexdigest(), size
            digest.update(raw)
            size += len(raw)


def load_index(path, artifact) -> ArtifactIndex:
    """The index at ``path``, refused unless it is the index of ``artifact``'s
    bytes (size and SHA-256, re-read here)."""
    data = json.loads(Path(path).read_text())
    if not isinstance(data, dict) or data.get("schema") != INDEX_SCHEMA:
        raise IndexMismatch(f"refused:index: {path} is not a {INDEX_SCHEMA} index")
    sha, size = artifact_digest(artifact)
    if (data.get("artifact_sha256"), data.get("artifact_bytes")) != (sha, size):
        raise IndexMismatch(f"refused:index: {path} indexes another artifact "
                            f"({str(data.get('artifact_sha256'))[:16]}..., "
                            f"{data.get('artifact_bytes')} bytes) than {artifact} "
                            f"({sha[:16]}..., {size} bytes)")
    return ArtifactIndex(sha, size, data["head"], [Entry.from_list(e) for e in data["entries"]])


def index_for(artifact, index_path=None) -> ArtifactIndex:
    """The verified index of ``artifact``: read from ``index_path``, or built."""
    return build_index(artifact) if index_path is None else load_index(index_path, artifact)


def read_records(path, entries: Iterable[Entry]) -> list:
    """The records of ``entries`` read by their byte spans, in the order
    given; each must carry its entry's event id."""
    out = []
    with open(path, "rb") as handle:
        for entry in entries:
            handle.seek(entry.offset)
            record = json.loads(handle.read(entry.length).decode("utf-8"))
            if not isinstance(record, dict) or record.get("event_id") != entry.event_id:
                raise IndexMismatch(f"refused:index: the bytes at {entry.offset} are not "
                                    f"record {entry.event_id}")
            out.append(record)
    return out
