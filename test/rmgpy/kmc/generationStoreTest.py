"""Tests for the node-local generation store.

Scope: these tests prove behaviour under process interruption (SIGKILL of a
claimant or of a claimant's parent), concurrent processes, hostile filesystem
objects and malformed bytes.  They do NOT prove power-loss durability: an
``fsync`` that is recorded in the right order is evidence of intent, not of what
a disk does when power is cut.  Wrong file ownership cannot be produced without
privileges, so it is simulated by replacing the module's ``_fstat`` indirection.
"""

import errno
import json
import multiprocessing
import os
import signal
import stat
import subprocess
import sys
import textwrap
import time
from pathlib import Path
from types import SimpleNamespace

import pytest

from rmgpy.kmc import generation_identity as gi, generation_store as gs

SHA = "a" * 64
SHA2 = "b" * 64
SHA3 = "c" * 64
PRODUCER = "1" * 64
OTHER_PRODUCER = "2" * 64
CTX = multiprocessing.get_context("fork")
REPO_ROOT = str(Path(__file__).resolve().parents[3])


@pytest.fixture(autouse=True)
def hard_timeout():
    """Every test dies after 90 s instead of hanging the suite."""

    def expired(signum, frame):
        raise TimeoutError("test exceeded its hard timeout")

    previous = signal.signal(signal.SIGALRM, expired)
    signal.alarm(90)
    yield
    signal.alarm(0)
    signal.signal(signal.SIGALRM, previous)


@pytest.fixture(autouse=True)
def unmapped_root_is_root(monkeypatch):
    """In a user namespace the real root is reported as uid 65534; treat that owner as root.

    Only the ancestors above the store are owned by it; the shim changes nothing on a normal host.
    """
    if os.stat("/").st_uid != 65534 or os.getuid() == 65534:
        return
    real = gs._fstat

    def shim(fd):
        result = real(fd)
        return StatShim(result, 0) if result.st_uid == 65534 else result

    monkeypatch.setattr(gs, "_fstat", shim)


@pytest.fixture
def root(tmp_path):
    return str(tmp_path / "store")


def make_store(root, producer=PRODUCER, **kwargs):
    return gs.GenerationStore(root, producer_identity=producer, **kwargs)


def publish(store, sha, payload):
    with store.claim(sha, timeout=30) as claim:
        claim.publish(payload)


def entry_path(root, sha):
    return Path(root) / "v1" / "gen" / sha[:2] / f"{sha}.entry"


def quarantine_files(root):
    return sorted(p.name for p in (Path(root) / "v1" / "quarantine").iterdir())


def read_key(root):
    return (Path(root) / "v1" / "key").read_bytes()


def rewrite(path, data):
    temp = Path(str(path) + ".rewrite")
    fd = os.open(temp, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    with os.fdopen(fd, "wb") as handle:
        handle.write(data)
    os.replace(temp, path)


def signed_entry(root, sha, envelope_bytes):
    return gs.encode_entry(read_key(root), sha, envelope_bytes)


def good_envelope_bytes(sha=SHA, producer=PRODUCER, payload=b"hello"):
    return gi.envelope_to_bytes(gi.build_envelope(sha, producer, payload))


class StatShim:
    """A stat result whose st_uid is overridden."""

    def __init__(self, real, uid):
        self._real = real
        self.st_uid = uid

    def __getattr__(self, name):
        return getattr(self._real, name)


# ---------------------------------------------------------------------------
# Identity validation comes first
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("bad", ["../x", "A" * 64, "a" * 63, "a" * 65, "", "g" * 64, None, 5, b"a" * 64])
def test_invalid_request_sha_is_refused_before_any_filesystem_call(root, monkeypatch, bad):
    store = make_store(root)

    def spy(*args, **kwargs):
        raise AssertionError("filesystem touched before request_sha validation")

    for name in ("open", "mkdir", "stat", "lstat", "listdir", "scandir", "rename", "replace", "unlink", "fsync"):
        monkeypatch.setattr(os, name, spy)
    with pytest.raises(gs.InvalidIdentityError):
        store.peek(bad)
    with pytest.raises(gs.InvalidIdentityError):
        store.claim(bad)


@pytest.mark.parametrize("bad", ["x", "A" * 64, None, "a" * 63])
def test_producer_identity_is_mandatory_and_validated_without_touching_the_filesystem(root, monkeypatch, bad):
    with pytest.raises(TypeError):
        gs.GenerationStore(root)  # type: ignore[call-arg]

    def spy(*args, **kwargs):
        raise AssertionError("filesystem touched while validating the producer identity")

    monkeypatch.setattr(os, "open", spy)
    with pytest.raises(gs.InvalidIdentityError):
        gs.GenerationStore(root, producer_identity=bad)


def test_the_root_must_be_an_absolute_normalised_path(tmp_path):
    for bad in ("relative/store", str(tmp_path) + "/../x", str(tmp_path) + "/store/", ""):
        store = make_store(bad)
        assert isinstance(store.peek(SHA), gs.Refused)
        with pytest.raises(gs.UntrustedStoreError):
            store.claim(SHA)


# ---------------------------------------------------------------------------
# Root trust
# ---------------------------------------------------------------------------


def test_a_fresh_store_round_trips_and_has_the_contract_layout(root):
    store = make_store(root)
    assert store.peek(SHA) == gs.Miss()
    assert not os.path.exists(root), "peek must not create anything"
    publish(store, SHA, b"hello")
    assert store.peek(SHA) == gs.Hit(b"hello")
    base = Path(root) / "v1"
    for relative in ("gen", "locks", "quarantine", "key", "init.lock", f"gen/aa/{SHA}.entry", f"locks/{SHA}.lock"):
        assert (base / relative).exists(), relative
    for directory in (Path(root), base, base / "gen", base / "gen" / "aa", base / "locks", base / "quarantine"):
        assert stat.S_IMODE(directory.stat().st_mode) == 0o700, directory
    for file in (base / "key", base / "init.lock", entry_path(root, SHA), base / "locks" / f"{SHA}.lock"):
        assert stat.S_IMODE(file.stat().st_mode) == 0o600, file
    assert (base / "key").stat().st_size == 32
    raw = entry_path(root, SHA).read_bytes()
    assert raw.startswith(b"KMCGS1 ") and raw[71:72] == b"\n" and len(raw) > 72


def test_symlinked_root_is_refused(tmp_path):
    real = tmp_path / "real"
    real.mkdir(mode=0o700)
    link = tmp_path / "link"
    link.symlink_to(real)
    store = make_store(str(link))
    outcome = store.peek(SHA)
    assert isinstance(outcome, gs.Refused) and "symbolic link" in outcome.reason
    with pytest.raises(gs.UntrustedStoreError):
        store.claim(SHA)
    assert list(real.iterdir()) == []


def test_symlinked_ancestor_is_refused(tmp_path):
    real = tmp_path / "real"
    real.mkdir(mode=0o700)
    (real / "store").mkdir(mode=0o700)
    link = tmp_path / "link"
    link.symlink_to(real)
    store = make_store(str(link / "store"))
    assert isinstance(store.peek(SHA), gs.Refused)
    with pytest.raises(gs.UntrustedStoreError):
        store.claim(SHA)


def test_other_writable_non_sticky_ancestor_is_refused_and_sticky_is_accepted(tmp_path):
    open_dir = tmp_path / "open"
    open_dir.mkdir()
    (open_dir / "store").mkdir(mode=0o700)
    os.chmod(open_dir, 0o777)
    store = make_store(str(open_dir / "store"))
    assert isinstance(store.peek(SHA), gs.Refused)
    with pytest.raises(gs.UntrustedStoreError):
        store.claim(SHA)
    os.chmod(open_dir, 0o1777)
    sticky = make_store(str(open_dir / "store"))
    assert sticky.peek(SHA) == gs.Miss()
    group_writable = tmp_path / "grp"
    group_writable.mkdir()
    (group_writable / "store").mkdir(mode=0o700)
    os.chmod(group_writable, 0o770)
    assert isinstance(make_store(str(group_writable / "store")).peek(SHA), gs.Refused)


def test_uid_owned_sticky_world_writable_store_root_is_refused(root):
    publish(make_store(root), SHA, b"x")
    os.chmod(root, 0o1777)
    store = make_store(root)
    outcome = store.peek(SHA)
    assert isinstance(outcome, gs.Refused) and "0700" in outcome.reason
    with pytest.raises(gs.UntrustedStoreError):
        store.claim(SHA)


@pytest.mark.parametrize("relative", ["", "v1", "v1/gen", "v1/gen/aa", "v1/locks", "v1/quarantine"])
def test_internal_directory_with_mode_0755_is_refused_and_nothing_is_quarantined(root, relative):
    publish(make_store(root), SHA, b"x")
    rewrite(entry_path(root, SHA), b"junk that would be quarantined on a trusted root")
    target = Path(root) / relative if relative else Path(root)
    os.chmod(target, 0o755)
    store = make_store(root)
    outcome = store.peek(SHA)
    assert isinstance(outcome, gs.Refused) and "0700" in outcome.reason
    with pytest.raises(gs.UntrustedStoreError):
        store.claim(SHA)
    os.chmod(target, 0o700)
    assert entry_path(root, SHA).exists()
    assert quarantine_files(root) == []


def test_wrong_owner_is_refused_via_a_fstat_shim(root, monkeypatch):
    publish(make_store(root), SHA, b"x")
    victim = (Path(root) / "v1" / "gen").stat().st_ino
    real = gs._fstat

    def shim(fd):
        result = real(fd)
        if result.st_ino == victim:
            return StatShim(result, result.st_uid + 1)
        return result

    monkeypatch.setattr(gs, "_fstat", shim)
    store = make_store(root)
    outcome = store.peek(SHA)
    assert isinstance(outcome, gs.Refused) and "owner" in outcome.reason
    with pytest.raises(gs.UntrustedStoreError):
        store.claim(SHA)


def test_wrong_owner_ancestor_is_refused_via_a_fstat_shim(root, monkeypatch):
    publish(make_store(root), SHA, b"x")
    parent = Path(root).parent.stat().st_ino
    real = gs._fstat

    def shim(fd):
        result = real(fd)
        if result.st_ino == parent:
            return StatShim(result, os.getuid() + 12345)
        return result

    monkeypatch.setattr(gs, "_fstat", shim)
    assert isinstance(make_store(root).peek(SHA), gs.Refused)


MOUNTINFO = r"""
22 1 8:2 / / rw,relatime shared:1 - ext4 /dev/sda2 rw
23 22 0:25 / /home rw,nosuid - tmpfs tmpfs rw
24 22 0:26 / /mnt/with\040space rw,noatime master:2 shared:3 - xfs /dev/sdb rw,attr2
25 22 0:27 / /data rw - nfs4 server:/export rw
26 22 0:28 / /ov rw - overlay overlay rw,lowerdir=/a:/b
27 22 0:29 / /f rw - fuse.sshfs sshfs rw
28 22 0:30 / /vault rw - btrfs /dev/sdc rw
29 22 0:31 / /vault/inner/sub rw - ext4 /dev/sdd rw
30 22 0:32 / /mnt/tab\011bed rw - ext4 /dev/sde rw
31 22 0:33 / /homework rw - ext4 /dev/sdf rw
"""


def test_mountinfo_parsing_handles_optional_fields_and_octal_escapes():
    mounts = gs.parse_mountinfo(MOUNTINFO)
    by_point = {m.mount_point: m.fstype for m in mounts}
    assert by_point["/"] == "ext4"
    assert by_point["/mnt/with space"] == "xfs"
    assert by_point["/mnt/tab\tbed"] == "ext4"
    assert by_point["/f"] == "fuse.sshfs"
    assert len(mounts) == 10


@pytest.mark.parametrize(
    "root,accepted",
    [
        ("/srv/store", True),
        ("/mnt/with space/store", True),
        ("/vault/store", True),
        ("/homework/store", True),
        ("/home/store", False),
        ("/data/store", False),
        ("/ov/store", False),
        ("/f/store", False),
        ("/vault", False),
    ],
)
def test_mountinfo_filesystem_allowlist(root, accepted):
    mounts = gs.parse_mountinfo(MOUNTINFO)
    reason = gs.filesystem_refusal(root, mounts)
    assert (reason is None) == accepted, reason
    if not accepted and root != "/vault":
        assert any(name in reason for name in ("tmpfs", "nfs4", "overlay", "fuse.sshfs"))


def test_a_mount_point_strictly_inside_the_root_refuses_the_root():
    mounts = gs.parse_mountinfo(MOUNTINFO)
    reason = gs.filesystem_refusal("/vault", mounts)
    assert reason is not None and "inside" in reason
    assert gs.filesystem_refusal("/vault/inner", mounts) is not None


def test_unreadable_or_empty_mountinfo_refuses_the_root():
    assert gs.filesystem_refusal("/srv/store", gs.parse_mountinfo("")) is not None


def test_a_root_on_a_disallowed_filesystem_is_refused_end_to_end(root, monkeypatch):
    publish(make_store(root), SHA, b"x")
    text = f"1 0 0:1 / / rw - ext4 /dev/x rw\n2 1 0:2 / {Path(root).parent} rw - tmpfs tmpfs rw\n"
    monkeypatch.setattr(gs, "_read_mountinfo", lambda: text)
    store = make_store(root)
    outcome = store.peek(SHA)
    assert isinstance(outcome, gs.Refused) and "tmpfs" in outcome.reason
    with pytest.raises(gs.UntrustedStoreError, match="tmpfs"):
        store.claim(SHA)


# ---------------------------------------------------------------------------
# Key
# ---------------------------------------------------------------------------


def _init_worker(root, sha, barrier):
    store = make_store(root)
    barrier.wait(30)
    publish(store, sha, b"payload-" + sha[:4].encode())


def test_two_processes_initialising_at_once_end_with_one_32_byte_key(root):
    Path(root).mkdir(mode=0o700)
    barrier = CTX.Barrier(2)
    procs = [CTX.Process(target=_init_worker, args=(root, sha, barrier)) for sha in (SHA, SHA2)]
    for proc in procs:
        proc.start()
    for proc in procs:
        proc.join(60)
        assert proc.exitcode == 0
    assert len(read_key(root)) == 32
    assert not (Path(root) / "v1" / "key.tmp").exists()
    fresh = make_store(root)
    assert fresh.peek(SHA) == gs.Hit(b"payload-aaaa")
    assert fresh.peek(SHA2) == gs.Hit(b"payload-bbbb")


def test_stale_key_tmp_is_removed_under_the_init_lock(root):
    Path(root).mkdir(mode=0o700)
    (Path(root) / "v1").mkdir(mode=0o700)
    outside = Path(root).parent / "outside"
    outside.write_bytes(b"precious")
    (Path(root) / "v1" / "key.tmp").symlink_to(outside)
    publish(make_store(root), SHA, b"x")
    assert not (Path(root) / "v1" / "key.tmp").exists()
    assert outside.read_bytes() == b"precious"
    assert len(read_key(root)) == 32


def test_key_failures_are_refused_and_never_quarantine(root):
    publish(make_store(root), SHA, b"x")
    key = Path(root) / "v1" / "key"
    original = key.read_bytes()

    def reset():
        if key.is_symlink() or key.exists():
            key.unlink()
        for extra in Path(root).parent.glob("alias*"):
            extra.unlink()
        rewrite(key, original)

    alias = Path(root).parent / "alias"
    os.link(key, alias)
    outcome = make_store(root).peek(SHA)
    assert isinstance(outcome, gs.Refused) and "key" in outcome.reason
    with pytest.raises(gs.UntrustedStoreError):
        make_store(root).claim(SHA)
    reset()
    assert make_store(root).peek(SHA) == gs.Hit(b"x")

    rewrite(key, original[:31])
    assert isinstance(make_store(root).peek(SHA), gs.Refused)
    with pytest.raises(gs.UntrustedStoreError):
        make_store(root).claim(SHA)

    rewrite(key, original + b"\0")
    assert isinstance(make_store(root).peek(SHA), gs.Refused)

    rewrite(key, original)
    os.chmod(key, 0o644)
    assert isinstance(make_store(root).peek(SHA), gs.Refused)
    os.chmod(key, 0o600)

    key.unlink()
    key.symlink_to(Path(root).parent / "elsewhere")
    (Path(root).parent / "elsewhere").write_bytes(original)
    assert isinstance(make_store(root).peek(SHA), gs.Refused)
    key.unlink()
    rewrite(key, original)
    assert quarantine_files(root) == []
    assert entry_path(root, SHA).exists()
    assert make_store(root).peek(SHA) == gs.Hit(b"x")


def test_a_replaced_key_makes_old_entries_unauthenticated(root):
    publish(make_store(root), SHA, b"x")
    key = Path(root) / "v1" / "key"
    rewrite(key, bytes(32))
    assert isinstance(make_store(root).peek(SHA), gs.Refused)


# ---------------------------------------------------------------------------
# Corruption and hostile entries
# ---------------------------------------------------------------------------


def _seed(root, payload=b"hello", store=None):
    store = store or make_store(root)
    publish(store, SHA, payload)
    return entry_path(root, SHA).read_bytes()


def _flip(data):
    index = len(data) - 20
    return data[:index] + bytes([data[index] ^ 0x01]) + data[index + 1 :]


def _surrogate_envelope():
    return (
        b'{"payload":"\\ud800","payload_sha256":"' + b"0" * 64 + b'","producer_identity":"' + PRODUCER.encode()
        + b'","request_sha":"' + SHA.encode() + b'","schema_version":1}'
    )


def _duplicate_envelope():
    good = good_envelope_bytes()
    return good.replace(b'"schema_version":1', b'"schema_version":1,"schema_version":1')


CORRUPTIONS = {
    "truncated": lambda data, root: data[: len(data) // 2],
    "truncated inside header": lambda data, root: data[:30],
    "empty": lambda data, root: b"",
    "flipped byte": lambda data, root: _flip(data),
    "bad hmac": lambda data, root: data[:7] + b"0" * 64 + data[71:],
    "bad magic": lambda data, root: b"XXXXXXX" + data[7:],
    "no newline": lambda data, root: data[:71] + b"X" + data[72:],
    "valid hmac, duplicate keys": lambda data, root: signed_entry(root, SHA, _duplicate_envelope()),
    "valid hmac, lone surrogate": lambda data, root: signed_entry(root, SHA, _surrogate_envelope()),
    "valid hmac, oversize": lambda data, root: signed_entry(root, SHA, good_envelope_bytes(payload=b"z" * 5000)),
    "valid hmac, nesting bomb": lambda data, root: signed_entry(root, SHA, b"[" * 5000 + b"]" * 5000),
    "valid hmac, not an object": lambda data, root: signed_entry(root, SHA, b"[]"),
    "valid hmac, bound to another request": lambda data, root: signed_entry(
        root, SHA, good_envelope_bytes(sha=SHA2)
    ),
}


@pytest.mark.parametrize("name", sorted(CORRUPTIONS))
def test_corrupt_entry_is_refused_by_peek_and_quarantined_under_a_claim(root, name):
    store = make_store(root, max_entry_bytes=2048)
    data = _seed(root, store=store)
    corrupt = CORRUPTIONS[name](data, root)
    rewrite(entry_path(root, SHA), corrupt)

    reader = make_store(root, max_entry_bytes=2048)
    outcome = reader.peek(SHA)
    assert isinstance(outcome, gs.Refused) and outcome.reason
    assert entry_path(root, SHA).exists(), "peek must not mutate"
    assert quarantine_files(root) == []

    with reader.claim(SHA, timeout=30) as claim:
        assert isinstance(claim.outcome, gs.Refused)
        names = quarantine_files(root)
        assert len(names) == 2
        moved = [n for n in names if not n.endswith(".reason")][0]
        assert moved.startswith(SHA + ".") and moved.endswith(f".{os.getpid()}")
        assert (Path(root) / "v1" / "quarantine" / moved).read_bytes() == corrupt
        reason = (Path(root) / "v1" / "quarantine" / (moved + ".reason")).read_text()
        assert reason.strip()
        assert not entry_path(root, SHA).exists()
        claim.publish(b"replacement")
    assert reader.peek(SHA) == gs.Hit(b"replacement")


def test_entry_copied_under_another_sha_is_refused(root):
    data = _seed(root)
    store = make_store(root)
    publish(store, SHA2, b"other")
    rewrite(entry_path(root, SHA2), data)
    outcome = make_store(root).peek(SHA2)
    assert isinstance(outcome, gs.Refused)
    with make_store(root).claim(SHA2, timeout=30) as claim:
        assert isinstance(claim.outcome, gs.Refused)
        assert len(quarantine_files(root)) == 2
    assert make_store(root).peek(SHA) == gs.Hit(b"hello")


def test_entry_for_another_producer_is_refused_and_quarantined(root):
    _seed(root)
    other = make_store(root, producer=OTHER_PRODUCER)
    outcome = other.peek(SHA)
    assert isinstance(outcome, gs.Refused) and "producer" in outcome.reason
    with other.claim(SHA, timeout=30) as claim:
        assert isinstance(claim.outcome, gs.Refused)
        assert len(quarantine_files(root)) == 2
    assert make_store(root).peek(SHA) == gs.Miss()


def test_symlink_entry_to_a_valid_authenticated_entry_is_refused(root):
    data = _seed(root)
    target = Path(root).parent / "target.entry"
    rewrite(target, data)
    entry = entry_path(root, SHA)
    entry.unlink()
    entry.symlink_to(target)
    outcome = make_store(root).peek(SHA)
    assert isinstance(outcome, gs.Refused) and "symbolic link" in outcome.reason
    with make_store(root).claim(SHA, timeout=30) as claim:
        assert isinstance(claim.outcome, gs.Refused)
        assert len(quarantine_files(root)) == 2
    assert target.read_bytes() == data, "the symlink target must be untouched"
    assert make_store(root).peek(SHA) == gs.Miss()


def test_fifo_entry_is_refused_without_blocking(root):
    _seed(root)
    entry = entry_path(root, SHA)
    entry.unlink()
    os.mkfifo(entry, 0o600)
    started = time.monotonic()
    outcome = make_store(root).peek(SHA)
    assert isinstance(outcome, gs.Refused) and "regular" in outcome.reason
    with make_store(root).claim(SHA, timeout=30) as claim:
        assert isinstance(claim.outcome, gs.Refused)
        assert len(quarantine_files(root)) == 2
    assert time.monotonic() - started < 10


def test_hard_linked_entry_is_refused_and_the_other_link_survives(root):
    data = _seed(root)
    alias = Path(root).parent / "alias.entry"
    os.link(entry_path(root, SHA), alias)
    outcome = make_store(root).peek(SHA)
    assert isinstance(outcome, gs.Refused) and "link" in outcome.reason
    with make_store(root).claim(SHA, timeout=30) as claim:
        assert isinstance(claim.outcome, gs.Refused)
    assert alias.read_bytes() == data
    assert len(quarantine_files(root)) == 2


def test_entry_with_the_wrong_mode_is_refused(root):
    _seed(root)
    os.chmod(entry_path(root, SHA), 0o644)
    outcome = make_store(root).peek(SHA)
    assert isinstance(outcome, gs.Refused) and "0600" in outcome.reason


def test_a_directory_in_place_of_an_entry_is_refused(root):
    _seed(root)
    entry = entry_path(root, SHA)
    entry.unlink()
    entry.mkdir(mode=0o700)
    assert isinstance(make_store(root).peek(SHA), gs.Refused)


@pytest.mark.skipif(os.getuid() == 0, reason="root can read mode 0000 files")
def test_unexpected_errors_propagate_instead_of_becoming_a_miss(root):
    _seed(root)
    os.chmod(entry_path(root, SHA), 0)
    try:
        with pytest.raises(PermissionError):
            make_store(root).peek(SHA)
        with pytest.raises(PermissionError):
            make_store(root).claim(SHA, timeout=30).__enter__()
    finally:
        os.chmod(entry_path(root, SHA), 0o600)


def test_quarantine_rechecks_under_the_claim_and_spares_a_valid_replacement(root, monkeypatch):
    data = _seed(root)
    entry = entry_path(root, SHA)
    rewrite(entry, _flip(data))
    seen = []

    def replace_with_valid(claim):
        seen.append(claim)
        rewrite(entry, data)

    monkeypatch.setattr(gs, "_before_quarantine", replace_with_valid)
    with make_store(root).claim(SHA, timeout=30) as claim:
        assert claim.outcome == gs.Hit(b"hello")
    assert seen, "the quarantine path was not reached"
    assert quarantine_files(root) == []
    assert entry.read_bytes() == data


def test_untouched_valid_entries_are_never_quarantined(root):
    _seed(root)
    for _ in range(3):
        with make_store(root).claim(SHA, timeout=30) as claim:
            assert claim.outcome == gs.Hit(b"hello")
    assert quarantine_files(root) == []


# ---------------------------------------------------------------------------
# Publication
# ---------------------------------------------------------------------------


def _bucket_listing(root, sha=SHA):
    return sorted(p.name for p in (Path(root) / "v1" / "gen" / sha[:2]).iterdir())


def _failing_fsync(kind, errno_value=errno.EIO, after_replace=None):
    real_fstat = os.fstat
    real_fsync = os.fsync

    def shim(fd):
        mode = real_fstat(fd).st_mode
        wanted = stat.S_ISREG(mode) if kind == "file" else stat.S_ISDIR(mode)
        if wanted and (after_replace is None or after_replace["done"]):
            raise OSError(errno_value, "injected fsync failure")
        return real_fsync(fd)

    return shim


def _entry_bytes_for(root, sha, payload, producer=PRODUCER):
    return signed_entry(root, sha, gi.envelope_to_bytes(gi.build_envelope(sha, producer, payload)))


def test_fsync_file_failure_leaves_no_entry_and_preserves_the_previous_one(root, monkeypatch):
    store = make_store(root)
    publish(store, SHA, b"first")
    replacement = _entry_bytes_for(root, SHA, b"second")
    with store.claim(SHA, timeout=30) as claim:
        monkeypatch.setattr(gs, "_fsync", _failing_fsync("file"))
        with pytest.raises(OSError, match="injected"):
            claim._install(replacement)
        monkeypatch.undo()
    assert store.peek(SHA) == gs.Hit(b"first")
    assert _bucket_listing(root) == [f"{SHA}.entry"], "the temp file must be removed"


def test_fsync_file_failure_on_a_first_publication_leaves_no_entry(root, monkeypatch):
    store = make_store(root)
    with store.claim(SHA, timeout=30) as claim:
        monkeypatch.setattr(gs, "_fsync", _failing_fsync("file"))
        with pytest.raises(OSError, match="injected"):
            claim.publish(b"first")
        monkeypatch.undo()
        assert not isinstance(claim.outcome, gs.Hit)
    assert store.peek(SHA) == gs.Miss()
    assert _bucket_listing(root) == []


def test_rename_failure_preserves_the_previous_entry_and_removes_the_temp(root, monkeypatch):
    store = make_store(root)
    publish(store, SHA, b"first")
    replacement = _entry_bytes_for(root, SHA, b"second")

    def refuse(*args, **kwargs):
        raise OSError(errno.EXDEV, "injected rename failure")

    with store.claim(SHA, timeout=30) as claim:
        monkeypatch.setattr(gs, "_replace", refuse)
        with pytest.raises(OSError, match="injected"):
            claim._install(replacement)
        monkeypatch.undo()
    assert store.peek(SHA) == gs.Hit(b"first")
    assert _bucket_listing(root) == [f"{SHA}.entry"]


def test_post_rename_directory_fsync_failure_is_uncertain_and_deletes_nothing(root, monkeypatch):
    store = make_store(root)
    state = {"done": False}
    real_replace = os.replace

    def tracking_replace(*args, **kwargs):
        real_replace(*args, **kwargs)
        state["done"] = True

    with store.claim(SHA, timeout=30) as claim:
        monkeypatch.setattr(gs, "_replace", tracking_replace)
        monkeypatch.setattr(gs, "_fsync", _failing_fsync("dir", after_replace=state))
        with pytest.raises(gs.PublicationUncertainError) as caught:
            claim.publish(b"maybe")
        monkeypatch.undo()
        assert isinstance(caught.value.__cause__, OSError)
    assert state["done"]
    assert store.peek(SHA) == gs.Hit(b"maybe")
    assert _bucket_listing(root) == [f"{SHA}.entry"]


def test_publication_orders_file_fsync_before_rename_before_bucket_fsync(root, monkeypatch):
    store = make_store(root)
    log = []
    real_fsync, real_replace = os.fsync, os.replace

    def record_fsync(fd):
        info = os.fstat(fd)
        log.append(("fsync", info.st_ino, stat.S_ISREG(info.st_mode)))
        return real_fsync(fd)

    def record_replace(src, dst, *, src_dir_fd=None, dst_dir_fd=None):
        log.append(("replace", os.fstat(dst_dir_fd).st_ino, False))
        return real_replace(src, dst, src_dir_fd=src_dir_fd, dst_dir_fd=dst_dir_fd)

    monkeypatch.setattr(gs, "_fsync", record_fsync)
    monkeypatch.setattr(gs, "_replace", record_replace)
    publish(store, SHA, b"ordered")
    monkeypatch.undo()

    base = Path(root) / "v1"
    entry_ino = entry_path(root, SHA).stat().st_ino
    bucket_ino = (base / "gen" / SHA[:2]).stat().st_ino
    file_fsync = [i for i, e in enumerate(log) if e == ("fsync", entry_ino, True)]
    renames = [i for i, e in enumerate(log) if e[0] == "replace"]
    bucket_fsync = [i for i, e in enumerate(log) if e == ("fsync", bucket_ino, False)]
    assert len(file_fsync) == 1 and len(renames) == 1 and bucket_fsync, log
    assert log[renames[0]][1] == bucket_ino
    assert file_fsync[0] < renames[0] < bucket_fsync[-1]
    for directory in (base, base / "gen"):
        assert ("fsync", directory.stat().st_ino, False) in log, f"{directory} was not fsynced"
    first_publish_fsyncs = [i for i, e in enumerate(log) if e[0] == "fsync" and not e[2]]
    assert first_publish_fsyncs, "directory fsyncs missing"


def test_every_publisher_fsyncs_the_internal_chain_even_when_dirs_exist(root, monkeypatch):
    store = make_store(root)
    publish(store, SHA, b"one")
    log = []
    real_fsync = os.fsync

    def record(fd):
        info = os.fstat(fd)
        log.append(info.st_ino)
        return real_fsync(fd)

    monkeypatch.setattr(gs, "_fsync", record)
    publish(store, SHA2, b"two")
    monkeypatch.undo()
    base = Path(root) / "v1"
    for directory in (base, base / "gen", base / "gen" / SHA2[:2]):
        assert directory.stat().st_ino in log


def test_conflicting_bytes_are_refused_and_identical_bytes_are_a_no_op(root):
    store = make_store(root)
    publish(store, SHA, b"original")
    before = entry_path(root, SHA).stat()
    publish(store, SHA, b"original")
    after = entry_path(root, SHA).stat()
    assert (before.st_ino, before.st_mtime_ns) == (after.st_ino, after.st_mtime_ns)
    with pytest.raises(gs.ConflictError):
        publish(store, SHA, b"different")
    assert store.peek(SHA) == gs.Hit(b"original")
    assert _bucket_listing(root) == [f"{SHA}.entry"]


def test_short_writes_are_looped(root, monkeypatch):
    store = make_store(root)
    calls = []
    real_write = os.write

    def dribble(fd, data):
        calls.append(len(data))
        return real_write(fd, bytes(data)[:7])

    monkeypatch.setattr(gs, "_write", dribble)
    payload = b"short writes are looped\n" * 40
    publish(store, SHA, payload)
    monkeypatch.undo()
    assert len(calls) > len(payload) // 7
    assert store.peek(SHA) == gs.Hit(payload)


def test_oversize_payloads_and_non_bytes_are_rejected_before_writing(root):
    store = make_store(root, max_entry_bytes=1024)
    with store.claim(SHA, timeout=30) as claim:
        with pytest.raises(gs.EntryTooLargeError):
            claim.publish(b"z" * 5000)
        with pytest.raises(gs.StoreError):
            claim.publish("text")  # type: ignore[arg-type]
    assert not entry_path(root, SHA).exists()


def test_publishing_over_a_corrupt_entry_quarantines_it_first(root):
    data = _seed(root)
    store = make_store(root)
    with store.claim(SHA, timeout=30) as claim:
        rewrite(entry_path(root, SHA), _flip(data))
        claim.publish(b"fresh")
    assert len(quarantine_files(root)) == 2
    assert store.peek(SHA) == gs.Hit(b"fresh")


def test_only_claims_can_publish(root):
    store = make_store(root)
    assert not hasattr(store, "publish")


# ---------------------------------------------------------------------------
# Claims, locks, concurrency
# ---------------------------------------------------------------------------


def _claim_worker(root, sha, barrier, evidence_dir):
    store = make_store(root)
    barrier.wait(30)
    generated = False
    with store.claim(sha, timeout=60) as claim:
        if isinstance(claim.outcome, gs.Hit):
            seen = claim.outcome.payload
        else:
            generated = True
            time.sleep(0.3)
            seen = f"made-by-{os.getpid()}".encode()
            claim.publish(seen)
    final = store.peek(sha)
    assert isinstance(final, gs.Hit)
    (Path(evidence_dir) / str(os.getpid())).write_text(
        json.dumps({"generated": generated, "seen": seen.hex(), "final": final.payload.hex()})
    )


def test_four_processes_claiming_one_sha_generate_exactly_once(root, tmp_path):
    evidence = tmp_path / "evidence"
    evidence.mkdir()
    publish(make_store(root), SHA2, b"warm the directory tree")
    barrier = CTX.Barrier(4)
    procs = [CTX.Process(target=_claim_worker, args=(root, SHA, barrier, str(evidence))) for _ in range(4)]
    for proc in procs:
        proc.start()
    for proc in procs:
        proc.join(90)
        assert proc.exitcode == 0
    records = [json.loads(p.read_text()) for p in evidence.iterdir()]
    assert len(records) == 4
    assert sum(r["generated"] for r in records) == 1
    assert len({r["seen"] for r in records}) == 1
    assert len({r["final"] for r in records}) == 1
    assert records[0]["seen"] == records[0]["final"]


def test_sequential_callers_sharing_a_sha_generate_once(root):
    generations = []

    def get_or_generate(store):
        with store.claim(SHA, timeout=30) as claim:
            if isinstance(claim.outcome, gs.Hit):
                return claim.outcome.payload
            generations.append(1)
            claim.publish(b"generated")
            return b"generated"

    assert get_or_generate(make_store(root)) == b"generated"
    assert get_or_generate(make_store(root)) == b"generated"
    assert len(generations) == 1


def test_nested_claims_raise_and_the_slot_is_released_afterwards(root):
    store = make_store(root)
    with store.claim(SHA, timeout=30):
        with pytest.raises(gs.NestedClaimError):
            store.claim(SHA2, timeout=1)
        with pytest.raises(gs.NestedClaimError):
            store.claim(SHA, timeout=1)
        with pytest.raises(gs.NestedClaimError):
            make_store(root).claim(SHA3, timeout=1)
    with pytest.raises(RuntimeError):
        with store.claim(SHA2, timeout=30):
            raise RuntimeError("generation failed")
    with store.claim(SHA3, timeout=30) as claim:
        assert isinstance(claim.outcome, gs.Miss)


def test_peek_is_allowed_while_a_claim_is_held(root):
    store = make_store(root)
    publish(store, SHA2, b"other")
    with store.claim(SHA, timeout=30):
        assert store.peek(SHA2) == gs.Hit(b"other")
        assert store.peek(SHA) == gs.Miss()


def test_a_claim_is_unusable_after_exit_and_in_a_forked_child(root):
    store = make_store(root)
    with store.claim(SHA, timeout=30) as claim:
        read, write = os.pipe()
        pid = os.fork()
        if pid == 0:
            os.close(read)
            code = 3
            try:
                try:
                    claim.publish(b"from the child")
                except gs.ClaimError:
                    code = 0
                try:
                    claim.check()
                except gs.ClaimError:
                    code += 10
                # a child must build its own store and its own claim
                child_store = make_store(root)
                with child_store.claim(SHA2, timeout=10) as own:
                    own.publish(b"child's own")
                    code += 100
            finally:
                os.write(write, str(code).encode())
                os._exit(0)
        os.close(write)
        report = os.read(read, 16).decode()
        os.close(read)
        os.waitpid(pid, 0)
        assert report == "110", report
    with pytest.raises(gs.ClaimError):
        claim.publish(b"after exit")
    with pytest.raises(gs.ClaimError):
        claim.check()
    assert store.peek(SHA) == gs.Miss()
    assert store.peek(SHA2) == gs.Hit(b"child's own")


def test_a_forked_child_does_not_keep_the_parents_lock_alive(root):
    store = make_store(root)
    ready = CTX.Event()
    release = CTX.Event()

    def child_main():
        ready.set()
        release.wait(30)

    with store.claim(SHA, timeout=30):
        child = CTX.Process(target=child_main)
        child.start()
        assert ready.wait(30)
    # the parent released; the still-running forked child must not hold the lock
    try:
        with make_store(root).claim(SHA, timeout=5) as claim:
            assert isinstance(claim.outcome, gs.Miss)
    finally:
        release.set()
        child.join(30)


def _hold_claim(root, sha, started, release, ino_file=None):
    store = make_store(root)
    with store.claim(sha, timeout=60):
        started.set()
        release.wait(60)


def test_claim_timeout_raises_and_leaves_no_claim_active(root):
    started, release = CTX.Event(), CTX.Event()
    holder = CTX.Process(target=_hold_claim, args=(root, SHA, started, release))
    holder.start()
    try:
        assert started.wait(30)
        store = make_store(root)
        began = time.monotonic()
        with pytest.raises(gs.ClaimTimeoutError):
            store.claim(SHA, timeout=0.4)
        assert 0.3 < time.monotonic() - began < 10
        with store.claim(SHA2, timeout=10) as claim:
            assert isinstance(claim.outcome, gs.Miss)
    finally:
        release.set()
        holder.join(30)


def test_hostile_lock_files_are_refused(root):
    store = make_store(root)
    publish(store, SHA2, b"x")
    locks = Path(root) / "v1" / "locks"
    (locks / f"{SHA}.lock").symlink_to(Path(root).parent / "elsewhere")
    with pytest.raises(gs.UntrustedStoreError):
        store.claim(SHA, timeout=5)
    assert not (Path(root).parent / "elsewhere").exists()
    (locks / f"{SHA}.lock").unlink()
    os.mkfifo(locks / f"{SHA}.lock", 0o600)
    with pytest.raises(gs.UntrustedStoreError):
        store.claim(SHA, timeout=5)
    (locks / f"{SHA}.lock").unlink()
    first = locks / f"{SHA}.lock"
    fd = os.open(first, os.O_CREAT | os.O_WRONLY, 0o600)
    os.close(fd)
    os.link(first, Path(root).parent / "lockalias")
    with pytest.raises(gs.UntrustedStoreError):
        store.claim(SHA, timeout=5)


def test_lock_files_are_never_removed_or_replaced(root):
    store = make_store(root)
    publish(store, SHA, b"x")
    lock = Path(root) / "v1" / "locks" / f"{SHA}.lock"
    inode = lock.stat().st_ino
    for _ in range(3):
        with store.claim(SHA, timeout=30):
            pass
    rewrite(entry_path(root, SHA), b"corrupt")
    with store.claim(SHA, timeout=30) as claim:
        assert isinstance(claim.outcome, gs.Refused)
    assert lock.stat().st_ino == inode and lock.stat().st_size == 0


def _timed_claim(root, sha, name, events_dir, release, waiting):
    store = make_store(root)
    waiting.set()
    with store.claim(sha, timeout=60) as claim:
        (Path(events_dir) / f"{name}.acquired").write_text(f"{time.monotonic()}|{type(claim.outcome).__name__}")
        release.wait(60)
        (Path(events_dir) / f"{name}.released").write_text(f"{time.monotonic()}")


def test_quarantine_does_not_disturb_the_lock_so_a_third_contender_still_queues(root, tmp_path):
    publish(make_store(root), SHA, b"x")
    rewrite(entry_path(root, SHA), b"corrupt entry")
    lock = Path(root) / "v1" / "locks" / f"{SHA}.lock"
    inode_before = lock.stat().st_ino
    events = tmp_path / "events"
    events.mkdir()
    h_release, w_release, x_release = CTX.Event(), CTX.Event(), CTX.Event()
    h_wait, w_wait, x_wait = CTX.Event(), CTX.Event(), CTX.Event()
    procs = {}
    try:
        procs["H"] = CTX.Process(target=_timed_claim, args=(root, SHA, "H", str(events), h_release, h_wait))
        procs["H"].start()
        deadline = time.monotonic() + 30
        while not (events / "H.acquired").exists():
            assert time.monotonic() < deadline
            time.sleep(0.02)
        assert (events / "H.acquired").read_text().endswith("Refused")
        assert len(quarantine_files(root)) == 2, "H quarantined while holding the lock"

        procs["W"] = CTX.Process(target=_timed_claim, args=(root, SHA, "W", str(events), w_release, w_wait))
        procs["W"].start()
        assert w_wait.wait(30)
        time.sleep(0.4)
        procs["X"] = CTX.Process(target=_timed_claim, args=(root, SHA, "X", str(events), x_release, x_wait))
        procs["X"].start()
        assert x_wait.wait(30)
        time.sleep(0.4)
        assert not (events / "W.acquired").exists() and not (events / "X.acquired").exists()

        h_release.set()
        deadline = time.monotonic() + 30
        while not ((events / "W.acquired").exists() or (events / "X.acquired").exists()):
            assert time.monotonic() < deadline
            time.sleep(0.02)
        time.sleep(0.6)
        acquired = [n for n in ("W", "X") if (events / f"{n}.acquired").exists()]
        assert len(acquired) == 1, f"both contenders hold the lock: {acquired}"
        first = acquired[0]
        other = "X" if first == "W" else "W"
        assert (events / f"{first}.acquired").read_text().endswith("Miss"), "the quarantine is not repeated"
        assert lock.stat().st_ino == inode_before, "the lock file was replaced"
        assert not (events / f"{other}.acquired").exists()
        {"W": w_release, "X": x_release}[first].set()
        deadline = time.monotonic() + 30
        while not (events / f"{other}.acquired").exists():
            assert time.monotonic() < deadline
            time.sleep(0.02)
        assert float((events / f"{first}.released").read_text()) <= float(
            (events / f"{other}.acquired").read_text().split("|")[0]
        )
        assert lock.stat().st_ino == inode_before
    finally:
        for event in (h_release, w_release, x_release):
            event.set()
        for proc in procs.values():
            proc.join(30)
            if proc.is_alive():
                proc.kill()


# ---------------------------------------------------------------------------
# Process death
# ---------------------------------------------------------------------------


def _generate_then_hang(root, sha, ready, hang_before_rename=False):
    store = make_store(root)
    with store.claim(sha, timeout=60) as claim:
        if hang_before_rename:

            def hang(*args, **kwargs):
                ready.set()
                time.sleep(600)

            gs._replace = hang
            claim.publish(b"never lands")
        ready.set()
        time.sleep(600)


def test_sigkill_of_the_claiming_child_releases_the_lock_and_the_next_claimant_generates(root):
    publish(make_store(root), SHA2, b"warm")
    ready = CTX.Event()
    child = CTX.Process(target=_generate_then_hang, args=(root, SHA, ready))
    child.start()
    try:
        assert ready.wait(30)
        with pytest.raises(gs.ClaimTimeoutError):
            make_store(root).claim(SHA, timeout=0.3)
        os.kill(child.pid, signal.SIGKILL)
        child.join(30)
    finally:
        if child.is_alive():
            child.kill()
    assert not entry_path(root, SHA).exists()
    store = make_store(root)
    with store.claim(SHA, timeout=15) as claim:
        assert isinstance(claim.outcome, gs.Miss)
        claim.publish(b"second attempt")
    assert store.peek(SHA) == gs.Hit(b"second attempt")


def test_sigkill_between_temp_write_and_rename_leaves_no_entry_and_the_next_claimant_publishes(root):
    publish(make_store(root), SHA2, b"warm")
    ready = CTX.Event()
    child = CTX.Process(target=_generate_then_hang, args=(root, SHA, ready, True))
    child.start()
    try:
        assert ready.wait(30)
        os.kill(child.pid, signal.SIGKILL)
        child.join(30)
    finally:
        if child.is_alive():
            child.kill()
    assert not entry_path(root, SHA).exists()
    assert make_store(root).peek(SHA) == gs.Miss()
    store = make_store(root)
    with store.claim(SHA, timeout=15) as claim:
        assert isinstance(claim.outcome, gs.Miss)
        claim.publish(b"recovered")
    assert store.peek(SHA) == gs.Hit(b"recovered")


PARENT_SCRIPT = textwrap.dedent(
    """
    import subprocess, sys, time
    child = subprocess.Popen([sys.executable, "-c", sys.argv[1]] + sys.argv[2:], start_new_session=True)
    open(sys.argv[5] + ".parent_ready", "w").write(str(child.pid))
    time.sleep(600)
    """
)

CHILD_SCRIPT = textwrap.dedent(
    """
    import os, sys, time
    sys.path.insert(0, {repo!r})
    from rmgpy.kmc import generation_store as gs
    root, sha, producer, marker = sys.argv[1:5]
    if os.stat("/").st_uid == 65534 and os.getuid() != 65534:  # unmapped root in a user namespace
        class Shim:
            def __init__(self, real):
                self._real = real
                self.st_uid = 0 if real.st_uid == 65534 else real.st_uid
            def __getattr__(self, name):
                return getattr(self._real, name)
        real_fstat = gs._fstat
        gs._fstat = lambda fd: Shim(real_fstat(fd))
    store = gs.GenerationStore(root, producer_identity=producer)
    with store.claim(sha, timeout=60) as claim:
        open(marker + ".child_ready", "w").write(str(os.getpid()))
        while not os.path.exists(marker + ".go"):
            time.sleep(0.05)
        claim.publish(b"published by the surviving child")
    """
)


def test_sigkill_of_the_parent_does_not_release_an_independent_childs_claim(root, tmp_path):
    publish(make_store(root), SHA2, b"warm")
    marker = str(tmp_path / "marker")
    child_script = CHILD_SCRIPT.format(repo=REPO_ROOT)
    parent = subprocess.Popen(
        [sys.executable, "-c", PARENT_SCRIPT, child_script, root, SHA, PRODUCER, marker], start_new_session=True
    )
    child_pid = None
    try:
        deadline = time.monotonic() + 60
        while not (os.path.exists(marker + ".child_ready") and os.path.exists(marker + ".parent_ready")):
            assert time.monotonic() < deadline, "child never claimed"
            time.sleep(0.05)
        child_pid = int(open(marker + ".parent_ready").read())
        os.kill(parent.pid, signal.SIGKILL)
        parent.wait(30)
        store = make_store(root)
        with pytest.raises(gs.ClaimTimeoutError):
            store.claim(SHA, timeout=0.5)
        open(marker + ".go", "w").close()
        with store.claim(SHA, timeout=60) as claim:
            assert claim.outcome == gs.Hit(b"published by the surviving child")
    finally:
        if parent.poll() is None:
            parent.kill()
        if child_pid:
            try:
                os.kill(child_pid, signal.SIGKILL)
            except ProcessLookupError:
                pass


# ---------------------------------------------------------------------------
# Imports are out of scope
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", ["import_entry", "import_entries", "migrate"])
def test_import_apis_refuse_with_a_clear_error(root, name):
    store = make_store(root)
    with pytest.raises(gs.ImportsNotSupportedError, match="imports not supported"):
        getattr(store, name)(root, SHA)
    assert not os.path.exists(root)
