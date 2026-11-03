"""A node-local, authenticated, durable, locked store for generation envelopes.

The store keeps one file per request identity.  It handles envelope bytes and
returns the verified envelope payload; it never decodes reactions and it is not
wired into compilation.

Threat model
------------
Node-local cache for one user on one machine.

* TRUSTED: every process running as this uid.  Such a process can read the key
  and forge entries; authentication here proves only possession of the local
  key, never correct chemistry.
* DEFENDED AGAINST: other local users, crashes and process death, concurrent
  compilers, accidental corruption, hostile filesystem objects placed by other
  users (symlinks, FIFOs, hard links, wrong modes), and malformed bytes.
* OUT OF SCOPE: imports of foreign artifacts (every import or migration API
  refuses), shared or network filesystems (the root's filesystem type must be
  ext4, xfs or btrfs), and power-loss proofs.  The ``fsync`` ordering is the
  conventional one, but the tests prove process-interruption behaviour only.

Layout (every directory is uid-owned mode 0700, every file mode 0600)::

    <root>/v1/gen/<sha[0:2]>/<sha>.entry     authenticated entries
    <root>/v1/locks/<sha>.lock               per-request lock files
    <root>/v1/quarantine/                    refused entries and ``.reason`` files
    <root>/v1/key                            32 random bytes
    <root>/v1/init.lock                      serialises key creation

Entry format: ``b"KMCGS1 " + hex(hmac) + b"\\n" + envelope_bytes`` with
``hmac = HMAC-SHA256(key, b"kmc-generation-store/v1\\0" + request_sha + b"\\0" +
envelope_bytes)``.  Reads are bounded by the store cap.

Trust is bound to descriptors, not path strings.  The root is walked from ``/``
one component at a time with ``O_NOFOLLOW | O_DIRECTORY``; every later open,
create, rename and unlink goes through those descriptors.  An untrusted root
(or key) disables reads and writes and never quarantines.

Outcomes
--------
``peek(sha)`` and ``Claim.outcome`` are ``Hit(payload)``, ``Miss()`` or
``Refused(reason)``.  ``Refused`` is returned only for the specific checks in
this module (untrusted root or key, hostile entry objects, failed
authentication or envelope verification); every other exception propagates.
``peek`` never mutates anything.  Quarantine happens only inside a claim on a
trusted root, after re-reading and re-checking the entry under the lock, so a
reader never quarantines a valid replacement.

Claims
------
``store.claim(sha, timeout=None)`` takes an exclusive ``flock`` on a lock file
that is never renamed, truncated or unlinked.  A process holds at most ONE claim
at a time (a second claim, for any sha, raises ``NestedClaimError``), which
removes lock-order deadlocks.  A claim records ``os.getpid()`` and refuses use in
another pid and after exit.  There is no descriptor-passing API: a child that
generates must construct its own store and its own claim, so the parent must not
hold a claim across a spawn or fork that expects the child to proceed.  (Forked
children close the inherited lock descriptor without unlocking it.)  Lock files
are advisory locks on a local filesystem; a crashed holder releases on exit.

Known limits
------------
* Stale ``<sha>.entry.tmp.<pid>`` files left by a killed publisher are harmless
  but are not swept; only a same-pid leftover is replaced.
* The quarantine directory is unbounded.
* The identity comparison before a quarantine rename is a check, not an atomic
  compare-and-rename; the lock excludes every cooperating process, and a
  same-uid writer that ignores the lock is trusted by the threat model.

Prerequisites for reuse (6a-4)
------------------------------
* Generation mutates global state that a cache hit would skip: the atom-id
  counter (``rmgpy/molecule/molecule.py`` ~2891) and the missing-reverse counts
  (``rmgpy/data/kinetics/family.py`` ~2076).
* The database snapshot is verified, not materialized, so loading can race
  edits between verification and use.
* ``generation_identity.UNPINNED_ITEMS`` lists inputs the identities do not pin.
* Cache size limits (``max_entry_bytes``) must be reconciled with the codec and
  parity-driver decode limits.
* Attribution replay bypasses the generation hook (``compiler.py`` ~1048).
"""

from __future__ import annotations

import errno
import fcntl
import hashlib
import hmac
import os
import posixpath
import re
import secrets
import stat
import threading
import time
from dataclasses import dataclass
from typing import Any, Optional

from rmgpy.kmc import generation_identity as gi

ALLOWED_FILESYSTEMS = ("ext4", "xfs", "btrfs")
MAGIC = b"KMCGS1 "
HEADER_LENGTH = len(MAGIC) + 64 + 1
DEFAULT_MAX_ENTRY_BYTES = gi.DEFAULT_MAX_ENVELOPE_BYTES
_MAC_DOMAIN = b"kmc-generation-store/v1\0"
_HEX64 = re.compile(r"[0-9a-f]{64}")
_HEX64_BYTES = re.compile(rb"[0-9a-f]{64}")
_READ_FLAGS = os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK | os.O_NOCTTY | os.O_CLOEXEC
_DIR_FLAGS = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC


class StoreError(Exception):
    """Base class for store errors."""


class UntrustedStoreError(StoreError):
    """The root, the key or a lock file cannot be trusted; nothing was read or written."""


class PublicationUncertainError(StoreError):
    """The rename happened but a later step failed: the entry may be visible, durability is unknown."""


class ConflictError(StoreError):
    """A different payload is already published for this request."""


class NestedClaimError(StoreError):
    """This process already holds a claim."""


class ClaimError(StoreError):
    """A claim was used after it ended, or in a process other than the one that took it."""


class ClaimTimeoutError(StoreError, TimeoutError):
    """The lock was not acquired before the timeout."""


class ImportsNotSupportedError(StoreError):
    """Imports and migrations of foreign artifacts are not supported."""


class EntryTooLargeError(StoreError, ValueError):
    """The envelope exceeds the store cap."""


class InvalidIdentityError(StoreError, ValueError):
    """A request sha or producer identity is not 64 lowercase hex characters."""


@dataclass(frozen=True)
class Hit:
    payload: bytes


@dataclass(frozen=True)
class Miss:
    pass


@dataclass(frozen=True)
class Refused:
    reason: str


class _Refusal(Exception):
    """Internal: a specific, expected refusal with a human-readable reason."""


# Thin indirections so the tests can observe or fail exactly one operation.
_fstat = os.fstat
_fsync = os.fsync
_replace = os.replace
_write = os.write
_flock = fcntl.flock


def _before_quarantine(claim: "Claim") -> None:
    """Test seam between observing a refused entry and quarantining it."""


def _validate_identity(value: Any, name: str) -> str:
    if type(value) is not str or not _HEX64.fullmatch(value):
        raise InvalidIdentityError(f"{name} must be exactly 64 lowercase hex characters")
    return value


# ---------------------------------------------------------------------------
# Filesystem allowlist (/proc/self/mountinfo)
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Mount:
    mount_point: str
    fstype: str


_OCTAL = re.compile(rb"\\([0-7]{3})")


def _unescape(field: str) -> str:
    raw = _OCTAL.sub(lambda m: bytes([int(m.group(1), 8) & 0xFF]), field.encode("utf-8", "surrogateescape"))
    return os.fsdecode(raw)


def parse_mountinfo(text: str) -> list[Mount]:
    """Parse ``/proc/self/mountinfo`` (proc_pid_mountinfo(5)): optional fields, octal escapes."""
    mounts = []
    for line in text.splitlines():
        if not line.strip():
            continue
        fields = line.split(" ")
        try:
            separator = fields.index("-", 6)
        except ValueError:
            raise ValueError(f"malformed mountinfo line: {line!r}") from None
        if len(fields) < separator + 3:
            raise ValueError(f"malformed mountinfo line: {line!r}")
        mounts.append(Mount(_unescape(fields[4]), fields[separator + 1]))
    return mounts


def _read_mountinfo() -> str:
    with open("/proc/self/mountinfo", "rb") as handle:
        return handle.read().decode("utf-8", "surrogateescape")


def filesystem_refusal(root: str, mounts: list[Mount]) -> Optional[str]:
    """Return why ``root`` must be refused, or None if its filesystem is allowed."""

    def contains(point: str, path: str) -> bool:
        return point == "/" or path == point or path.startswith(point.rstrip("/") + "/")

    containing = [(len(m.mount_point), index, m) for index, m in enumerate(mounts) if contains(m.mount_point, root)]
    if not containing:
        return "no mount point contains the root; its filesystem type is unknown"
    best = max(containing, key=lambda item: item[:2])[2]
    if best.fstype not in ALLOWED_FILESYSTEMS:
        return (
            f"filesystem type {best.fstype!r} (mounted at {best.mount_point!r}) is not allowed; "
            f"allowed: {', '.join(ALLOWED_FILESYSTEMS)}"
        )
    prefix = root.rstrip("/") + "/"
    for mount in mounts:
        if mount.mount_point != root and mount.mount_point.startswith(prefix):
            return f"mount point {mount.mount_point!r} lies inside the store root"
    return None


# ---------------------------------------------------------------------------
# Descriptor-bound directory tree
# ---------------------------------------------------------------------------


def _identity(st: Any) -> tuple:
    return (st.st_dev, st.st_ino, st.st_size, st.st_mtime_ns, stat.S_IFMT(st.st_mode))


def _check_ancestor(st: Any, name: str) -> None:
    if st.st_uid not in (0, os.getuid()):
        raise _Refusal(f"ancestor {name!r} has owner uid {st.st_uid}; expected this user or root")
    if st.st_mode & (stat.S_IWGRP | stat.S_IWOTH) and not st.st_mode & stat.S_ISVTX:
        raise _Refusal(f"ancestor {name!r} is group- or other-writable and not sticky")


def _check_private_dir(st: Any, name: str) -> None:
    if st.st_uid != os.getuid():
        raise _Refusal(f"directory {name!r} has owner uid {st.st_uid}; expected uid {os.getuid()}")
    mode = stat.S_IMODE(st.st_mode)
    if mode != 0o700:
        raise _Refusal(f"directory {name!r} has mode {mode:04o}; expected 0700")


def _describe_non_directory(parent: int, name: str) -> str:
    try:
        info = os.stat(name, dir_fd=parent, follow_symlinks=False)
    except OSError:
        return f"{name!r} cannot be opened as a directory"
    if stat.S_ISLNK(info.st_mode):
        return f"{name!r} is a symbolic link"
    return f"{name!r} is not a directory"


def _open_dir_at(parent: int, name: str, create: bool) -> Optional[int]:
    """Open directory ``name`` under ``parent`` without following links; None if absent and not creating."""
    for _ in range(3):
        try:
            return os.open(name, _DIR_FLAGS, dir_fd=parent)
        except FileNotFoundError:
            if not create:
                return None
            try:
                os.mkdir(name, 0o700, dir_fd=parent)
            except FileExistsError:
                pass
        except OSError as error:
            if error.errno in (errno.ELOOP, errno.ENOTDIR):
                raise _Refusal(_describe_non_directory(parent, name)) from error
            raise
    raise _Refusal(f"{name!r} kept disappearing while it was being created")


def _read_up_to(fd: int, limit: int) -> bytes:
    chunks = []
    total = 0
    while total < limit:
        chunk = os.read(fd, min(1 << 20, limit - total))
        if not chunk:
            break
        chunks.append(chunk)
        total += len(chunk)
    return b"".join(chunks)


def _write_all(fd: int, data: bytes) -> None:
    view = memoryview(data)
    while len(view):
        written = _write(fd, view)
        if written <= 0:
            raise OSError(errno.EIO, "write made no progress")
        view = view[written:]


class _Tree:
    """The open, validated descriptors of one store."""

    def __init__(self) -> None:
        self.chain: list[int] = []  # "/" ... store root
        self.v1 = self.gen = self.locks = self.quarantine = -1

    @property
    def root(self) -> int:
        return self.chain[-1]

    def close(self) -> None:
        fds = list(self.chain) + [fd for fd in (self.v1, self.gen, self.locks, self.quarantine) if fd >= 0]
        self.chain = []
        self.v1 = self.gen = self.locks = self.quarantine = -1
        for fd in fds:
            try:
                os.close(fd)
            except OSError:
                pass


def _check_root_path(root: str) -> list[str]:
    if not root.startswith("/") or "\0" in root:
        raise _Refusal("the store root must be an absolute path")
    parts = root.split("/")[1:]
    if not parts or any(part in ("", ".", "..") for part in parts):
        raise _Refusal("the store root must be a normalised absolute path with no '.', '..' or empty components")
    return parts


def _open_tree(root: str, create: bool) -> Optional[_Tree]:
    """Walk to the root and open the internal directories; None if something is absent and not creating."""
    parts = _check_root_path(root)
    tree = _Tree()
    try:
        fd = os.open("/", os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC)
        tree.chain.append(fd)
        _check_ancestor(_fstat(fd), "/")
        for index, name in enumerate(parts):
            last = index == len(parts) - 1
            if last:
                try:
                    mounts = parse_mountinfo(_read_mountinfo())
                except (OSError, ValueError) as error:
                    raise _Refusal(f"cannot read the mount table: {error}") from error
                reason = filesystem_refusal(root, mounts)
                if reason is not None:
                    raise _Refusal(reason)
            child = _open_dir_at(tree.chain[-1], name, create and last)
            if child is None:
                tree.close()
                return None
            tree.chain.append(child)
            info = _fstat(child)
            if last:
                _check_private_dir(info, name)
            else:
                _check_ancestor(info, name)
        for attribute, name in (("v1", "v1"),):
            child = _open_dir_at(tree.root, name, create)
            if child is None:
                tree.close()
                return None
            setattr(tree, attribute, child)
            _check_private_dir(_fstat(child), name)
        for attribute, name in (("gen", "gen"), ("locks", "locks"), ("quarantine", "quarantine")):
            child = _open_dir_at(tree.v1, name, create)
            if child is None:
                tree.close()
                return None
            setattr(tree, attribute, child)
            _check_private_dir(_fstat(child), name)
        return tree
    except BaseException:
        tree.close()
        raise


def _open_bucket(tree: _Tree, request_sha: str, create: bool) -> Optional[int]:
    name = request_sha[:2]
    fd = _open_dir_at(tree.gen, name, create)
    if fd is None:
        return None
    try:
        _check_private_dir(_fstat(fd), f"gen/{name}")
    except BaseException:
        os.close(fd)
        raise
    return fd


def _check_private_file(st: Any, what: str, size: Optional[int] = None) -> None:
    if not stat.S_ISREG(st.st_mode):
        raise _Refusal(f"{what} is not a regular file")
    if st.st_uid != os.getuid():
        raise _Refusal(f"{what} has the wrong owner (uid {st.st_uid})")
    mode = stat.S_IMODE(st.st_mode)
    if mode != 0o600:
        raise _Refusal(f"{what} has mode {mode:04o}; expected 0600")
    if st.st_nlink != 1:
        raise _Refusal(f"{what} has {st.st_nlink} hard links")
    if size is not None and st.st_size != size:
        raise _Refusal(f"{what} has size {st.st_size}; expected {size}")


def _open_lock_file(dir_fd: int, name: str) -> int:
    try:
        fd = os.open(name, os.O_CREAT | os.O_RDWR | os.O_NOFOLLOW | os.O_CLOEXEC, 0o600, dir_fd=dir_fd)
    except OSError as error:
        if error.errno == errno.ELOOP:
            raise _Refusal(f"lock file {name!r} is a symbolic link") from error
        raise
    try:
        _check_private_file(_fstat(fd), f"lock file {name!r}")
    except BaseException:
        os.close(fd)
        raise
    return fd


# ---------------------------------------------------------------------------
# Key
# ---------------------------------------------------------------------------


def _read_key(v1: int) -> Optional[bytes]:
    try:
        fd = os.open("key", _READ_FLAGS, dir_fd=v1)
    except FileNotFoundError:
        return None
    except OSError as error:
        if error.errno == errno.ELOOP:
            raise _Refusal("key is a symbolic link") from error
        raise
    try:
        _check_private_file(_fstat(fd), "key", size=32)
        data = _read_up_to(fd, 33)
    except _Refusal as error:
        raise _Refusal(str(error)) from None
    finally:
        os.close(fd)
    if len(data) != 32:
        raise _Refusal(f"key has {len(data)} readable bytes; expected 32")
    return data


def _ensure_key(tree: _Tree) -> bytes:
    lock = _open_lock_file(tree.v1, "init.lock")
    try:
        _flock(lock, fcntl.LOCK_EX)
        try:
            os.unlink("key.tmp", dir_fd=tree.v1)
        except FileNotFoundError:
            pass
        key = _read_key(tree.v1)
        if key is not None:
            return key
        key = secrets.token_bytes(32)
        fd = os.open("key.tmp", os.O_CREAT | os.O_EXCL | os.O_WRONLY | os.O_NOFOLLOW | os.O_CLOEXEC, 0o600, dir_fd=tree.v1)
        try:
            _write_all(fd, key)
            _fsync(fd)
        finally:
            os.close(fd)
        os.rename("key.tmp", "key", src_dir_fd=tree.v1, dst_dir_fd=tree.v1)
        _fsync(tree.v1)
        for ancestor in reversed(tree.chain):
            _fsync(ancestor)
        return key
    finally:
        os.close(lock)


# ---------------------------------------------------------------------------
# Entries
# ---------------------------------------------------------------------------


def encode_entry(key: bytes, request_sha: str, envelope_bytes: bytes) -> bytes:
    """Frame ``envelope_bytes`` with the authenticated header."""
    mac = hmac.new(key, _MAC_DOMAIN + request_sha.encode("ascii") + b"\0" + envelope_bytes, hashlib.sha256)
    return MAGIC + mac.hexdigest().encode("ascii") + b"\n" + envelope_bytes


# ---------------------------------------------------------------------------
# Claims
# ---------------------------------------------------------------------------

_ACTIVE: Optional["Claim"] = None
_ACTIVE_LOCK = threading.Lock()


def _after_fork_in_child() -> None:
    """A forked child never inherits a claim: close the descriptor without unlocking."""
    global _ACTIVE, _ACTIVE_LOCK
    _ACTIVE_LOCK = threading.Lock()
    claim, _ACTIVE = _ACTIVE, None
    if claim is not None and claim._lock_fd >= 0:
        try:
            os.close(claim._lock_fd)
        except OSError:
            pass
        claim._lock_fd = -1
        claim._closed = True


os.register_at_fork(after_in_child=_after_fork_in_child)


def _acquire_lock(fd: int, timeout: Optional[float]) -> None:
    if timeout is None:
        while True:
            try:
                _flock(fd, fcntl.LOCK_EX)
                return
            except InterruptedError:
                continue
    end = time.monotonic() + timeout
    delay = 0.005
    while True:
        try:
            _flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
            return
        except BlockingIOError:
            pass
        except InterruptedError:
            continue
        remaining = end - time.monotonic()
        if remaining <= 0:
            raise ClaimTimeoutError(f"the lock was not acquired within {timeout} s")
        time.sleep(min(delay, remaining))
        delay = min(delay * 2, 0.05)


class Claim:
    """Exclusive ownership of one request sha.  Use as a context manager; created by ``store.claim``."""

    def __init__(self, store: "GenerationStore", request_sha: str) -> None:
        self._store = store
        self._sha = request_sha
        self._pid = os.getpid()
        self._lock_fd = -1
        self._closed = False
        self._tree: Optional[_Tree] = None
        self._key = b""
        self.outcome: Any = None
        self.quarantined: list[str] = []

    # -- lifecycle ----------------------------------------------------------

    @property
    def request_sha(self) -> str:
        return self._sha

    def _live(self) -> None:
        if self._pid != os.getpid():
            raise ClaimError("a claim cannot be used in a different process; build your own store and claim")
        if self._closed:
            raise ClaimError("the claim has ended")

    def release(self) -> None:
        global _ACTIVE
        if self._closed:
            return
        self._closed = True
        try:
            if self._pid == os.getpid() and self._lock_fd >= 0:
                try:
                    _flock(self._lock_fd, fcntl.LOCK_UN)
                except OSError:
                    pass
                os.close(self._lock_fd)
        finally:
            self._lock_fd = -1
            with _ACTIVE_LOCK:
                if _ACTIVE is self:
                    _ACTIVE = None

    def __enter__(self) -> "Claim":
        self._live()
        return self

    def __exit__(self, *exc_info: Any) -> None:
        self.release()

    # -- reading, quarantine -----------------------------------------------

    def check(self) -> Any:
        """Re-read the entry under the lock; a refused entry is quarantined first (re-checked)."""
        self._live()
        for _ in range(8):
            bucket = self._bucket(create=False)
            if bucket is None:
                return Miss()
            try:
                outcome, observed = self._store._read_entry(bucket, self._sha, self._key)
                if not isinstance(outcome, Refused):
                    return outcome
                _before_quarantine(self)
                if self._quarantine(bucket, observed, outcome.reason):
                    return outcome
            finally:
                os.close(bucket)
        raise StoreError("the entry kept changing while the claim was held")

    def _bucket(self, create: bool) -> Optional[int]:
        try:
            return _open_bucket(self._tree, self._sha, create)
        except _Refusal as refusal:
            raise UntrustedStoreError(str(refusal)) from refusal

    def _quarantine(self, bucket: int, observed: tuple, reason: str) -> bool:
        name = f"{self._sha}.entry"
        try:
            current = os.stat(name, dir_fd=bucket, follow_symlinks=False)
        except FileNotFoundError:
            return False
        if _identity(current) != observed:
            return False
        target = f"{self._sha}.{time.time_ns()}.{os.getpid()}"
        quarantine = self._tree.quarantine
        os.rename(name, target, src_dir_fd=bucket, dst_dir_fd=quarantine)
        fd = os.open(target + ".reason", os.O_CREAT | os.O_EXCL | os.O_WRONLY | os.O_NOFOLLOW | os.O_CLOEXEC, 0o600, dir_fd=quarantine)
        try:
            _write_all(fd, (reason + "\n").encode("utf-8", "replace"))
            _fsync(fd)
        finally:
            os.close(fd)
        _fsync(bucket)
        _fsync(quarantine)
        self.quarantined.append(target)
        return True

    # -- publication ----------------------------------------------------------

    def publish(self, payload: bytes) -> None:
        """Publish ``payload`` for this request (inside the claim only)."""
        self._live()
        if not isinstance(payload, (bytes, bytearray)):
            raise StoreError("payload must be bytes")
        payload = bytes(payload)
        envelope_bytes = gi.envelope_to_bytes(gi.build_envelope(self._sha, self._store._producer, payload))
        if len(envelope_bytes) > self._store._cap:
            raise EntryTooLargeError(f"envelope is {len(envelope_bytes)} bytes; the cap is {self._store._cap}")
        current = self.check()
        if isinstance(current, Hit):
            if current.payload == payload:
                return
            raise ConflictError(f"a different payload is already published for {self._sha}")
        self._install(encode_entry(self._key, self._sha, envelope_bytes))

    def _install(self, entry_bytes: bytes) -> None:
        """Write the entry atomically: temp, fsync file, rename, fsync bucket."""
        self._live()
        tree = self._tree
        bucket = self._bucket(create=True)
        try:
            for directory in (tree.v1, tree.gen, bucket):
                _fsync(directory)
            temp = f"{self._sha}.entry.tmp.{os.getpid()}"
            final = f"{self._sha}.entry"
            flags = os.O_CREAT | os.O_EXCL | os.O_WRONLY | os.O_NOFOLLOW | os.O_CLOEXEC
            try:
                fd = os.open(temp, flags, 0o600, dir_fd=bucket)
            except FileExistsError:
                os.unlink(temp, dir_fd=bucket)  # a leftover of a dead process that had this pid
                fd = os.open(temp, flags, 0o600, dir_fd=bucket)
            renamed = False
            try:
                try:
                    _write_all(fd, entry_bytes)
                    _fsync(fd)
                finally:
                    os.close(fd)
                _replace(temp, final, src_dir_fd=bucket, dst_dir_fd=bucket)
                renamed = True
            except BaseException:
                if not renamed:
                    try:
                        os.unlink(temp, dir_fd=bucket)
                    except OSError:
                        pass
                raise
            try:
                _fsync(bucket)
            except Exception as error:
                raise PublicationUncertainError(
                    f"the entry for {self._sha} was renamed into place but the directory fsync failed: {error}"
                ) from error
        finally:
            os.close(bucket)


# ---------------------------------------------------------------------------
# The store
# ---------------------------------------------------------------------------


class GenerationStore:
    """Handle on one store root, for one expected producer identity.  Nothing touches the filesystem until used."""

    def __init__(self, root: Any, *, producer_identity: str, max_entry_bytes: int = DEFAULT_MAX_ENTRY_BYTES) -> None:
        self._producer = _validate_identity(producer_identity, "producer_identity")
        if type(max_entry_bytes) is not int or max_entry_bytes <= 0:
            raise ValueError("max_entry_bytes must be a positive integer")
        self._root = os.fspath(root)
        self._cap = max_entry_bytes
        self._tree: Optional[_Tree] = None
        self._key: Optional[bytes] = None

    def close(self) -> None:
        if self._tree is not None:
            self._tree.close()
            self._tree = None

    def __enter__(self) -> "GenerationStore":
        return self

    def __exit__(self, *exc_info: Any) -> None:
        self.close()

    # -- internal -------------------------------------------------------------

    def _open(self, create: bool) -> Optional[_Tree]:
        if self._tree is None:
            self._tree = _open_tree(self._root, create)
        return self._tree

    def _read_entry(self, bucket: int, request_sha: str, key: bytes) -> tuple[Any, Optional[tuple]]:
        name = f"{request_sha}.entry"
        try:
            fd = os.open(name, _READ_FLAGS, dir_fd=bucket)
        except FileNotFoundError:
            return Miss(), None
        except OSError as error:
            if error.errno in (errno.ELOOP, errno.ENXIO):
                try:
                    observed = _identity(os.stat(name, dir_fd=bucket, follow_symlinks=False))
                except FileNotFoundError:
                    return Miss(), None
                what = "a symbolic link" if error.errno == errno.ELOOP else "not a regular file"
                return Refused(f"entry is {what}"), observed
            raise
        try:
            info = _fstat(fd)
            observed = _identity(info)
            try:
                _check_private_file(info, "entry")
            except _Refusal as refusal:
                return Refused(str(refusal)), observed
            limit = HEADER_LENGTH + self._cap
            if info.st_size > limit:
                return Refused("entry exceeds the size cap"), observed
            data = _read_up_to(fd, limit + 1)
        finally:
            os.close(fd)
        if len(data) > limit:
            return Refused("entry exceeds the size cap"), observed
        return self._authenticate(data, request_sha, key), observed

    def _authenticate(self, data: bytes, request_sha: str, key: bytes) -> Any:
        if len(data) < HEADER_LENGTH:
            return Refused("entry is truncated")
        if data[: len(MAGIC)] != MAGIC or data[HEADER_LENGTH - 1 : HEADER_LENGTH] != b"\n":
            return Refused("entry header is malformed")
        presented = data[len(MAGIC) : HEADER_LENGTH - 1]
        if not _HEX64_BYTES.fullmatch(presented):
            return Refused("entry header is malformed")
        body = data[HEADER_LENGTH:]
        expected = hmac.new(key, _MAC_DOMAIN + request_sha.encode("ascii") + b"\0" + body, hashlib.sha256)
        if not hmac.compare_digest(presented, expected.hexdigest().encode("ascii")):
            return Refused("entry authentication failed")
        try:
            envelope = gi.envelope_from_bytes(body, max_bytes=self._cap)
            payload = gi.verify_envelope(envelope, request_sha=request_sha, producer_identity=self._producer)
        except gi.EnvelopeError as error:
            return Refused(f"envelope rejected: {error}")
        return Hit(payload)

    # -- public ---------------------------------------------------------------

    def peek(self, request_sha: str) -> Any:
        """Unlocked read.  May return ``Refused``; never mutates anything."""
        _validate_identity(request_sha, "request_sha")
        try:
            tree = self._open(create=False)
            if tree is None:
                return Miss()
            if self._key is None:
                self._key = _read_key(tree.v1)
                if self._key is None:
                    return Miss()
            bucket = _open_bucket(tree, request_sha, create=False)
            if bucket is None:
                return Miss()
            try:
                outcome, _ = self._read_entry(bucket, request_sha, self._key)
            finally:
                os.close(bucket)
            return outcome
        except _Refusal as refusal:
            return Refused(str(refusal))

    def claim(self, request_sha: str, timeout: Optional[float] = None) -> Claim:
        """Take the exclusive claim on ``request_sha`` and check its entry (quarantining a refused one)."""
        global _ACTIVE
        _validate_identity(request_sha, "request_sha")
        if timeout is not None and (isinstance(timeout, bool) or not isinstance(timeout, (int, float)) or timeout < 0):
            raise ValueError("timeout must be None or a non-negative number of seconds")
        claim = Claim(self, request_sha)
        with _ACTIVE_LOCK:
            if _ACTIVE is not None and _ACTIVE._pid == os.getpid() and not _ACTIVE._closed:
                raise NestedClaimError("this process already holds a claim; a process may hold at most one")
            _ACTIVE = claim
        try:
            try:
                tree = self._open(create=True)
                if self._key is None:
                    self._key = _ensure_key(tree)
                claim._tree = tree
                claim._key = self._key
                claim._lock_fd = _open_lock_file(tree.locks, f"{request_sha}.lock")
            except _Refusal as refusal:
                raise UntrustedStoreError(str(refusal)) from refusal
            _acquire_lock(claim._lock_fd, timeout)
            claim.outcome = claim.check()
            return claim
        except BaseException:
            claim.release()
            raise

    def _refuse_imports(self, *args: Any, **kwargs: Any) -> None:
        raise ImportsNotSupportedError(
            "imports not supported: the store serves only entries this user wrote for exactly this request"
        )

    import_entry = _refuse_imports
    import_entries = _refuse_imports
    migrate = _refuse_imports
