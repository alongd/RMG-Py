"""Stub of the generation store (red phase): every operation raises NotImplementedError."""

from __future__ import annotations

import fcntl
import os
from dataclasses import dataclass
from typing import Optional

ALLOWED_FILESYSTEMS = ("ext4", "xfs", "btrfs")
HEADER_LENGTH = 72
DEFAULT_MAX_ENTRY_BYTES = 64 << 20


class StoreError(Exception):
    pass


class UntrustedStoreError(StoreError):
    pass


class PublicationUncertainError(StoreError):
    pass


class ConflictError(StoreError):
    pass


class NestedClaimError(StoreError):
    pass


class ClaimError(StoreError):
    pass


class ClaimTimeoutError(StoreError, TimeoutError):
    pass


class ImportsNotSupportedError(StoreError):
    pass


class EntryTooLargeError(StoreError, ValueError):
    pass


class InvalidIdentityError(StoreError, ValueError):
    pass


@dataclass(frozen=True)
class Hit:
    payload: bytes


@dataclass(frozen=True)
class Miss:
    pass


@dataclass(frozen=True)
class Refused:
    reason: str


# Indirections that the tests replace to observe or fail one operation.
_fstat = os.fstat
_fsync = os.fsync
_replace = os.replace
_write = os.write
_flock = fcntl.flock


def _before_quarantine(claim) -> None:
    pass


def _read_mountinfo() -> str:
    raise NotImplementedError


def parse_mountinfo(text: str) -> list:
    raise NotImplementedError


def filesystem_refusal(root: str, mounts: list) -> Optional[str]:
    raise NotImplementedError


def encode_entry(key: bytes, request_sha: str, envelope_bytes: bytes) -> bytes:
    raise NotImplementedError


class Claim:
    def __init__(self, *args, **kwargs):
        raise NotImplementedError

    def check(self):
        raise NotImplementedError

    def _install(self, entry_bytes):
        raise NotImplementedError


class GenerationStore:
    def __init__(self, root, *, producer_identity, max_entry_bytes=DEFAULT_MAX_ENTRY_BYTES):
        raise NotImplementedError

    def peek(self, request_sha):
        raise NotImplementedError

    def claim(self, request_sha, timeout=None):
        raise NotImplementedError

    def import_entry(self, *args, **kwargs):
        raise NotImplementedError

    def import_entries(self, *args, **kwargs):
        raise NotImplementedError

    def migrate(self, *args, **kwargs):
        raise NotImplementedError
