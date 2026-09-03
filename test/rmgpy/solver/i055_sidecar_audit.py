#!/usr/bin/env python3
"""Audit a polymer_pools.json payload with CKMG's own liveness auditor.

Runs under ck_env. CKMG is imported READ-ONLY -- nothing in /home/alon/Code/CKMG
is written, and no CKMG policy, whitelist or tolerance is touched. Each argv
entry is 'label=path'.

Requires, for every CONFIGURED pool whose unzip channel carries A > 0:
  * pool.mass_loss_channels['unzip'] is True  -- CKMG still counts the channel
    as removing pool mass (its rule: unzip counts iff A > 0);
  * pool.structurally_dead is False.
"""

import json
import sys

from ckmg.utils.sidecar import audit_sidecar_liveness


def audit(label, path):
    with open(path) as f:
        raw = json.load(f)
    rec = audit_sidecar_liveness(raw)
    print(f"  [{label}] {path}")
    print(f"    schema_version        = {rec.schema_version!r}")
    print(f"    auditable             = {rec.auditable!r}")
    print(f"    stale_topology        = {rec.stale_topology!r}")
    print(f"    dead_configured_pools = {rec.dead_configured_pools!r}")

    if not rec.auditable:
        print(f"    FAIL: not auditable ({rec.not_auditable_reason!r})")
        return False

    # Which pools does the RAW payload declare with unzip A > 0? Those are the
    # ones this fix must keep alive; a pool without the channel is irrelevant.
    with_unzip = []
    for block in raw.get("pools", []):
        try:
            a = float(((block.get("channels") or {}).get("unzip") or {})
                      .get("A", 0.0) or 0.0)
        except (TypeError, ValueError):
            a = 0.0
        if a > 0.0:
            with_unzip.append((block.get("label"), a))

    if not with_unzip:
        print("    FAIL: payload declares no pool with unzip A > 0, so this "
              "check would be vacuous")
        return False

    ok = True
    for lab, a in with_unzip:
        pool = rec.pools.get(lab)
        if pool is None:
            print(f"    FAIL: pool {lab!r} (unzip A={a:g}) absent from the audit")
            ok = False
            continue
        chans = dict(pool.mass_loss_channels)
        print(f"    pool {lab!r}: unzip A={a:g}  configured={pool.configured!r}")
        print(f"      mass_loss_channels      = {chans!r}")
        print(f"      live_mass_loss_channels = {pool.live_mass_loss_channels!r}")
        print(f"      reachable_terminal_mass_loss = "
              f"{pool.reachable_terminal_mass_loss!r}")
        print(f"      structurally_dead       = {pool.structurally_dead!r}")
        if chans.get("unzip") is not True:
            print(f"      FAIL: mass_loss_channels['unzip'] is "
                  f"{chans.get('unzip')!r}, not True")
            ok = False
        else:
            print("      OK: mass_loss_channels['unzip'] is True")
        if pool.structurally_dead is not False:
            print(f"      FAIL: structurally_dead is "
                  f"{pool.structurally_dead!r}, not False")
            ok = False
        else:
            print("      OK: structurally_dead is False")
    return ok


def main():
    all_ok = True
    for arg in sys.argv[1:]:
        label, _, path = arg.partition("=")
        all_ok &= audit(label, path)
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
