#!/usr/bin/env python3
"""
Round 87 probe -- reproduce the four findings before any of them is fixed.

Run against the engine worktree with nothing else on the path::

    python docs/i221-saturation/probes/round87_provenance_probe.py

Every check below has a control that fails when the harness is wrong, because a
probe that cannot fail in the direction the defect lies reports PASS and proves
nothing. `REPRODUCED` means the defect is present as reported; `NOT REPRODUCED`
means the review's claim did not survive contact with the code and the claim is
what needs revisiting.
"""

import inspect
import math
import os
import sys
import tempfile
import traceback

from rmgpy.data.base import saturate_for_estimation
from rmgpy.data.kinetics.library import LibraryReaction
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.kinetics.quarantine import (
    QUARANTINE_FILENAME,
    check_quarantine,
    load_family_quarantine,
)
from rmgpy.exceptions import (
    AtomTypeError,
    DatabaseError,
    InvalidAdjacencyListError,
    QuarantinedKineticsError,
    SaturatedStructureError,
)
from rmgpy.kinetics import Marcus
from rmgpy.molecule import Molecule
from rmgpy.species import Species

FAMILY = "Fake_Quarantined_Family"

MANIFEST = """
name = "Fake_Quarantined_Family/quarantine"
state = "QUARANTINED FOR TESTING"
appliesToKineticsClass = "Marcus"
reason = "a reason that must reach the error message"
"""

results = []
controls = []


def record(name, reproduced, detail):
    results.append((name, reproduced, detail))
    flag = "REPRODUCED" if reproduced else "NOT REPRODUCED"
    print("\n[{0}] {1}\n    {2}".format(flag, name, detail.replace("\n", "\n    ")))


def control(name, ok, detail):
    """A control must hold BOTH before and after the fix, or the probe proves nothing."""
    controls.append((name, ok, detail))
    print("\n[{0}] control: {1}\n    {2}".format(
        "HOLDS" if ok else "BROKEN", name, detail.replace("\n", "\n    ")))


def make_marcus(comment=""):
    return Marcus(
        A=(1.73e06, "m^3/(mol*s)"), n=2,
        lmbd_i_coefs=[21824.5, -0.0341626, -0.0013254, 4.92966e-07],
        beta=(1.2e10, "1/m"), wr=(0, "kJ/mol"), wp=(0, "kJ/mol"),
        lmbd_o=(0, "J/mol"), comment=comment,
    )


def species():
    return ([Species(label="Lip", molecule=[Molecule(smiles="[Li+]")]),
             Species(label="CH3", molecule=[Molecule(smiles="[CH3]")])],
            [Species(label="CH3Li", molecule=[Molecule(smiles="C[Li]")])])


def register(quarantine):
    """Register the synthetic family in the kinetics database singleton."""
    class _Family(object):
        def __init__(self, q):
            self.quarantine = q
            self.label = FAMILY

    class _Kinetics(object):
        def __init__(self, families):
            self.families = families

    class _Database(object):
        def __init__(self, families):
            self.kinetics = _Kinetics(families)

    import rmgpy.data.rmg
    rmgpy.data.rmg.database = _Database({FAMILY: _Family(quarantine),
                                         "Ordinary_Family": _Family(None)})


def refused(reaction, kinetics):
    try:
        check_quarantine(reaction, stage="the round-87 probe", kinetics=kinetics)
    except QuarantinedKineticsError:
        return True
    return False


# --------------------------------------------------------------------------
# 1. HIGH -- library/seed provenance bypasses the gate
# --------------------------------------------------------------------------
def probe_provenance_bypass(tmpdir):
    with open(os.path.join(tmpdir, QUARANTINE_FILENAME), "w") as f:
        f.write(MANIFEST)
    register(load_family_quarantine(FAMILY, tmpdir))

    reactants, products = species()
    # The provenance RMG itself writes for an estimated rate, verbatim in shape.
    comment = ("Estimated using template [Root_2R->C] for rate rule [Root_2R->C]\n"
               "Euclidian distance = 0\n"
               "family: {0}".format(FAMILY))

    template = TemplateReaction(reactants=reactants, products=products,
                                family=FAMILY, reversible=True)
    template.template = ["Root_2R->C"]
    control = refused(template, make_marcus(comment))

    copied = LibraryReaction(reactants=reactants, products=products,
                            library="copied_seed", kinetics=make_marcus(comment),
                            reversible=True)
    bypassed = not refused(copied, copied.kinetics)

    collide = LibraryReaction(reactants=reactants, products=products,
                             library=FAMILY, kinetics=make_marcus(""), reversible=True)
    false_positive = refused(collide, collide.kinetics)

    record("HIGH: a quarantined rate copied into a library is admitted",
           control and bypassed,
           "control (TemplateReaction of the quarantined family) refused: {0}\n"
           "LibraryReaction(library='copied_seed') carrying the SAME Marcus rate and the\n"
           "  same 'family: {1}' provenance comment admitted: {2}\n"
           "  its .family attribute is {3!r} -- a LIBRARY name in the slot the gate reads\n"
           "  as a FAMILY name (library.py:100, self.family = library)".format(
               control, FAMILY, bypassed, copied.family))

    record("HIGH (second direction): a library named like the family is refused",
           false_positive,
           "LibraryReaction(library='{0}') with NO family provenance at all was refused: {1}\n"
           "  Same slot, opposite error: the gate cannot tell a library label from a family\n"
           "  label, so it both misses real matches and invents false ones.".format(
               FAMILY, false_positive))


# --------------------------------------------------------------------------
# 2. MEDIUM -- the engine pin is weaker than it claims
# --------------------------------------------------------------------------
def loads(tmpdir, extra, tag):
    directory = os.path.join(tmpdir, tag)
    os.makedirs(directory)
    with open(os.path.join(directory, QUARANTINE_FILENAME), "w") as f:
        f.write(MANIFEST + extra)
    try:
        load_family_quarantine(FAMILY, directory)
        return True, "loaded"
    except DatabaseError as exc:
        return False, "refused: {0}".format(str(exc).split("\n")[0][:110])


def probe_pin_strength(tmpdir):
    absent, detail_c = loads(
        tmpdir, 'requiresEngineModule = "rmgpy.no_such_module_at_all"\n', "control")
    control("a module no engine provides is refused", not absent,
            "expected refusal, got: {0}".format(detail_c))

    satisfied, detail_s = loads(
        tmpdir,
        'requiresEngineModule = "rmgpy.data.kinetics.quarantine"\n'
        'requiresEngineSymbol = "check_quarantine"\n', "satisfied")
    control("a requirement this engine DOES meet still loads", satisfied,
            "expected a clean load, got: {0}".format(detail_s))

    pi, detail_pi = loads(
        tmpdir, 'requiresEngineModule = "math"\nrequiresEngineSymbol = "pi"\n', "pi")
    record("MEDIUM: a non-callable satisfies requiresEngineSymbol", pi,
           "manifest requiring math.pi as its gate {0}; math.pi is {1!r}, callable={2}".format(
               detail_pi, math.pi, callable(math.pi)))

    orphan, detail_o = loads(
        tmpdir, 'requiresEngineSymbol = "utterly_absent_symbol"\n', "orphan")
    record("MEDIUM: requiresEngineSymbol without a module is silently ignored", orphan,
           "manifest declaring only a symbol {0} -- _check_engine_requirements returns\n"
           "early on 'if not module_name', so the symbol is never looked up".format(detail_o))

    commit, detail_k = loads(
        tmpdir, 'requiresEngineCommit = "0000000000000000000000000000000000000000"\n',
        "commit")
    record("MEDIUM: requiresEngineCommit is declarable and never checked", commit,
           "manifest pinning an all-zero commit {0}".format(detail_k))

    # Symbol existence does not prove the gate is REACHED. The discriminating case is a
    # module that exists and does NOT bind the symbol: `os` is never going to hold
    # `check_quarantine`, so a manifest naming it as a call site must be refused.
    unwired, detail_w = loads(
        tmpdir,
        'requiresEngineModule = "rmgpy.data.kinetics.quarantine"\n'
        'requiresEngineSymbol = "check_quarantine"\n'
        'requiresEngineCallSites = ("os",)\n', "unwired")
    record("MEDIUM: symbol existence does not prove the gate is reached", unwired,
           "a manifest requiring the gate to be reached from a module that does not bind\n"
           "it {0}".format(detail_w))

    wired, detail_r = loads(
        tmpdir,
        'requiresEngineModule = "rmgpy.data.kinetics.quarantine"\n'
        'requiresEngineSymbol = "check_quarantine"\n'
        'requiresEngineCallSites = ("rmgpy.rmg.model",)\n', "wired")
    control("the real call site satisfies the wiring check", wired,
            "rmgpy.rmg.model binds the same check_quarantine object: {0}".format(detail_r))


# --------------------------------------------------------------------------
# 3. MEDIUM -- the adjacency list in the error message cannot be re-parsed
# --------------------------------------------------------------------------
def probe_message_reparse():
    metastable = Molecule().from_adjacency_list("multiplicity 3\n1 Ar u2 p3 c0\n")
    try:
        saturate_for_estimation(metastable, "transport data")
        record("MEDIUM: the saturated adjacency list re-parses", False,
               "saturation did not raise at all -- the premise of the whole fix is gone")
        return
    except SaturatedStructureError as exc:
        message = str(exc)

    block = message.split("its saturated form was:\n")[1].split("\nand the underlying")[0]

    # The saturated form of a metastable can never satisfy atom typing -- that IS the
    # failure being reported. So the readable-back contract is: it parses with
    # raise_atomtype_exception=False, and the ONLY thing that stops it otherwise is the
    # AtomTypeError the message already quotes. An InvalidAdjacencyListError instead means
    # the block is malformed in a way that has nothing to do with the species.
    try:
        Molecule().from_adjacency_list(block, raise_atomtype_exception=False)
        malformed, detail = False, "the saturated form re-parses when atom typing is relaxed"
    except Exception as exc:
        malformed = True
        detail = ("re-parsing the saturated form with raise_atomtype_exception=False raises\n"
                  "{0}: {1}\nthe block itself was:\n{2}".format(
                      type(exc).__name__, str(exc).split("\n")[0], block.strip()))

    try:
        Molecule().from_adjacency_list(block)
        residual = "nothing -- it parsed even with atom typing enforced"
        honest = False
    except AtomTypeError:
        residual = "AtomTypeError, which is the failure the message is reporting"
        honest = True
    except Exception as exc:
        residual = "{0}, which is NOT what the message reports".format(type(exc).__name__)
        honest = False

    record("MEDIUM: the saturated adjacency list in the message cannot be re-parsed",
           malformed or not honest,
           "{0}\ninsisting on an atom type instead gives {1}".format(detail, residual))

    original = message.split("The species was:\n")[1].split("\nits saturated form")[0]
    try:
        Molecule().from_adjacency_list(original)
        species_ok = True
    except Exception:
        species_ok = False
    control("the SPECIES block in the same message re-parses", species_ok,
            "the unsaturated species is an ordinary molecule and must round-trip")


# --------------------------------------------------------------------------
# 4. LOW -- the shared helper is not universal
# --------------------------------------------------------------------------
def probe_helper_coverage():
    import rmgpy.qm.main
    import rmgpy.data.transport
    import rmgpy.data.solvation
    import rmgpy.data.thermo
    import rmgpy.tools.uncertainty

    routed, raw = [], []
    for module in (rmgpy.data.transport, rmgpy.data.solvation, rmgpy.data.thermo,
                   rmgpy.tools.uncertainty, rmgpy.qm.main):
        source = inspect.getsource(module)
        name = module.__name__
        if "saturate_for_estimation(" in source:
            routed.append(name)
        if ".saturate_radicals()" in source:
            raw.append(name)

    record("LOW: the shared helper is not universal", bool(raw),
           "routed through saturate_for_estimation: {0}\n"
           "still calling saturate_radicals() directly: {1}".format(
               ", ".join(routed) or "none", ", ".join(raw) or "none"))


def main():
    tmpdir = tempfile.mkdtemp(prefix="round87-probe-")
    print("probe scratch: {0}".format(tmpdir))
    for fn, args in ((probe_provenance_bypass, (tmpdir,)),
                     (probe_pin_strength, (tmpdir,)),
                     (probe_message_reparse, ()),
                     (probe_helper_coverage, ())):
        try:
            fn(*args)
        except Exception:
            traceback.print_exc()
            record(fn.__name__, False, "the probe itself raised -- see traceback above")

    print("\n" + "=" * 72 + "\nRESULT\n" + "=" * 72)
    for name, ok, _ in controls:
        print("  {0:<14} control: {1}".format("HOLDS" if ok else "BROKEN", name))
    for name, reproduced, _ in results:
        print("  {0:<14} {1}".format("REPRODUCED" if reproduced else "not reproduced", name))

    broken = [n for n, ok, _ in controls if not ok]
    reproduced = [n for n, r, _ in results if r]
    print("\n{0} of {1} findings reproduced; {2} of {3} controls hold.".format(
        len(reproduced), len(results), len(controls) - len(broken), len(controls)))
    if broken:
        print("CONTROLS BROKEN -- this run measures nothing until they are restored.")
        return 2
    # Exit 1 while any finding is still present, so the same probe is the red before the
    # fix and the green after it.
    return 1 if reproduced else 0


if __name__ == "__main__":
    sys.exit(main())
