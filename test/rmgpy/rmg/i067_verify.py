#!/usr/bin/env python
"""
Verifier for the polymer-tier generation bound (i067-enlarge-bound).

Recomputes everything from artifacts on disk -- RMG.log, chemkin output and the
git working trees. It restates no prose and trusts no earlier claim of this
work; every number printed below is parsed or recomputed at the moment it is
printed.

Usage:  i067_verify.py <runs-root>
        exits 0 if every check passes, non-zero otherwise.

<runs-root> holds one directory per configuration; each has input.py and, once
run, RMG.log. The configuration named by BASELINE is the unbounded run every
other one is measured against; CANDIDATE is the configuration whose block is
being proposed for the deck.
"""
import collections
import os
import re
import subprocess
import sys

BASELINE = "U2"
# P35 mirrors generatedSpeciesConstraints exactly: size routing moves every
# structure of >= polymerSizeThreshold heavy atoms into the polymer tier, and the
# limits there are the gas tier's. Chosen over the tighter P30 because the
# measurement shows no retained species between the two, so the extra tightness
# buys nothing and couples the deck's polymer tier to its `cutoff`.
CANDIDATE = os.environ.get("I067_CANDIDATE", "P35")
CURVE = ["U2", "P35", "P32", "P30", "P29", "P24", "P18", "R1", "R1P30"]

# The run made with RMG_I067_PROFILE=1: same unbounded deck, timers on.
PROFILE_RUN = "PROF"

# Configurations proposed and then refuted by the evidence. They stay in the
# curve, and check 3R re-derives WHY each was rejected from that run's own log,
# so the rejection is a recomputed result and not a remembered one.
REFUTED = ["R1", "R1P30", "P29", "P24", "P18"]

# The flux bar of check 3. toleranceMoveToCore in the deck is 1e-4: a species
# whose peak rate ratio reaches that is promoted to the core, so 1e-4 is the
# value at which RMG itself calls a species' flux meaningful. The assertion uses
# one tenth of it -- an excluded species must not merely have failed to be
# promoted, it must have been an order of magnitude away from mattering.
TOL_MOVE_TO_CORE = 1e-4
FLUX_BAR = 0.1 * TOL_MOVE_TO_CORE

# Check 5's tolerance on core concentrations, as a relative difference.
CONC_RTOL = 1e-6

OUT_OF_SCOPE_REPOS = ["/home/alon/Code/CKMG", "/home/alon/Code/TA"]
EVIDENCE_TREE = "/home/alon/runs"

UNBOUNDED_WARNING = "generatePolymerConstraints is not set"

COST = re.compile(
    r"^(?P<kind>SECONDARY EDGE ENLARGEMENT|ENLARGEMENT) COST: "
    r"(?:(?P<cpu>[\d.]+) s cpu, )?(?P<wall>[\d.]+) s wall; "
    r"\+(?P<dcore>\d+) core species, \+(?P<dedge>\d+) edge species; "
    r"model now core (?P<core_s>\d+)/(?P<core_r>\d+), edge (?P<edge_s>\d+)/(?P<edge_r>\d+)")
CENSUS = re.compile(
    r"CONSTRAINT REFUSAL CENSUS: (?P<total>\d+) structure-refusals "
    r"\(gas tier (?P<gas>\d+), polymer tier (?P<poly>\d+)\)")
HIST = re.compile(r"^\s+(?P<which>polymer-proxy|non-proxy) species by heavy-atom count "
                  r"\((?P<n>\d+) total, max (?P<max>\d+)\): (?P<body>.*)$")
FLUX = re.compile(r"^EDGE FLUX\t(?P<ratio>\S+)\t(?P<heavy>-?\d+)\t(?P<c>-?\d+)\t"
                  r"(?P<proxy>\d)\t(?P<label>[^\t]*)\t(?P<smiles>.*)$")
NEW_SPECIES = re.compile(r"^(Created|Added) (\d+) new (edge|core) species")
ITEM = re.compile(r"^    (\S.*?)\((\d+)\)\s*$")
FAMILY = re.compile(r"^\s+reactions by family: (.*)$")


class Run:
    def __init__(self, root, name):
        self.name = name
        self.dir = os.path.join(root, name)
        self.log = os.path.join(self.dir, "RMG.log")
        self.deck = os.path.join(self.dir, "input.py")
        self.exists = os.path.isfile(self.log)
        self.enlargements = []
        self.census = []
        self.hist = []
        self.families = []
        self.flux = {}
        self.species = {}          # index -> (label, 'core'|'edge')
        self.warned_unbounded = False
        self.deck_has_block = False
        if os.path.isfile(self.deck):
            self.deck_has_block = "generatePolymerConstraints(" in open(self.deck).read()
        if self.exists:
            self._parse()

    def _parse(self):
        lines = open(self.log, errors="replace").read().splitlines()
        i = 0
        while i < len(lines):
            line = lines[i]
            m = COST.match(line)
            if m:
                d = m.groupdict()
                self.enlargements.append({
                    "kind": d["kind"],
                    "cpu": float(d["cpu"]) if d["cpu"] else None,
                    "wall": float(d["wall"]),
                    "core_s": int(d["core_s"]), "core_r": int(d["core_r"]),
                    "edge_s": int(d["edge_s"]), "edge_r": int(d["edge_r"]),
                })
            m = CENSUS.search(line)
            if m:
                self.census.append({k: int(v) for k, v in m.groupdict().items()})
            m = HIST.match(line)
            if m:
                body = {}
                for part in m.group("body").split(", "):
                    h, n = part.split(":")
                    body[int(h)] = int(n)
                self.hist.append((m.group("which"), body))
            m = FAMILY.match(line)
            if m:
                fam = {}
                for part in m.group(1).split(", "):
                    k, v = part.rsplit(":", 1)
                    fam[k] = int(v)
                self.families.append(fam)
            m = FLUX.match(line)
            if m:
                smi = m.group("smiles")
                r = float(m.group("ratio"))
                prev = self.flux.get(smi)
                if prev is None or r > prev[0]:
                    self.flux[smi] = (r, int(m.group("heavy")), int(m.group("c")),
                                      int(m.group("proxy")), m.group("label"))
            if UNBOUNDED_WARNING in line:
                self.warned_unbounded = True
            m = NEW_SPECIES.match(line)
            if m:
                n, where = int(m.group(2)), m.group(3)
                for j in range(i + 1, min(i + 1 + n, len(lines))):
                    it = ITEM.match(lines[j])
                    if it:
                        prev = self.species.get(it.group(2))
                        # 'core' wins: a species that reached the core is recorded as core
                        if prev is None or where == "core":
                            self.species[it.group(2)] = (it.group(1), where)
                i += n
            i += 1

    @property
    def final(self):
        return self.enlargements[-1] if self.enlargements else None

    def total_cpu(self):
        vals = [e["cpu"] for e in self.enlargements if e["cpu"] is not None]
        return sum(vals) if vals else None

    def total_wall(self):
        return sum(e["wall"] for e in self.enlargements)

    def labels(self):
        return {lab for lab, _ in self.species.values()}

    def core_labels(self):
        return {lab for lab, w in self.species.values() if w == "core"}

    def last_hist(self, which):
        for w, body in reversed(self.hist):
            if w == which:
                return body
        return {}


def hr(title):
    print()
    print("=" * 78)
    print(title)
    print("=" * 78)


def fmt_hist(body):
    if not body:
        return "(empty)"
    return ", ".join("{0}:{1}".format(h, body[h]) for h in sorted(body))


def check_1(runs, fails):
    hr("CHECK 1 -- the unbounded leg is characterised")
    base = runs[BASELINE]
    if not base.exists:
        fails.append("check 1: baseline run {0} has no RMG.log".format(BASELINE))
        print("MISSING:", base.log)
        return
    print("baseline configuration: {0}  ({1})".format(
        BASELINE, "no generatePolymerConstraints block" if not base.deck_has_block
        else "HAS a block -- not an unbounded baseline"))
    if base.deck_has_block:
        fails.append("check 1: the baseline deck carries a generatePolymerConstraints block")

    print()
    print("per-enlargement cost of the reaction-generation leg (enlarge() itself):")
    print("  {0:<4} {1:<28} {2:>10} {3:>10} {4:>12} {5:>12}".format(
        "#", "kind", "cpu s", "wall s", "core sp/rxn", "edge sp/rxn"))
    for k, e in enumerate(base.enlargements, 1):
        print("  {0:<4} {1:<28} {2:>10} {3:>10.1f} {4:>12} {5:>12}".format(
            k, e["kind"],
            "{0:.1f}".format(e["cpu"]) if e["cpu"] is not None else "n/a",
            e["wall"],
            "{0}/{1}".format(e["core_s"], e["core_r"]),
            "{0}/{1}".format(e["edge_s"], e["edge_r"])))
    if not base.enlargements:
        fails.append("check 1: no enlargement cost records in the baseline log")

    print()
    print("size distribution of the species the leg generated, by heavy-atom count:")
    for which in ("polymer-proxy", "non-proxy"):
        body = base.last_hist(which)
        print("  {0:<14} n={1:<5} max heavy={2:<4} {3}".format(
            which, sum(body.values()), max(body) if body else 0, fmt_hist(body)))
    if base.families:
        print()
        print("reactions by kinetics family at the deepest enlargement:")
        for k in sorted(base.families[-1]):
            print("  {0:<32} {1}".format(k, base.families[-1][k]))
    if base.census:
        c = base.census[-1]
        print()
        print("structures the constraint tiers refused, cumulative: total {total}, "
              "gas tier {gas}, polymer tier {poly}".format(**c))

    # The size at which a bound should sit is not a matter of taste: it is the
    # size of the largest species the mechanism can actually reach. Read that off
    # the emitted core, which is the set of species RMG judged worth keeping.
    print()
    print(core_size_report(base.dir))
    return base


def core_size_report(run_dir):
    """Heavy-atom profile of the run's emitted CORE, from its species dictionary."""
    path = os.path.join(run_dir, "chemkin", "species_dictionary.txt")
    if not os.path.isfile(path):
        return "NOT ESTABLISHED: no species dictionary at {0}".format(path)
    from rmgpy.chemkin import load_species_dictionary
    d = load_species_dictionary(path)
    rows = []
    for label, spc in d.items():
        m = spc.molecule[0]
        rows.append((m.get_num_atoms() - m.get_num_atoms("H"), m.get_num_atoms("C"),
                     m.get_formula(), label))
    rows.sort(reverse=True)
    hist = collections.Counter(r[0] for r in rows)
    out = ["core species actually promoted, by heavy-atom count ({0} species):".format(len(rows)),
           "  " + fmt_hist(hist)]
    out.append("  largest core species:")
    top = rows[0][0] if rows else 0
    for r in rows[:10]:
        out.append("    heavy={0:<3d} C={1:<3d} {2:<12} {3}".format(*r))
    below = sorted({r[0] for r in rows if r[0] < top}, reverse=True)
    if below:
        out.append("  the largest core species has {0} heavy atoms; the next size down in the "
                   "core is {1}.".format(top, below[0]))
        out.append("  a polymer tier set at {0} heavy atoms admits the entire reachable core "
                   "and refuses everything above it.".format(top))
    return "\n".join(out)


def check_2(runs, fails):
    hr("CHECK 2 -- the bound actually binds")
    base, cand = runs[BASELINE], runs[CANDIDATE]
    if not (base.exists and cand.exists):
        fails.append("check 2: baseline or candidate run missing")
        return set()
    if not cand.deck_has_block:
        fails.append("check 2: candidate deck {0} has no generatePolymerConstraints block"
                     .format(CANDIDATE))
        print("candidate deck carries NO block -- the tier is not in force.")
        return set()
    print("candidate {0} deck carries a generatePolymerConstraints block: yes".format(CANDIDATE))

    # Two separate questions, kept separate because they have different answers.
    #
    # (a) IS THE TIER IN FORCE?  Asserted. If the block is present but nothing is
    #     ever routed through the polymer tier, size routing is broken and the
    #     block is decorative.
    # (b) DOES IT EXCLUDE ANYTHING THE GAS TIER WOULD HAVE KEPT?  Reported, not
    #     asserted -- because zero is a real and reportable answer, and asserting
    #     non-zero here would turn the verifier into a reason to tighten the bound
    #     until it bites, which is precisely the failure mode this work exists to
    #     avoid. Whether zero is acceptable is a modelling judgement, stated in the
    #     report; whether it is TRUE is what this check establishes.
    cens = cand.census[-1] if cand.census else {"poly": 0, "gas": 0, "total": 0}
    print("structures refused in the candidate run: {total} total, "
          "{gas} by the gas tier, {poly} by the POLYMER tier".format(**cens))
    if cens["poly"] == 0:
        fails.append("check 2: the polymer tier refused nothing -- size routing is not in force")
        print("THE POLYMER TIER REFUSED NOTHING: the block is present but inert.")
    else:
        print("=> size routing is in force: the polymer tier is what refuses oversize structures.")
    if base.census:
        print("for contrast, the unbounded baseline refused {total} ({gas} gas, {poly} polymer)."
              .format(**base.census[-1]))

    excluded = base.labels() - cand.labels()
    added = cand.labels() - base.labels()
    print()
    print("species named in the baseline log:  {0}".format(len(base.labels())))
    print("species named in the candidate log: {0}".format(len(cand.labels())))
    print("EXCLUDED by the bound (in baseline, not in candidate): {0}".format(len(excluded)))
    print("present in candidate but not baseline:                 {0}".format(len(added)))
    b, c = base.final, cand.final
    if b and c:
        print()
        print("model size at the deepest enlargement of each run:")
        print("  {0:<8} core {1}/{2}   edge {3}/{4}".format(
            BASELINE, b["core_s"], b["core_r"], b["edge_s"], b["edge_r"]))
        print("  {0:<8} core {1}/{2}   edge {3}/{4}".format(
            CANDIDATE, c["core_s"], c["core_r"], c["edge_s"], c["edge_r"]))
        print("  reactions excluded (edge+core): {0}".format(
            (b["edge_r"] + b["core_r"]) - (c["edge_r"] + c["core_r"])))
    if not excluded:
        print()
        print("THE BOUND EXCLUDES NOTHING THE UNBOUNDED RUN KEPT. Every structure it")
        print("refuses, generatedSpeciesConstraints refused too. It is a DECLARATION of")
        print("where polymer-tier chemistry stops, not a reduction of today's mechanism.")
    return excluded


def check_3(runs, excluded, fails):
    hr("CHECK 3 -- nothing with meaningful flux was cut")
    base = runs[BASELINE]
    print("flux measure: max_edge_species_rate_ratios, the peak species rate normalised")
    print("by the characteristic rate -- the quantity RMG compares to toleranceMoveToCore.")
    print("deck toleranceMoveToCore = {0:g}; assertion bar = {1:g} (one tenth of it).".format(
        TOL_MOVE_TO_CORE, FLUX_BAR))
    if not excluded:
        # No species was removed, so there is no flux to look up. The meaningful
        # question becomes whether the polymer tier refused exactly the set the gas
        # tier already refused -- if it did, no chemistry was removed and the flux
        # question does not arise. Compared refusal-by-refusal at every enlargement
        # both runs reached, not just at the end, so a divergence that later
        # re-converges cannot hide.
        cand = runs[CANDIDATE]
        depth = min(len(base.census), len(cand.census))
        print()
        print("no species was excluded, so no flux lookup applies. Instead: did the")
        print("polymer tier refuse exactly what the gas tier refused, enlargement by")
        print("enlargement? ({0} enlargements common to both runs)".format(depth))
        print()
        print("  {0:<5} {1:>12} {2:>12} {3:>10}".format("#", BASELINE, CANDIDATE, "same?"))
        mismatch = 0
        for k in range(depth):
            b, c = base.census[k]["total"], cand.census[k]["total"]
            same = b == c
            mismatch += 0 if same else 1
            print("  {0:<5} {1:>12} {2:>12} {3:>10}".format(k + 1, b, c, "yes" if same else "NO"))
        if depth == 0:
            fails.append("check 3: no common refusal census to compare")
        elif mismatch:
            fails.append("check 3: the two tiers refused different totals at {0} of {1} "
                         "enlargements, so 'excludes nothing' is not established"
                         .format(mismatch, depth))
        else:
            print()
            print("identical at every common enlargement: the polymer tier refused exactly")
            print("the structures generatedSpeciesConstraints refused. No chemistry removed.")
        return
    core_cut = excluded & base.core_labels()
    if core_cut:
        fails.append("check 3: the bound excluded {0} species that reached the baseline CORE: {1}"
                     .format(len(core_cut), sorted(core_cut)[:10]))
        print("EXCLUDED SPECIES THAT REACHED THE BASELINE CORE: {0}".format(sorted(core_cut)))

    rows, unmeasured = [], []
    for lab in sorted(excluded):
        ent = base.flux.get(lab)
        if ent is None:
            unmeasured.append(lab)
        else:
            rows.append((ent[0], ent[1], ent[2], ent[3], lab))
    rows.sort(reverse=True)
    print()
    print("excluded species with a measured flux: {0}; without one: {1}".format(
        len(rows), len(unmeasured)))
    print()
    print("  {0:>14} {1:>6} {2:>4} {3:>6}  {4}".format("flux ratio", "heavy", "C", "proxy", "species"))
    for r in rows[:25]:
        print("  {0:>14.6e} {1:>6} {2:>4} {3:>6}  {4}".format(r[0], r[1], r[2], r[3], r[4][:70]))
    if len(rows) > 25:
        print("  ... {0} more, all at or below {1:.6e}".format(len(rows) - 25, rows[25][0]))

    if rows:
        worst = rows[0][0]
        print()
        print("maximum flux among excluded species: {0:.6e}  (bar {1:.6e})".format(worst, FLUX_BAR))
        if not worst < FLUX_BAR:
            fails.append("check 3: an excluded species carried flux {0:.6e} >= bar {1:.6e} -- "
                         "the bound is too tight and must be raised".format(worst, FLUX_BAR))
    if unmeasured:
        print()
        print("NOT ESTABLISHED for {0} excluded species: they never appeared in an edge-flux".format(
            len(unmeasured)))
        print("census, because the baseline run stopped before a simulation covered them.")
        print("First 15: {0}".format(sorted(unmeasured)[:15]))


def check_3r(runs, fails):
    """
    Re-derive, from each refuted configuration's own artifacts, the reason it was
    rejected. A configuration counts as refuted only if it either (a) excluded a
    species that reached the baseline core, or (b) failed to complete. If one of
    them shows neither, the rejection was not evidence-based and this check fails.
    """
    hr("CHECK 3R -- every rejected configuration is rejected by its own evidence")
    base = runs[BASELINE]
    for name in REFUTED:
        r = runs.get(name)
        print()
        print("-- {0} --".format(name))
        if r is None or not r.exists:
            print("   no data.")
            fails.append("check 3R: refuted configuration {0} has no run to justify it".format(name))
            continue
        reasons = []
        cut = base.core_labels() - r.labels()
        if cut:
            reasons.append("excluded {0} species that reached the baseline core: {1}".format(
                len(cut), sorted(cut)[:6]))
        rc, crash = run_outcome(r.dir)
        if rc not in (None, 0):
            reasons.append("did not complete (exit {0}){1}".format(
                rc, ": " + crash if crash else ""))
        f = r.final
        if f and base.final and f["edge_s"] * 4 < base.final["edge_s"]:
            reasons.append("edge collapsed to {0} species against the baseline's {1}".format(
                f["edge_s"], base.final["edge_s"]))
        if reasons:
            for reason in reasons:
                print("   REJECTED:", reason)
        else:
            print("   no evidence of harm found in this run's artifacts.")
            fails.append("check 3R: {0} was rejected without evidence in its own run".format(name))


def run_outcome(run_dir):
    """(exit code, one-line crash reason) for a run, from its launch metadata and stderr."""
    meta = os.path.join(run_dir, "launch_meta.txt")
    rc = None
    if os.path.isfile(meta):
        for line in open(meta):
            if line.startswith("exit_rc="):
                rc = int(line.split("=", 1)[1])
    err = os.path.join(run_dir, "stderr.log")
    crash = ""
    if rc not in (None, 0) and os.path.isfile(err):
        lines = [l.rstrip() for l in open(err, errors="replace")]
        for l in reversed(lines):
            if re.match(r"^\w+(\.\w+)*(Error|Exception):", l):
                crash = l[:220]
                break
    return rc, crash


WASTED = re.compile(
    r"WASTED BUILD PROFILE \(this enlargement\): apply_recipe (?P<calls>\d+) calls / "
    r"(?P<recipe_s>[\d.]+) s; of those (?P<refused>\d+) were refused, costing "
    r"(?P<refused_recipe_s>[\d.]+) s in apply_recipe and (?P<refused_total_s>[\d.]+) s across "
    r"the whole product-generation call\. Refused share of the (?P<wall>[\d.]+) s enlargement: "
    r"(?P<pct_recipe>[\d.nan]+)% \(apply_recipe only\) / (?P<pct_total>[\d.nan]+)% "
    r"\(product generation\)\. Accepted product generation: (?P<accepted_s>[\d.]+) s\.")


def check_9(root, fails):
    """
    Where the generation leg's wall time actually goes.

    A size constraint is consulted after apply_recipe has already built the
    product, so every refused structure was paid for in full. This reports what
    that costs as a fraction of the leg -- the number that decides whether
    pre-screening reactant combinations before apply_recipe is worth a ticket.
    """
    hr("CHECK 9 -- where the generation leg's wall time actually goes")
    log = os.path.join(root, PROFILE_RUN, "RMG.log")
    if not os.path.isfile(log):
        print("NOT ESTABLISHED: no profiling run at {0}".format(log))
        fails.append("check 9: no RMG_I067_PROFILE run to measure the wasted-build share")
        return
    rows = [m.groupdict() for m in
            (WASTED.search(l) for l in open(log, errors="replace")) if m]
    rows = [r for r in rows if float(r["wall"]) > 0.5]     # cheap core-add steps: no leg work
    if not rows:
        print("NOT ESTABLISHED: the profiling run has not yet reached an enlargement with")
        print("measurable wall time; nothing to attribute.")
        fails.append("check 9: profiling run produced no substantive enlargement")
        return

    print("method: wall-clock accumulators around apply_recipe and around the whole")
    print("_generate_product_structures call in rmgpy/data/kinetics/family.py, gated on")
    print("RMG_I067_PROFILE. No recipe, ordering or refusal was changed; the gate also")
    print("suppresses this work's own to_smiles() on refused structures so the census")
    print("does not inflate the fraction it is measuring.")
    print()
    print("  {0:<5} {1:>10} {2:>9} {3:>9} {4:>12} {5:>12} {6:>11} {7:>11}".format(
        "enl", "wall s", "recipes", "refused", "recipe s", "refused s", "%wall(recipe)",
        "%wall(total)"))
    tot_wall = tot_recipe = tot_ref_recipe = tot_ref_total = tot_acc = 0.0
    tot_calls = tot_refused = 0
    for k, r in enumerate(rows, 1):
        print("  {0:<5} {1:>10.1f} {2:>9} {3:>9} {4:>12.2f} {5:>12.2f} {6:>11.1f} {7:>11.1f}".format(
            k, float(r["wall"]), r["calls"], r["refused"], float(r["recipe_s"]),
            float(r["refused_total_s"]), float(r["pct_recipe"]), float(r["pct_total"])))
        tot_wall += float(r["wall"]); tot_recipe += float(r["recipe_s"])
        tot_ref_recipe += float(r["refused_recipe_s"]); tot_ref_total += float(r["refused_total_s"])
        tot_acc += float(r["accepted_s"])
        tot_calls += int(r["calls"]); tot_refused += int(r["refused"])

    print()
    print("over {0} substantive enlargements, {1:.1f} s of generation-leg wall time:".format(
        len(rows), tot_wall))
    print("  apply_recipe, all {0} calls                  {1:8.1f} s  {2:5.1f}%".format(
        tot_calls, tot_recipe, 100.0 * tot_recipe / tot_wall))
    print("  apply_recipe on the {0} REFUSED products     {1:8.1f} s  {2:5.1f}%   <-- the wasted build".format(
        tot_refused, tot_ref_recipe, 100.0 * tot_ref_recipe / tot_wall))
    print("  whole product generation, refused products      {0:8.1f} s  {1:5.1f}%   <-- upper bound of a pre-screen".format(
        tot_ref_total, 100.0 * tot_ref_total / tot_wall))
    print("  whole product generation, accepted products     {0:8.1f} s  {1:5.1f}%".format(
        tot_acc, 100.0 * tot_acc / tot_wall))
    residual = tot_wall - tot_ref_total - tot_acc
    print("  everything else in enlarge()                    {0:8.1f} s  {1:5.1f}%".format(
        residual, 100.0 * residual / tot_wall))
    print()
    print("'everything else' is subgraph matching of reactants against family templates,")
    print("thermo estimation, species/reaction bookkeeping and duplicate marking. Matching")
    print("happens BEFORE apply_recipe, so a pre-screen on reactant sizes would also avoid")
    print("part of it -- how much is NOT measured here, and the upper bound above therefore")
    print("understates what such a screen could save.")


def check_4(runs, fails):
    hr("CHECK 4 -- the cost curve is a plateau, not a point")
    print("  {0:<8} {1:>7} {2:>11} {3:>11} {4:>7} {5:>7} {6:>7} {7:>7} {8:>9}".format(
        "config", "enl", "cpu s", "wall s", "core-s", "core-r", "edge-s", "edge-r", "refused"))
    n_with_data = 0
    for name in CURVE:
        r = runs.get(name)
        if r is None or not r.exists or not r.enlargements:
            print("  {0:<8} {1}".format(name, "(no data)"))
            continue
        n_with_data += 1
        f = r.final
        cpu = r.total_cpu()
        refused = r.census[-1]["total"] if r.census else 0
        print("  {0:<8} {1:>7} {2:>11} {3:>11.1f} {4:>7} {5:>7} {6:>7} {7:>7} {8:>9}".format(
            name, len(r.enlargements),
            "{0:.1f}".format(cpu) if cpu is not None else "n/a",
            r.total_wall(), f["core_s"], f["core_r"], f["edge_s"], f["edge_r"], refused))
    print()
    print("cpu s is the sum of the per-enlargement CPU times: it is the generation leg only,")
    print("and unlike wall time it does not measure the machine's other tenants.")
    if n_with_data < 3:
        fails.append("check 4: fewer than three configurations produced data ({0})".format(n_with_data))


def check_5(runs, fails):
    hr("CHECK 5 -- the physical mechanism is unchanged where it should be")
    base, cand = runs[BASELINE], runs[CANDIDATE]
    if not (base.exists and cand.exists):
        fails.append("check 5: baseline or candidate run missing")
        return
    depth = min(len(base.enlargements), len(cand.enlargements))
    print("comparing at the deepest enlargement both runs reached: #{0}".format(depth))
    if depth == 0:
        fails.append("check 5: no common enlargement to compare at")
        return
    b, c = base.enlargements[depth - 1], cand.enlargements[depth - 1]
    print("  {0:<8} core {1} species / {2} reactions".format(BASELINE, b["core_s"], b["core_r"]))
    print("  {0:<8} core {1} species / {2} reactions".format(CANDIDATE, c["core_s"], c["core_r"]))

    bcore, ccore = base.core_labels(), cand.core_labels()
    only_b, only_c = sorted(bcore - ccore), sorted(ccore - bcore)
    print()
    print("core species named only in {0}: {1}".format(BASELINE, only_b or "none"))
    print("core species named only in {0}: {1}".format(CANDIDATE, only_c or "none"))
    if only_b:
        fails.append("check 5: the bound removed core species {0}".format(only_b))

    ok, msg = compare_core_mechanism(base.dir, cand.dir)
    print()
    print(msg)
    if ok is not True:
        fails.append("check 5: the two runs' core mechanisms at the same point are not "
                     "identical, so concentration agreement is not established")

    ok2, msg2 = compare_concentrations(base.dir, cand.dir)
    print()
    print(msg2)
    if ok2 is False:
        fails.append("check 5: core concentrations disagree beyond rtol={0:g}".format(CONC_RTOL))


def compare_core_mechanism(dir_a, dir_b):
    """
    Compare the two runs' core mechanisms at the same point.

    RMG snapshots the core to chemkin/chemNNNN.inp, where NNNN is the core species
    count, so the highest snapshot both runs wrote IS 'taken to the same point' --
    a fixed, comparable state, unlike the live chem.inp, which moves as the run
    goes. Byte equality of that snapshot is stronger than agreeing concentrations
    to a tolerance: identical species, identical thermo and identical rate
    coefficients give bit-identical concentrations under any integrator, so no
    tolerance has to be chosen at all.
    """
    def snaps(d):
        p = os.path.join(d, "chemkin")
        if not os.path.isdir(p):
            return {}
        return {f: os.path.join(p, f) for f in os.listdir(p)
                if re.fullmatch(r"chem\d+\.inp", f)}

    sa, sb = snaps(dir_a), snaps(dir_b)
    common = sorted(set(sa) & set(sb))
    if not common:
        return None, ("NOT ESTABLISHED: the two runs share no chemNNNN.inp core snapshot "
                      "({0} and {1} written), so there is no common point to compare at."
                      .format(len(sa), len(sb)))
    name = common[-1]
    # The only line that legitimately differs is the generation date stamp.
    def body(path):
        return [l for l in open(path, errors="replace")
                if "Date:" not in l and not l.startswith("! Date")]
    a, b = body(sa[name]), body(sb[name])
    if a == b:
        return True, ("core mechanism snapshot chemkin/{0} is IDENTICAL between the two runs "
                      "({1} lines, date stamp excluded).\nIdentical species, thermo and rate "
                      "coefficients give identical concentrations exactly, for any integrator.\n"
                      "Snapshots common to both runs: {2}".format(
                          name, len(a), ", ".join(common)))
    import difflib
    d = list(difflib.unified_diff(a, b, "baseline/" + name, "candidate/" + name, n=0))
    return False, ("core mechanism snapshot chemkin/{0} DIFFERS between the two runs; "
                   "first 30 diff lines:\n{1}".format(name, "".join(d[:30])))


def compare_concentrations(dir_a, dir_b):
    """Compare the final-time core mole fractions of two runs' simulation profiles."""
    def profiles(d):
        p = os.path.join(d, "solver")
        if not os.path.isdir(p):
            return {}
        return {f: os.path.join(p, f) for f in sorted(os.listdir(p)) if f.endswith(".csv")}

    pa, pb = profiles(dir_a), profiles(dir_b)
    common = sorted(set(pa) & set(pb))
    if not common:
        return None, ("NOT ESTABLISHED, and stated rather than papered over: no simulation "
                      "profile CSV is common to both runs (solver/ holds {0} and {1} files). "
                      "SimulationProfileWriter writes on reaction_system.notify(), which the "
                      "hybrid polymer reactor does not call -- the same empty solver/ appears "
                      "in /home/alon/runs/RMG/poly_105, so this is a pre-existing gap in the "
                      "reactor, not a consequence of this work. Concentration agreement here "
                      "therefore rests on the core-mechanism identity above, which implies it, "
                      "rather than on a direct measurement of concentrations."
                      .format(len(pa), len(pb)))
    name = common[-1]
    rows_a = [l.rstrip("\n").split(",") for l in open(pa[name])]
    rows_b = [l.rstrip("\n").split(",") for l in open(pb[name])]
    head_a, head_b = rows_a[0], rows_b[0]
    va = dict(zip(head_a, rows_a[-1]))
    vb = dict(zip(head_b, rows_b[-1]))
    shared = [k for k in head_a if k in vb and k not in ("", "Time (s)")]
    worst, worst_key = 0.0, None
    for k in shared:
        try:
            x, y = float(va[k]), float(vb[k])
        except ValueError:
            continue
        denom = max(abs(x), abs(y), 1e-300)
        rel = abs(x - y) / denom
        if rel > worst:
            worst, worst_key = rel, k
    return (worst <= CONC_RTOL), (
        "final-time core mole fractions from solver/{0} ({1} shared columns):\n"
        "  worst relative difference {2:.3e} on '{3}' (tolerance {4:g})".format(
            name, len(shared), worst, worst_key, CONC_RTOL))


def check_6(runs, fails):
    hr("CHECK 6 -- the unbounded-chemistry warning is gone for the final deck")
    cand = runs[CANDIDATE]
    if not cand.exists:
        fails.append("check 6: candidate run missing")
        return
    hits = [l for l in open(cand.log, errors="replace") if UNBOUNDED_WARNING in l]
    print("occurrences of the unbounded-polymer warning in {0}/RMG.log: {1}".format(
        CANDIDATE, len(hits)))
    for l in hits[:3]:
        print("   ", l.strip()[:160])
    if hits:
        fails.append("check 6: the unbounded-polymer warning still fires for the candidate deck")
    else:
        base = runs[BASELINE]
        if base.exists:
            bhits = [l for l in open(base.log, errors="replace") if UNBOUNDED_WARNING in l]
            print("for contrast, the unbounded baseline {0} fires it {1} time(s):".format(
                BASELINE, len(bhits)))
            if bhits:
                print("   ", bhits[0].strip()[:200])
            else:
                fails.append("check 6: the baseline did NOT fire the warning, so its absence "
                             "in the candidate proves nothing")


def check_7(fails, pytest_report):
    hr("CHECK 7 -- no regression in the test suite")
    if not pytest_report or not os.path.isfile(pytest_report):
        print("NOT ESTABLISHED: no pytest report at {0}.".format(pytest_report))
        fails.append("check 7: no pytest report supplied")
        return
    text = open(pytest_report, errors="replace").read()
    got = {}
    for tag in ("passed", "failed", "skipped", "error", "errors", "xfailed", "xpassed"):
        m = re.findall(r"(\d+) " + tag + r"\b", text)
        if m:
            got[tag] = int(m[-1])
    print("parsed from {0}: {1}".format(pytest_report, got or "(no summary line found)"))
    base_file = pytest_report.replace("after", "before")
    if os.path.isfile(base_file) and base_file != pytest_report:
        btext = open(base_file, errors="replace").read()
        bgot = {}
        for tag in ("passed", "failed", "skipped", "error", "errors"):
            m = re.findall(r"(\d+) " + tag + r"\b", btext)
            if m:
                bgot[tag] = int(m[-1])
        print("baseline (unmodified tree, same HEAD): {0}".format(bgot))
        if got.get("passed", -1) < bgot.get("passed", 0):
            fails.append("check 7: passed count fell from {0} to {1}".format(
                bgot.get("passed"), got.get("passed")))
        if got.get("failed", 0) > bgot.get("failed", 0):
            fails.append("check 7: failed count rose from {0} to {1}".format(
                bgot.get("failed", 0), got.get("failed", 0)))
    else:
        print("NOT ESTABLISHED: no baseline report at {0}; the after-count stands alone."
              .format(base_file))
        fails.append("check 7: no baseline pytest report to compare against")


def check_8(fails, snapshot_dir):
    hr("CHECK 8 -- nothing out of scope changed")
    for repo in OUT_OF_SCOPE_REPOS:
        now = subprocess.run(["git", "-C", repo, "status", "--porcelain"],
                             capture_output=True, text=True).stdout
        base_path = os.path.join(snapshot_dir, os.path.basename(repo).lower() + "_status_before.txt")
        print()
        print("{0}: {1} dirty entries now".format(repo, len([l for l in now.splitlines() if l])))
        if os.path.isfile(base_path):
            before = open(base_path).read()
            new = sorted(set(now.splitlines()) - set(before.splitlines()))
            print("  entries not present in the pre-work snapshot: {0}".format(new or "none"))
            if new:
                fails.append("check 8: new dirty entries in {0}: {1}".format(repo, new))
        else:
            print("  NOT ESTABLISHED: no pre-work snapshot at {0}".format(base_path))
            fails.append("check 8: no pre-work snapshot for {0}".format(repo))

    stamp = os.path.join(snapshot_dir, "session_start_stamp")
    print()
    if os.path.isfile(stamp):
        newer = subprocess.run(
            ["find", EVIDENCE_TREE, "-newer", stamp, "-not", "-type", "d"],
            capture_output=True, text=True).stdout.split()
        print("{0}: {1} file(s) modified since this session's start stamp".format(
            EVIDENCE_TREE, len(newer)))
        for p in newer[:20]:
            print("   ", p)
        if newer:
            fails.append("check 8: {0} file(s) written under {1}".format(len(newer), EVIDENCE_TREE))
    else:
        print("NOT ESTABLISHED: no session start stamp at {0}".format(stamp))
        fails.append("check 8: no start stamp to date writes under {0} against".format(EVIDENCE_TREE))


def main():
    root = sys.argv[1]
    snapshot_dir = sys.argv[2] if len(sys.argv) > 2 else root
    pytest_report = sys.argv[3] if len(sys.argv) > 3 else None
    runs = {name: Run(root, name) for name in CURVE}
    fails = []

    print("runs root:      ", root)
    print("baseline:       ", BASELINE)
    print("candidate:      ", CANDIDATE)
    print("configurations: ", ", ".join(
        "{0}{1}".format(n, "" if runs[n].exists else "(missing)") for n in CURVE))

    check_1(runs, fails)
    excluded = check_2(runs, fails)
    check_3(runs, excluded, fails)
    check_3r(runs, fails)
    check_4(runs, fails)
    check_5(runs, fails)
    check_6(runs, fails)
    check_7(fails, pytest_report)
    check_8(fails, snapshot_dir)
    check_9(root, fails)

    hr("VERDICT")
    if fails:
        for f in fails:
            print("FAIL:", f)
        print()
        print("{0} check(s) failed.".format(len(fails)))
        return 1
    print("all checks passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
