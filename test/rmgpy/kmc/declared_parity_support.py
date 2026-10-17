"""Test support (I-159): parity of ``DeclaredChainProvider`` with the legacy
polystyrene and linear-alkane providers.  No product code.

The legacy providers are the reference: every legacy spec is mapped to a
``ChainSpec`` on a declaration (``legacy_to_chain``,
``legacy_alkane_to_chain``) and the two builds are compared as graphs whose
node labels carry the atom attributes, the repeat coordinate (legacy
coordinates mapped through the position map pi, which is computed once by
isomorphism of the DP1 builds, never hand-written), the ports (id, tag,
physicality) and the caps.  Isomorphism of those graphs maps coordinates
through pi, ports to ports and caps to caps (up to the cap automorphism)."""

import hashlib
import json
from dataclasses import dataclass

import networkx as nx
from networkx.algorithms import isomorphism

from rmgpy.kmc import boundaries as b
from rmgpy.kmc import declared_chain as dc
from rmgpy.kmc import windows as win
from rmgpy.polymer_spec import LinearHomopolymerSpec

import interior_radical_specs as irs
from melt_engine_support import pe_mapping
from meltStateTest import ps_mapping

E = dc.EditOp


# --------------------------------------------------------------------------
# declarations
# --------------------------------------------------------------------------
def ps_declaration_mapping():
    """The PS declaration of ``meltStateTest.ps_mapping`` with its styrene
    small species also declared as the ``monomer_product`` (ps_mapping
    declares none; the legacy monomer build needs one)."""
    mapping = ps_mapping()
    mapping["monomer_product"] = json.loads(
        json.dumps(mapping["small_species"]["styrene"]["graph"]))
    return mapping


def ps_declaration():
    return LinearHomopolymerSpec.from_mapping(ps_declaration_mapping())


def pe_declaration():
    return LinearHomopolymerSpec.from_mapping(pe_mapping({10: 1}, transport=False))


# --------------------------------------------------------------------------
# graph comparison
# --------------------------------------------------------------------------
def _port_marks(metadata, prefix=()):
    marks = {}
    for port in metadata.ports:
        marks.setdefault(port.atom_id, []).append(
            prefix + ("port", port.port_id, port.tag, port.physicality))
        for cap in port.cap_atom_ids:
            marks.setdefault(cap, []).append(prefix + ("cap", port.port_id))
    return marks


def build_graph(build, coord_map=None):
    """networkx graph of a ProxyBuild: label = (attributes, mapped coordinate,
    port and cap marks); edge label = bond order."""
    rows = dict(b.rows_of_molecule(build.molecule))
    coords = {a: (u, k) for a, u, k in build.metadata.coordinates}
    marks = _port_marks(build.metadata)
    graph = nx.Graph()
    for a, row in rows.items():
        coord = coords.get(a)
        if coord is not None and coord_map is not None:
            coord = coord_map(coord)
        graph.add_node(a, label=repr((tuple(row[:5]), coord, tuple(sorted(marks.get(a, ()))))))
    for a, row in rows.items():
        for c, order in row[5]:
            if a < c:
                graph.add_edge(a, c, label=repr(float(order)))
    return graph


def _matcher(g1, g2):
    return isomorphism.GraphMatcher(g1, g2, node_match=lambda x, y: x["label"] == y["label"],
                                    edge_match=lambda x, y: x["label"] == y["label"])


def graphs_isomorphic(g1, g2) -> bool:
    if sorted(d["label"] for _n, d in g1.nodes(data=True)) != \
            sorted(d["label"] for _n, d in g2.nodes(data=True)):
        return False
    if sorted(d["label"] for _a, _c, d in g1.edges(data=True)) != \
            sorted(d["label"] for _a, _c, d in g2.edges(data=True)):
        return False
    return _matcher(g1, g2).is_isomorphic()


def builds_match(legacy, declared, coord_map=None) -> bool:
    return graphs_isomorphic(build_graph(legacy, coord_map), build_graph(declared))


def witness_graph(witness, coord_map=None):
    before, after = witness.reactant_rows(), witness.product_rows()
    support = b.edit_support(witness)
    roles, coords, marks = {}, {}, {}
    for participant in witness.participants:
        for a in participant.atom_ids:
            roles[a] = participant.index
        for a, u, k in participant.metadata.coordinates:
            coords[a] = (u, k) if coord_map is None else coord_map((u, k))
        for a, items in _port_marks(participant.metadata, (participant.index,)).items():
            marks.setdefault(a, []).extend(items)
    graph = nx.Graph()
    for a in before:
        graph.add_node(a, label=repr((tuple(before[a][:5]), tuple(after[a][:5]), roles.get(a),
                                      coords.get(a), tuple(sorted(marks.get(a, ()))),
                                      a in support)))
    bb, ba = witness.reactant_bonds(), witness.product_bonds()
    for pair in set(bb) | set(ba):
        graph.add_edge(*pair, label=repr((bb.get(pair, 0.0), ba.get(pair, 0.0))))
    return graph


def witnesses_match(legacy, declared, coord_map=None) -> bool:
    return graphs_isomorphic(witness_graph(legacy, coord_map), witness_graph(declared))


def fragment_graph(fragment, coord_map=None):
    coords = {i: (u, k) for i, u, k in fragment.coordinates}
    graph = nx.Graph()
    graph.add_node("old", label="old")
    for i, atom in enumerate(fragment.atoms):
        coord = coords.get(i)
        if coord is not None and coord_map is not None:
            coord = coord_map(coord)
        graph.add_node(i, label=repr((tuple(atom), coord, i == fragment.new_port,
                                      i in fragment.new_caps)))
    for i, j, order in fragment.bonds:
        graph.add_edge(i, j, label=repr(float(order)))
    for i, order in fragment.attachments:
        graph.add_edge("old", i, label=repr(float(order)))
    return graph


def fragments_match(legacy, declared, coord_map=None) -> bool:
    return graphs_isomorphic(fragment_graph(legacy, coord_map), fragment_graph(declared))


# --------------------------------------------------------------------------
# the PS position map (pi) and the PS adapter
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class PsContext:
    """Declared provider, legacy provider, pi (legacy position -> declared
    position) and the key of every legacy position."""

    provider: object
    legacy: object
    pi: tuple
    keys: tuple
    phase: int

    def coord_map(self, coord):
        unit, position = coord
        return unit, self.pi[position]


def derive_position_map(legacy_build, declared_build) -> dict:
    """pi from one isomorphism of the DP1 builds (ports and caps marked,
    coordinates stripped)."""
    g1 = build_graph(legacy_build, coord_map=lambda c: None)
    g2 = build_graph(declared_build, coord_map=lambda c: None)
    matcher = _matcher(g1, g2)
    if not matcher.is_isomorphic():
        raise AssertionError("the DP1 builds are not isomorphic")
    mapping = matcher.mapping
    legacy = {a: (u, k) for a, u, k in legacy_build.metadata.coordinates}
    declared = {a: (u, k) for a, u, k in declared_build.metadata.coordinates}
    pi = {}
    for a, (u, k) in legacy.items():
        du, dk = declared[mapping[a]]
        assert du == u, (a, u, du)
        pi[k] = dk
    return pi


def ps_context(spec=None) -> PsContext:
    provider = dc.DeclaredChainProvider(spec or ps_declaration())
    legacy = b.PROVIDERS["polystyrene"]
    pi = derive_position_map(legacy.build("dp1", b.PsChainSpec(1)),
                             provider.build("dp1", dc.ChainSpec(1)))
    geometry = provider.geometry
    keys = tuple(geometry.keys[pi[k]] for k in range(len(pi)))
    phase = next(i for i, ph in enumerate(geometry.phases, 1)
                 if set(ph.bond) == {keys[0], keys[1]})
    return PsContext(provider, legacy, tuple(pi[k] for k in range(len(pi))), keys, phase)


class AdapterRefusal(b.BoundaryError):
    """The legacy spec names a site the declared chain does not have."""


END_OF = {"methyl": "declared", "benzyl": "declared", "primary_radical": "radical",
          "benzylic_radical": "radical", "vinylidene": "declared", "styryl": "declared"}
MONOMER = "monomer"


def legacy_to_chain(spec, ctx: PsContext):
    """The ChainSpec of a PsChainSpec (every field), or MONOMER."""
    if spec.monomer:
        return MONOMER
    K = ctx.keys
    units = spec.units + int(spec.tail_shift)
    head = END_OF[spec.head_end] if spec.head_end else \
        "radical" if spec.radical_end else "artificial"
    tail = END_OF[spec.tail_end] if spec.tail_end else \
        "radical" if spec.benzylic_end else "artificial"
    edits = []
    for u in spec.radical_units:
        edits += [E("remove_h", u, K[1]), E("set_radical", u, K[1], 1)]
    for u in spec.secondary_radical_units:
        edits += [E("remove_h", u, K[0]), E("set_radical", u, K[0], 1)]
    for u in spec.pendant_loss_units:
        edits += [E("delete_atom", u, K[k]) for k in range(2, 8)]
        edits.append(E("set_radical", u, K[1], 1))
    if spec.ring_radical_unit is not None:
        u = spec.ring_radical_unit
        edits += [E("remove_h", u, K[5]), E("set_radical", u, K[5], 1)]
    backbone = []
    for u in range(spec.units):
        if not (u == 0 and spec.head_shift):
            backbone.append((u, 0))
        backbone.append((u, 1))
    if spec.tail_shift:
        backbone.append((spec.units, 0))
    double = []
    for end, pair in ((spec.head_end, backbone[:2]), (spec.tail_end, backbone[-2:])):
        if end in b.ALKENE_ENDS:
            if len(pair) < 2:
                raise AdapterRefusal(f"{end}: one backbone carbon")
            if sorted(pair) not in double:
                double.append(sorted(pair))
    for (u1, k1), (u2, k2) in double:
        edits += [E("remove_h", u1, K[k1]), E("remove_h", u2, K[k2]),
                  E("set_bond_order", u1, K[k1], 2, other=(u2, K[k2]))]
    return dc.ChainSpec(units, head=head, tail=tail,
                        head_phase=ctx.phase if spec.head_shift else 0,
                        tail_phase=ctx.phase if spec.tail_shift else 0, edits=tuple(edits))


def declared_ps_build(ctx: PsContext, spec, name="d"):
    chain = legacy_to_chain(spec, ctx)
    if chain == MONOMER:
        return ctx.provider.monomer(name)
    return ctx.provider.build(name, chain)


def ladder_shapes(max_dp=3):
    """Every species shape ``windows.ladder_plan`` builds for DP <= max_dp
    (intact, one-carbon, phase A/B scission, pendant loss, C-H cleavage)."""
    S = b.PsChainSpec
    out = []
    for dp in range(1, max_dp + 1):
        for head in b.HEAD_ENDS:
            for tail in b.TAIL_ENDS:
                out.append(S(dp, head_end=head, tail_end=tail))
        if dp == 1:
            for head in b.TAIL_ENDS:
                for tail in b.TAIL_ENDS:
                    out.append(S(1, head_shift=True, head_end=head, tail_end=tail))
            for head in b.HEAD_ENDS:
                for tail in b.HEAD_ENDS:
                    out.append(S(0, tail_shift=True, head_end=head, tail_end=tail))
        for head in b.HEAD_ENDS:
            for tail in b.HEAD_ENDS:
                out.append(S(dp, head_end=head, tail_end=tail, tail_shift=True))
        for head in b.TAIL_ENDS:
            for tail in b.TAIL_ENDS:
                out.append(S(dp, head_end=head, tail_end=tail, head_shift=True))
        for u in range(dp):
            if 2 * u + 1 < 2 * dp - 1:
                out.append(S(dp, head_end="methyl", tail_end="benzyl", pendant_loss_units=(u,)))
                out.append(S(dp, head_end="methyl", tail_end="benzyl", radical_units=(u,)))
            if u >= 1:
                out.append(S(dp, head_end="methyl", tail_end="benzyl",
                             secondary_radical_units=(u,)))
    return out


def ps_catalog() -> dict:
    """group -> [PsChainSpec, ...] (union of the legacy catalogs, each spec
    once, first group wins)."""
    seen, groups = set(), {}

    def add(group, spec):
        if repr(spec) not in seen:
            seen.add(repr(spec))
            groups.setdefault(group, []).append(spec)
    for units in (3, 5):
        for spec in b.ps_structure_specs(units).values():
            add("structure", spec)
    for kind in win.PARTICIPANT_KINDS:
        for dp in range(1, 5):
            for spec in win.kind_specs(kind, dp):
                add("kind", spec)
    for make in win.LONG_SPECS.values():
        for spec in make():
            add("long", spec)
    for make in win.ONE_CARBON.values():
        add("one_carbon", make())
    for spec in ladder_shapes(3):
        add("ladder", spec)
    return groups


def ps_parity(ctx: PsContext, specs):
    """(mismatches, refusals): a spec the legacy provider refuses must be
    refused by the declared side too."""
    mismatches, refused = [], []
    for spec in specs:
        try:
            legacy = ctx.legacy.build("l", spec)
        except b.BoundaryError as error:
            try:
                declared_ps_build(ctx, spec)
            except b.BoundaryError:
                refused.append((spec, str(error)))
                continue
            mismatches.append((spec, f"legacy refuses ({error}); declared builds"))
            continue
        try:
            declared = declared_ps_build(ctx, spec)
        except b.BoundaryError as error:
            mismatches.append((spec, f"declared refuses: {error}"))
            continue
        if spec.monomer:
            # the monomer positions live in the monomer's own key space: only
            # the unit is compared (no ports on either side)
            unit_only = lambda c: (c[0],)  # noqa: E731
            ok = graphs_isomorphic(build_graph(legacy, unit_only),
                                   build_graph(declared, unit_only))
        else:
            ok = builds_match(legacy, declared, ctx.coord_map)
        if not ok:
            mismatches.append((spec, "not marked-isomorphic"))
    return mismatches, refused


# --------------------------------------------------------------------------
# the alkane adapter
# --------------------------------------------------------------------------
ALKANE_END = {None: "artificial", "methyl": "declared", "primary_radical": "radical"}


def legacy_alkane_to_chain(spec, key):
    """``("alkane", n[, head, tail[, radical]])`` on the PE declaration
    (primitive period one CH2, declaration key ``key``)."""
    n = spec[1]
    head = spec[2] if len(spec) > 2 else None
    tail = spec[3] if len(spec) > 3 else None
    edits = ()
    if len(spec) > 4:
        r = spec[4]
        value = 1 + int(r == 0 and head == "primary_radical") + \
            int(r == n - 1 and tail == "primary_radical")
        edits = (E("remove_h", r, key), E("set_radical", r, key, value))
    return dc.ChainSpec(n, head=ALKANE_END[head], tail=ALKANE_END[tail], edits=edits)


def legacy_alkane_build(spec):
    n = spec[1]
    head = spec[2] if len(spec) > 2 else None
    tail = spec[3] if len(spec) > 3 else None
    if len(spec) > 4:
        return irs.build_alkyl_radical(n, spec[4], head, tail)
    return b.PROVIDERS["linear_alkane"].build(n, None, head, tail)


# --------------------------------------------------------------------------
# witnesses without RMG generation
# --------------------------------------------------------------------------
def identity_witness(build):
    molecules, participants = b.witness_reactants([build])
    rows = tuple(sorted(b.rows_of_molecule(molecules[0])))
    return b.make_witness(rows, rows, participants)


def cut_witness(build, first, second):
    """A homolysis-shaped witness: the bond between the heavy atoms at
    coordinates ``first`` and ``second`` is broken on the product side, each
    end gains a radical."""
    molecules, participants = b.witness_reactants([build])
    rows = dict(b.rows_of_molecule(molecules[0]))
    at = {(u, k): a for a, u, k in participants[0].metadata.coordinates}
    x, y = at[first], at[second]
    products = {}
    for a, row in rows.items():
        edges = tuple(e for e in row[5] if {a, e[0]} != {x, y})
        radical = row[3] + int(a in (x, y))
        products[a] = (row[0], row[1], row[2], radical, row[4], edges)
    return b.make_witness(tuple(sorted(rows.items())), tuple(sorted(products.items())),
                          participants)


# --------------------------------------------------------------------------
# determinism fingerprint (run in subprocesses under two hash seeds)
# --------------------------------------------------------------------------
def build_fingerprint(build) -> str:
    rows = tuple(sorted(b.rows_of_molecule(build.molecule)))
    return repr((rows, build.canonical_metadata(), len(build.construction_order)))


def fingerprint_chains(geometry):
    keys = [k for k in geometry.keys]
    return (dc.ChainSpec(1), dc.ChainSpec(3, head="declared", tail="radical"),
            dc.ChainSpec(4, head_phase=min(1, len(geometry.phases)),
                         edits=(E("remove_h", 2, keys[0]), E("set_radical", 2, keys[0], 1))))


def fingerprint(spec=None) -> str:
    parts = []
    for declaration in ((spec,) if spec is not None else (ps_declaration(), pe_declaration())):
        provider = dc.DeclaredChainProvider(declaration)
        parts.append(provider.name)
        for chain in fingerprint_chains(provider.geometry):
            parts.append(build_fingerprint(provider.build("f", chain)))
        parts.append(repr(provider.geometry.atom_classes))
    return hashlib.sha256("\n".join(parts).encode()).hexdigest()
