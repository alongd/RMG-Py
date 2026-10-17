"""Test support (I-150): a linear alkyl radical with the radical on an interior
carbon as a window-provider spec, ``("alkane", n, head_end, tail_end, radical)``
(``radical``: 0-based carbon of the radical site).  ``rmgpy.kmc.windows`` knows
alkane specs only with end chemistries; ``install()`` teaches ``build_spec``,
``longer_specs`` and ``event_specs`` (module attributes, process-wide) to carry
the interior radical, shifting it with a head extension.  No product code."""

from rmgpy.kmc import boundaries as b
from rmgpy.kmc import windows as win

_ORIGINAL = {}


def build_alkyl_radical(carbons, radical, head_end=None, tail_end=None, name=None):
    from rmgpy.kmc.rmg_moves import molecule_from_rows
    ends = {0: int(head_end == "primary_radical"), carbons - 1: int(tail_end == "primary_radical")}
    ends[radical] = ends.get(radical, 0) + 1
    neighbours = {i: [] for i in range(carbons)}
    for i in range(carbons - 1):
        neighbours[i].append(i + 1)
        neighbours[i + 1].append(i)
    h_of, next_id = {}, carbons
    for i in range(carbons):
        count = 4 - len(neighbours[i]) - ends.get(i, 0)
        h_of[i] = tuple(range(next_id, next_id + count))
        next_id += count
    rows = {}
    for i in range(carbons):
        edges = [(j, 1.0) for j in neighbours[i]] + [(h, 1.0) for h in h_of[i]]
        rows[i] = ("C", None, 0, ends.get(i, 0), 0, tuple(sorted(edges)))
        for h in h_of[i]:
            rows[h] = ("H", None, 0, 0, 0, ((i, 1.0),))
    molecule = molecule_from_rows(rows)

    def port(port_id, atom, tag, end):
        if end is None:
            return b.Port(port_id, atom, (h_of[atom][-1],), tag, "artificial")
        return b.Port(port_id, atom, (), tag, "physical")
    ports = (port("head", 0, "backbone_head", head_end),
             port("tail", carbons - 1, "backbone_tail", tail_end))
    metadata = b.BoundaryMetadata("linear_alkane", ports,
                                  tuple((i, i, 0) for i in range(carbons)))
    metadata.check(rows, name or f"C{carbons}")
    return b.ProxyBuild(name or f"C{carbons}r{radical}", molecule, metadata, tuple(range(next_id)))


def _interior(spec):
    return isinstance(spec, tuple) and spec[0] == "alkane" and len(spec) > 4


def install():
    if _ORIGINAL:
        return
    _ORIGINAL.update(build_spec=win.build_spec, longer_specs=win.longer_specs,
                     event_specs=win.event_specs)

    def build_spec(spec, name="p"):
        if _interior(spec):
            return build_alkyl_radical(spec[1], spec[4], spec[2], spec[3])
        return _ORIGINAL["build_spec"](spec, name)

    def longer_specs(specs, ports, k=1):
        out, shift = _ORIGINAL["longer_specs"](specs, ports, k)
        for i, spec in enumerate(specs):
            if _interior(spec):
                out[i] = tuple(out[i][:4]) + (spec[4] + shift.get(i, 0),)
        return out, shift

    def event_specs(specs, assignment):
        longer, shift = win.longer_specs(specs, [p for p, c in assignment if c == win.INTERIOR])
        for (index, port_id), choice in assignment:
            if choice == win.INTERIOR:
                continue
            spec = longer[index]
            if isinstance(spec, tuple) and spec[0] == "alkane":
                ends = list(spec[2:4]) + [None] * (4 - len(spec))
                ends[0 if port_id == "head" else 1] = choice
                longer[index] = ("alkane", spec[1], *ends) + tuple(spec[4:5])
                continue
            from dataclasses import replace
            if not isinstance(spec, b.PsChainSpec):
                return None, shift
            longer[index] = replace(spec, **{"head_end" if port_id == "head" else "tail_end":
                                             choice})
        return longer, shift

    win.build_spec, win.longer_specs, win.event_specs = build_spec, longer_specs, event_specs
    win.builds_of = lambda specs: [win.build_spec(spec, f"p{k}") for k, spec in enumerate(specs)]
