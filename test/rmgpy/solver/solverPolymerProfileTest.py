"""Polymer decks must write simulation profiles.

The HybridPolymerReactor wrapper is what RMG.register_listeners attaches
SimulationProfileWriter to, but the engine it builds is what calls
notify() inside simulate(). Observers attached to the wrapper must
therefore reach the engine."""

import csv
import os

from rmgpy.molecule import Molecule
from rmgpy.quantity import Quantity
from rmgpy.rmg.listener import SimulationProfileWriter
from rmgpy.rmg.polymer_input import HybridPolymerReactor, PolymerPhase, PolymerPool
from rmgpy.rmg.settings import ModelSettings, SimulatorSettings
from rmgpy.species import Species


def _spc(smiles, label):
    s = Species(molecule=[Molecule().from_smiles(smiles)])
    s.label = label
    return s


def _reactor():
    n2 = _spc("N#N", "N2")
    p2 = _spc("CC", "P2")
    mu0, mu1, mu2 = _spc("CO", "poly_mu0"), _spc("C=O", "poly_mu1"), _spc("C#N", "poly_mu2")
    for m in (mu0, mu1, mu2):
        m.is_moment_dummy = True
    pool = PolymerPool(label="poly", xs=2, monomer=Molecule().from_smiles("C=C"),
                       explicit_map={2: p2}, mu_species=[mu0, mu1, mu2],
                       proxy_species=_spc("CCCC", "poly"))
    phase = PolymerPhase(density=Quantity(1000.0, "kg/m^3"),
                         initial_moments={"poly": (1.0, 5.0, 30.0)},
                         initial_explicit={p2: 0.3, n2: 1.0},
                         pools=[pool], mass_transfer=[])
    reactor = HybridPolymerReactor(temperature=(800.0, "K"), pressure=(1.0e5, "Pa"),
                                   initialMoles={n2: 1.0}, polymerPhase=phase,
                                   terminationTime=(1.0, "s"))
    return reactor, [n2, p2, mu0, mu1, mu2]


def _simulate(reactor, core):
    ms = ModelSettings(tol_keep_in_edge=0.0, tol_move_to_core=1.0e-3,
                       tol_interrupt_simulation=1.0e8)
    reactor.simulate(core, [], [], [], surface_species=[], surface_reactions=[],
                     model_settings=ms,
                     simulator_settings=SimulatorSettings())


def test_profile_writer_attached_to_wrapper_writes_csv(tmp_path):
    reactor, core = _reactor()
    os.makedirs(tmp_path / "solver")
    reactor.initialize_model(core, [], [], [])
    # attached AFTER the engine exists: must still reach it
    reactor.attach(SimulationProfileWriter(str(tmp_path), 0, core))
    _simulate(reactor, core)

    with open(tmp_path / "solver" / "simulation_1_5.csv") as f:
        rows = list(csv.reader(f))
    header, data = rows[0], rows[1:]
    assert header[0] == "Time (s)"
    assert {"poly_mu0", "poly_mu1", "poly_mu2"} <= set(header)
    assert len(data) > 1
    assert float(data[-1][0]) > float(data[0][0])
    assert len(data[0]) == len(header)


def test_observer_attached_before_engine_build_reaches_engine(tmp_path):
    reactor, core = _reactor()
    os.makedirs(tmp_path / "solver")
    reactor.attach(SimulationProfileWriter(str(tmp_path), 0, core))
    _simulate(reactor, core)
    assert (tmp_path / "solver" / "simulation_1_5.csv").stat().st_size > 0
