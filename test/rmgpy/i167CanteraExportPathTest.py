#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2026 Prof. William H. Green (whgreen@mit.edu),           #
# Prof. Richard H. West (r.west@neu.edu) and the RMG Team (rmg_dev@mit.edu)   #
#                                                                             #
# Permission is hereby granted, free of charge, to any person obtaining a     #
# copy of this software and associated documentation files (the 'Software'),  #
# to deal in the Software without restriction, including without limitation   #
# the rights to use, copy, modify, merge, publish, distribute, sublicense,    #
# and/or sell copies of the Software, and to permit persons to whom the       #
# Software is furnished to do so, subject to the following conditions:        #
#                                                                             #
# The above copyright notice and this permission notice shall be included in  #
# all copies or substantial portions of the Software.                         #
#                                                                             #
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR  #
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,    #
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE #
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER      #
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     #
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER         #
# DEALINGS IN THE SOFTWARE.                                                   #
#                                                                             #
###############################################################################

"""
The end-of-run Cantera export path, for plasma mechanisms and for ordinary ones.

A plasma mechanism's Chemkin file carries a ``TDEP/<electron>/`` auxiliary line on
every reaction whose rate is evaluated at the electron temperature.  Cantera's
``ck2yaml`` has no ``TDEP`` handler: it falls through to parsing the species name
as a float and raises.  RMG used to run that translation unconditionally, collect
the failure, and escalate it at the end of ``execute`` -- so a plasma run whose
model generated correctly, and every one of whose files was written, exited 1.

Worse, and less visible: the translation and the direct writers' end-of-run
notes-stripping shared one ``try``, so the translation's failure also skipped the
step that produces ``cantera2/chem.yaml``.  A plasma run therefore ended with *no*
usable Cantera file at all -- not the degraded translated one, and not the native
one either.

Every assertion below is made against files on disk, written by the real writers
and read back, and the model under test is built here rather than generated, so
these tests do not depend on which species a flux criterion happens to promote.
"""

import logging
import os

import pytest
import yaml

from rmgpy.chemkin import kinetics_has_plasma_rate, save_chemkin_files, write_kinetics_entry
from rmgpy.exceptions import MechanismWriterError
from rmgpy.kinetics import Arrhenius, BadnellRRArrhenius, MultiArrhenius, VoronovEIArrhenius
from rmgpy.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.rmg.main import RMG, WriterConfig
from rmgpy.rmg.model import CoreEdgeReactionModel, ReactionModel
from rmgpy.species import Species
from rmgpy.thermo import NASA, NASAPolynomial
from rmgpy.transport import TransportData
from rmgpy.yaml_cantera2 import save_cantera_files


def _make_species(label, index, molecule):
    """An RMG Species carrying the thermo and transport both writers need."""
    spc = Species(label=label, molecule=[molecule])
    spc.index = index
    coeffs = [2.5, 0.0, 0.0, 0.0, 0.0, -745.375, -11.7246]
    spc.thermo = NASA(
        polynomials=[
            NASAPolynomial(coeffs=coeffs, Tmin=(200, 'K'), Tmax=(1000, 'K')),
            NASAPolynomial(coeffs=coeffs, Tmin=(1000, 'K'), Tmax=(6000, 'K')),
        ],
        Tmin=(200, 'K'), Tmax=(6000, 'K'),
    )
    spc.transport_data = TransportData(
        shapeIndex=0 if len(molecule.atoms) == 1 else 1,
        sigma=(3.0, 'angstrom'),
        epsilon=(100.0, 'K'),
        dipoleMoment=(0.0, 'De'),
        polarizability=(0.0, 'angstrom^3'),
        rotrelaxcollnum=0.0,
    )
    return spc


def _plasma_core():
    """
    A core model that genuinely contains electron-temperature-dependent kinetics:
    one Voronov electron-impact ionization and one Badnell radiative
    recombination, which are exactly the two rate laws that make the Chemkin
    writer emit a ``TDEP/`` line.  A plain Arrhenius reaction rides along as the
    control that must survive the translation unchanged in the non-plasma case.
    """
    e = _make_species('e', 1, Molecule().from_adjacency_list('1 e u0 p0 c-1'))
    o_atom = _make_species('O', 2, Molecule().from_adjacency_list('1 O u2 p2 c0'))
    o_cation = _make_species('O+', 3, Molecule().from_adjacency_list('1 O u3 p1 c+1'))
    o2 = _make_species('O2', 4, Molecule(smiles='[O][O]'))

    voronov = Reaction(
        index=1, reactants=[o_atom, e], products=[o_cation], electrons=2,
        reversible=False,
        kinetics=VoronovEIArrhenius(
            A=(3.59e-08, 'cm^3/(molecule*s)'), P=0.0, X=0.073, K=0.47,
            dE=(13.62, 'eV'), Tmin=(1.0e+04, 'K'), Tmax=(1.0e+06, 'K'),
        ),
    )
    badnell = Reaction(
        index=2, reactants=[o_cation], products=[o_atom], electrons=-1,
        reversible=False,
        kinetics=BadnellRRArrhenius(
            A=(8.0e-12, 'cm^3/(molecule*s)'), B=0.7,
            T0=(3.0e+02, 'K'), T1=(1.5e+06, 'K'),
            Tmin=(1.0e+04, 'K'), Tmax=(1.0e+06, 'K'),
        ),
    )
    thermal = Reaction(
        index=3, reactants=[o_atom, o_atom], products=[o2], electrons=0,
        reversible=False,
        kinetics=Arrhenius(A=(1.2e+11, 'cm^3/(mol*s)'), n=0.3, Ea=(5000.0, 'J/mol')),
    )
    return ReactionModel(species=[e, o_atom, o_cation, o2],
                         reactions=[voronov, badnell, thermal])


def _thermal_core():
    """The same harness with no plasma rate law anywhere -- the control model."""
    o_atom = _make_species('O', 1, Molecule().from_adjacency_list('1 O u2 p2 c0'))
    o2 = _make_species('O2', 2, Molecule(smiles='[O][O]'))
    thermal = Reaction(
        index=1, reactants=[o_atom, o_atom], products=[o2], electrons=0,
        reversible=False,
        kinetics=Arrhenius(A=(1.2e+11, 'cm^3/(mol*s)'), n=0.3, Ea=(5000.0, 'J/mol')),
    )
    return ReactionModel(species=[o_atom, o2], reactions=[thermal])


def _prepared_rmg(core, output_directory, cantera2_enabled=True):
    """
    An :class:`RMG` positioned exactly where ``execute`` leaves it just before the
    end-of-run Cantera export: the model is final, and both the Chemkin writer and
    (when enabled) the native Cantera writer have already written their
    per-iteration output. Nothing here is a stand-in for the code under test --
    the real writers run.

    ``cantera2_enabled=False`` reproduces the shipped defaults
    (``generateCanteraYAML1=False, generateCanteraYAML2=False``): the Chemkin
    writer runs, but no direct Cantera writer does.
    """
    rmg = RMG()
    rmg.output_directory = output_directory
    rmg.save_edge_species = False
    rmg.export_failures = []
    rmg.chemkin_writer_config = WriterConfig(save_interval=1)
    rmg.cantera2_writer_config = WriterConfig(save_interval=1) if cantera2_enabled else None
    rmg.cantera1_writer_config = None

    rmg.reaction_model = CoreEdgeReactionModel()
    rmg.reaction_model.core = core
    rmg.reaction_model.edge = ReactionModel(species=[], reactions=[])

    for sub in ('chemkin', 'cantera_from_ck', 'cantera2'):
        os.makedirs(os.path.join(output_directory, sub), exist_ok=True)

    save_chemkin_files(rmg, config=rmg.chemkin_writer_config)
    if cantera2_enabled:
        save_cantera_files(rmg, config=rmg.cantera2_writer_config)
    return rmg


@pytest.fixture()
def plasma_rmg(tmp_path):
    return _prepared_rmg(_plasma_core(), str(tmp_path))


@pytest.fixture()
def plasma_rmg_default_config(tmp_path):
    """A plasma run under the shipped defaults: no direct Cantera writer enabled."""
    return _prepared_rmg(_plasma_core(), str(tmp_path), cantera2_enabled=False)


@pytest.fixture()
def thermal_rmg(tmp_path):
    return _prepared_rmg(_thermal_core(), str(tmp_path))


def _phase(yaml_path):
    with open(yaml_path) as f:
        return yaml.safe_load(f)['phases'][0]


class TestPlasmaExportPath:
    """The plasma case: the translation is skipped and the native writer stands."""

    def test_the_mechanism_really_carries_tdep(self, plasma_rmg):
        """
        Guards every other test in this class.  If the Chemkin file stopped
        carrying ``TDEP/``, these tests would pass while exercising nothing, and
        that is precisely how a regression test rots into a decoration.
        """
        with open(os.path.join(plasma_rmg.output_directory, 'chemkin', 'chem.inp')) as f:
            text = f.read()
        assert 'TDEP/' in text, 'the model under test no longer produces a TDEP line'

    def test_export_records_no_failure(self, plasma_rmg):
        """
        The run must not be told that an export failed.  ``export_failures`` is
        what ``report_export_failures`` escalates into the non-zero exit.
        """
        plasma_rmg.generate_end_of_run_cantera_files()
        assert plasma_rmg.export_failures == [], (
            'a plasma export recorded a failure: '
            + '; '.join(step for step, _ in plasma_rmg.export_failures))

    def test_native_cantera_file_is_produced(self, plasma_rmg):
        """
        ``cantera2/chem.yaml`` is the authoritative artifact, and it is produced
        end-of-run by a step that used to share a ``try`` with the translation.
        """
        plasma_rmg.generate_end_of_run_cantera_files()
        path = os.path.join(plasma_rmg.output_directory, 'cantera2', 'chem.yaml')
        assert os.path.exists(path), 'the authoritative Cantera file was not written'

    def test_native_cantera_file_is_a_plasma(self, plasma_rmg):
        """The file that is produced is the plasma, not a degraded ideal gas."""
        plasma_rmg.generate_end_of_run_cantera_files()
        phase = _phase(os.path.join(plasma_rmg.output_directory, 'cantera2', 'chem.yaml'))
        assert phase['thermo'] == 'plasma'
        assert phase['transport'] == 'ionized-gas'
        assert 'electron-energy-distribution' in phase

    def test_translation_is_not_attempted(self, plasma_rmg):
        """No translated file, and no half-written one either."""
        plasma_rmg.generate_end_of_run_cantera_files()
        ck_dir = os.path.join(plasma_rmg.output_directory, 'cantera_from_ck')
        assert os.listdir(ck_dir) == [], (
            'the Chemkin-to-Cantera translation left files behind: '
            + repr(os.listdir(ck_dir)))

    def test_the_skip_is_logged_and_explained(self, plasma_rmg, caplog):
        """
        Skipping is not silence.  A run that quietly omits an output someone
        expected is the failure this whole path exists to remove, so the log must
        say that the step was skipped, why, and what stands in its place.
        """
        with caplog.at_level(logging.INFO):
            plasma_rmg.generate_end_of_run_cantera_files()
        text = caplog.text
        assert 'ck2yaml' in text
        assert 'TDEP' in text
        assert 'cantera2' in text
        assert 'skip' in text.lower()


class TestNonPlasmaExportPathIsUntouched:
    """The control: a model with no plasma rate law goes through unchanged."""

    def test_translation_still_runs(self, thermal_rmg):
        thermal_rmg.generate_end_of_run_cantera_files()
        assert thermal_rmg.export_failures == []
        assert os.path.exists(
            os.path.join(thermal_rmg.output_directory, 'cantera_from_ck', 'chem.yaml'))

    def test_native_writer_still_runs_alongside_it(self, thermal_rmg):
        thermal_rmg.generate_end_of_run_cantera_files()
        path = os.path.join(thermal_rmg.output_directory, 'cantera2', 'chem.yaml')
        assert os.path.exists(path)
        assert _phase(path)['thermo'] == 'ideal-gas'

    def test_nothing_is_logged_about_skipping(self, thermal_rmg, caplog):
        with caplog.at_level(logging.INFO):
            thermal_rmg.generate_end_of_run_cantera_files()
        assert 'TDEP' not in caplog.text


class TestDefaultConfigProducesNoArtifactAndSaysSo:
    """
    F3: under the shipped defaults (``generateCanteraYAML1=False,
    generateCanteraYAML2=False``) both direct Cantera writers are off, so a plasma
    run's skipped translation leaves no Cantera artifact at all. The run must say
    so and record an export failure -- not advertise an authoritative file it never
    wrote and then exit 0.
    """

    def test_records_an_export_failure(self, plasma_rmg_default_config):
        rmg = plasma_rmg_default_config
        rmg.generate_end_of_run_cantera_files()
        assert rmg.export_failures, (
            'a default-config plasma run produced no Cantera artifact yet recorded '
            'no export failure')

    def test_the_recorded_failure_fails_the_run(self, plasma_rmg_default_config):
        rmg = plasma_rmg_default_config
        rmg.generate_end_of_run_cantera_files()
        with pytest.raises(MechanismWriterError):
            rmg.report_export_failures()

    def test_the_message_names_the_missing_artifact_and_the_fix(self, plasma_rmg_default_config, caplog):
        rmg = plasma_rmg_default_config
        with caplog.at_level(logging.INFO):
            rmg.generate_end_of_run_cantera_files()
        text = caplog.text
        assert 'No Cantera artifact was produced' in text
        assert 'generateCanteraYAML2' in text
        # It must NOT claim the absent file is authoritative.
        assert 'authoritative' not in text

    def test_no_cantera_file_exists_on_disk(self, plasma_rmg_default_config):
        rmg = plasma_rmg_default_config
        rmg.generate_end_of_run_cantera_files()
        for sub in ('cantera_from_ck', 'cantera1', 'cantera2'):
            d = os.path.join(rmg.output_directory, sub)
            if os.path.isdir(d):
                assert os.listdir(d) == [], f'{sub}/ unexpectedly holds {os.listdir(d)!r}'


class TestCantera2EnabledStillReportsSuccess:
    """Verifier 2: the correct config must not regress. cantera2 on -> no failure,
    and the existing authoritative-artifact message stands."""

    def test_no_failure_and_message_kept(self, plasma_rmg, caplog):
        with caplog.at_level(logging.INFO):
            plasma_rmg.generate_end_of_run_cantera_files()
        assert plasma_rmg.export_failures == []
        assert 'authoritative' in caplog.text


class TestStaleTranslationArtifactIsRemoved:
    """
    F2: a plasma run must not leave -- and must never post-process -- a
    ``cantera_from_ck/chem.yaml`` it did not produce. Reusing an output directory
    from a previous non-plasma run leaves a stale translated file there; the skip
    path used to strip its transport notes so it looked freshly made and exit 0.
    """

    def test_planted_stale_chem_yaml_is_gone(self, plasma_rmg):
        stale = os.path.join(plasma_rmg.output_directory, 'cantera_from_ck', 'chem.yaml')
        with open(stale, 'w') as f:
            f.write('# stale translated mechanism from a previous non-plasma run\n')
        plasma_rmg.generate_end_of_run_cantera_files()
        assert not os.path.exists(stale), (
            'a stale cantera_from_ck/chem.yaml survived the plasma skip path')

    def test_planted_stale_chem_edge_yaml_is_gone(self, plasma_rmg):
        stale = os.path.join(plasma_rmg.output_directory, 'cantera_from_ck', 'chem_edge.yaml')
        with open(stale, 'w') as f:
            f.write('# stale edge translation from a previous run\n')
        plasma_rmg.generate_end_of_run_cantera_files()
        assert not os.path.exists(stale)

    def test_stale_file_is_removed_not_stripped_in_place(self, plasma_rmg):
        """The failure mode was ``strip_yaml_notes`` rewriting the stale file in
        place. Prove it is removed rather than rewritten: its sentinel never
        survives, in any form, at that path."""
        stale = os.path.join(plasma_rmg.output_directory, 'cantera_from_ck', 'chem.yaml')
        with open(stale, 'w') as f:
            f.write('# SENTINEL stale mechanism\n')
        plasma_rmg.generate_end_of_run_cantera_files()
        assert not os.path.exists(stale)

    def test_planted_stale_chem_annotated_yaml_is_gone(self, plasma_rmg):
        """F1, the enumerated gap. ``chem_annotated.yaml`` is the *other* file the
        translation writes (``generate_cantera_files_from_chemkin`` translates both
        ``chem.inp`` and ``chem_annotated.inp``). The skip cleanup used to delete
        ``chem.yaml`` and a ``chem_edge.yaml`` the translation never writes, and
        leave this one behind for a downstream reader to take as current output."""
        stale = os.path.join(plasma_rmg.output_directory, 'cantera_from_ck', 'chem_annotated.yaml')
        with open(stale, 'w') as f:
            f.write('# stale annotated translation from a previous non-plasma run\n')
        plasma_rmg.generate_end_of_run_cantera_files()
        assert not os.path.exists(stale), (
            'a stale cantera_from_ck/chem_annotated.yaml survived the plasma skip path')

    def test_every_genuinely_translated_file_is_gone_after_a_reused_dir(self, tmp_path):
        """Not a planted sentinel: a real thermal run writes the translation into a
        directory, then a plasma run reuses it. Every file the thermal translation
        actually produced must be gone -- none may be read as this run's output.
        This is the two-run reproduction the finding demands, made permanent."""
        d = str(tmp_path)
        thermal = _prepared_rmg(_thermal_core(), d, cantera2_enabled=True)
        thermal.generate_end_of_run_cantera_files()
        ck = os.path.join(d, 'cantera_from_ck')
        produced = sorted(os.listdir(ck))
        assert produced, 'the thermal run wrote no translated files -- the test would be vacuous'
        plasma = _prepared_rmg(_plasma_core(), d, cantera2_enabled=True)
        plasma.generate_end_of_run_cantera_files()
        survivors = [f for f in produced if os.path.exists(os.path.join(ck, f))]
        assert survivors == [], (
            'translated files from the previous run survived the skip: ' + repr(survivors))

    def test_stale_comparison_report_is_gone(self, plasma_rmg):
        """F2. A ``comparison_report.txt`` describes a comparison of the translated
        copy against a direct writer's file. The skip performs no translation and
        no comparison, so a prior run's report must not survive to describe a
        comparison this run did not perform."""
        report = os.path.join(plasma_rmg.output_directory, 'cantera2', 'comparison_report.txt')
        with open(report, 'w') as f:
            f.write('# stale comparison report from a previous run\n')
        plasma_rmg.generate_end_of_run_cantera_files()
        assert not os.path.exists(report), (
            'a stale cantera2/comparison_report.txt survived the plasma skip path')


def _writer_species():
    """An electron, an O atom, its cation and O2 -- enough for
    ``write_kinetics_entry`` to name an electron in a TDEP line and to write a
    charge-balanced thermal control reaction."""
    e = _make_species('e', 1, Molecule().from_adjacency_list('1 e u0 p0 c-1'))
    o_atom = _make_species('O', 2, Molecule().from_adjacency_list('1 O u2 p2 c0'))
    o_cation = _make_species('O+', 3, Molecule().from_adjacency_list('1 O u3 p1 c+1'))
    o2 = _make_species('O2', 4, Molecule(smiles='[O][O]'))
    return e, o_atom, o_cation, o2


class TestPredicateMatchesWriterUnderMultiWrappers:
    """
    F9: the predicate and the Chemkin writer are one rule
    (:func:`rmgpy.chemkin.kinetics_has_plasma_rate`), so a plasma rate law buried
    in a ``MultiArrhenius`` wrapper -- which the writer unwraps and emits ``TDEP/``
    for -- is also seen by ``mechanism_has_plasma_kinetics``. Both directions.
    """

    def _voronov(self):
        return VoronovEIArrhenius(
            A=(3.59e-08, 'cm^3/(molecule*s)'), P=0.0, X=0.073, K=0.47,
            dE=(13.62, 'eV'), Tmin=(1.0e+04, 'K'), Tmax=(1.0e+06, 'K'))

    def _plain(self):
        return Arrhenius(A=(1.2e+11, 'cm^3/(mol*s)'), n=0.3, Ea=(5000.0, 'J/mol'))

    def test_wrapped_plasma_subrate_predicate_and_writer_agree_true(self):
        e, o_atom, o_cation, _o2 = _writer_species()
        rxn = Reaction(index=1, reactants=[o_atom, e], products=[o_cation],
                       electrons=2, reversible=False,
                       kinetics=MultiArrhenius(arrhenius=[self._voronov()]))
        assert kinetics_has_plasma_rate(rxn.kinetics) is True
        entry = write_kinetics_entry(rxn, species_list=[e, o_atom, o_cation])
        assert 'TDEP/' in entry, 'writer did not emit TDEP for a wrapped plasma subrate'

    def test_wrapped_plain_subrate_predicate_and_writer_agree_false(self):
        e, o_atom, o_cation, o2 = _writer_species()
        rxn = Reaction(index=1, reactants=[o_atom, o_atom], products=[o2],
                       electrons=0, reversible=False,
                       kinetics=MultiArrhenius(arrhenius=[self._plain(), self._plain()]))
        assert kinetics_has_plasma_rate(rxn.kinetics) is False
        entry = write_kinetics_entry(rxn, species_list=[e, o_atom, o_cation, o2])
        assert 'TDEP/' not in entry, 'writer emitted TDEP for a non-plasma wrapped rate'

    def test_mechanism_predicate_sees_a_wrapped_plasma_reaction(self):
        e, o_atom, o_cation, _o2 = _writer_species()
        rxn = Reaction(index=1, reactants=[o_atom, e], products=[o_cation],
                       electrons=2, reversible=False,
                       kinetics=MultiArrhenius(arrhenius=[self._voronov()]))
        rmg = RMG()
        rmg.reaction_model = CoreEdgeReactionModel()
        rmg.reaction_model.core = ReactionModel(species=[e, o_atom, o_cation], reactions=[rxn])
        assert rmg.mechanism_has_plasma_kinetics() is True
