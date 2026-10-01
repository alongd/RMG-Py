#!/usr/bin/env python3

###############################################################################
#                                                                             #
# RMG - Reaction Mechanism Generator                                          #
#                                                                             #
# Copyright (c) 2002-2023 Prof. William H. Green (whgreen@mit.edu),           #
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

import logging
import os
import shutil
import textwrap
from unittest.mock import patch
import tempfile

import pandas as pd
import pytest

from rmgpy import get_path, settings
from rmgpy.data.rmg import RMGDatabase
from rmgpy.rmg.main import RMG, RMG_Memory, initialize_log, make_profile_graph
from rmgpy.rmg.model import CoreEdgeReactionModel

originalPath = get_path()


@pytest.mark.functional
def test_external_library_seed_restart_round_trip(tmp_path, monkeypatch):
    """A seed restart reopens its external library by absolute source path."""
    _external_library_restart_round_trip(tmp_path, monkeypatch, three_reactions=False)


@pytest.mark.functional
def test_external_library_three_reaction_core_and_edge_restart(tmp_path, monkeypatch):
    """Core and restart_edge share a single load of a three-reaction source."""
    _external_library_restart_round_trip(tmp_path, monkeypatch, three_reactions=True)


def _external_library_restart_round_trip(tmp_path, monkeypatch, three_reactions):
    library = tmp_path / "external" / "tiny-external"
    library.mkdir(parents=True)
    (library / "reactions.py").write_text(textwrap.dedent("""
        name = "tiny-external"
        entry(index=1, label="H2 + O2 <=> H + HO2",
              kinetics=Arrhenius(A=(1.0e6, 'm^3/(mol*s)'), n=0, Ea=(0, 'kJ/mol')))
    """))
    (library / "dictionary.txt").write_text(textwrap.dedent("""
        H2
        1 H u0 p0 c0 {2,S}
        2 H u0 p0 c0 {1,S}

        O2
        1 O u1 p2 c0 {2,S}
        2 O u1 p2 c0 {1,S}

        H
        1 H u1 p0 c0

        HO2
        1 O u1 p2 c0 {2,S}
        2 O u0 p2 c0 {1,S} {3,S}
        3 H u0 p0 c0 {2,S}
    """))
    if three_reactions:
        with (library / 'reactions.py').open('a') as f:
            f.write("entry(index=2, label='H2 + O2 <=> OH + OH', "
                    "kinetics=Arrhenius(A=(1.0e6, 'm^3/(mol*s)'), n=0, Ea=(0, 'kJ/mol')))\n")
            f.write("entry(index=3, label='H2O2 + H2 <=> H2O + H2O', "
                    "kinetics=Arrhenius(A=(1.0e6, 'm^3/(mol*s)'), n=0, Ea=(0, 'kJ/mol')))\n")
        with (library / 'dictionary.txt').open('a') as f:
            f.write(textwrap.dedent("""
                OH
                1 O u1 p2 c0 {2,S}
                2 H u0 p0 c0 {1,S}

                H2O2
                1 O u0 p2 c0 {2,S} {3,S}
                2 O u0 p2 c0 {1,S} {4,S}
                3 H u0 p0 c0 {1,S}
                4 H u0 p0 c0 {2,S}

                H2O
                1 O u0 p2 c0 {2,S} {3,S}
                2 H u0 p0 c0 {1,S}
                3 H u0 p0 c0 {1,S}
            """))
    source = os.path.realpath(str(library))
    first = tmp_path / "first"
    first.mkdir()
    input_file = first / "input.py"
    input_file.write_text(textwrap.dedent("""
        database(thermoLibraries=['primaryThermoLibrary'], reactionLibraries=[{source!r}],
                 seedMechanisms=[], kineticsDepositories=['training'],
                 kineticsFamilies=['H_Abstraction'], kineticsEstimator='rate rules')
        species(label='H2', reactive=True, structure=SMILES('[H][H]'))
        species(label='O2', reactive=True, structure=SMILES('[O][O]'))
        simpleReactor(temperature=(1000,'K'), pressure=(1,'bar'),
                      initialMoleFractions={{'H2': .67, 'O2': .33}}, terminationTime=(1e-12,'s'))
        simulator(atol=1e-16, rtol=1e-8)
        model(toleranceKeepInEdge=0, toleranceMoveToCore=1e-3, toleranceInterruptSimulation=1e-3)
        options(name='external-restart', generateOutputHTML=False, generatePlots=False,
                saveSimulationProfiles=False, saveEdgeSpecies=False)
    """.format(source=source)))
    RMG(input_file=str(input_file), output_directory=str(first)).execute()
    assert source in (first / 'seed' / 'seed' / 'reactions.py').read_text()

    restart_dir = tmp_path / "restart"
    restart_dir.mkdir()
    restart_input = first / "restart_from_seed.py"
    monkeypatch.chdir(restart_dir)
    restart = RMG(input_file=str(restart_input), output_directory=str(restart_dir))
    restart.execute()
    reactions = restart.database.kinetics.libraries['restart'].get_library_reactions()
    assert restart.database.kinetics.external_library_labels[source] == 'tiny-external'
    assert all(rxn.library != source for rxn in reactions if hasattr(rxn, 'library'))
    assert any(rxn.library == 'tiny-external' for rxn in reactions)

    seed_only_input = tmp_path / 'seed_only_restart.py'
    seed_only_input.write_text(
        restart_input.read_text()
        .replace("restartFromSeed(path='seed')", "restartFromSeed(path={!r})".format(str(first / 'seed')))
        .replace("reactionLibraries=[{!r}]".format(source), "reactionLibraries=[]")
    )
    seed_only = tmp_path / 'seed-only'
    seed_only.mkdir()
    restarted_from_seed = RMG(input_file=str(seed_only_input), output_directory=str(seed_only))
    if three_reactions:
        from rmgpy.data.kinetics.library import KineticsLibrary
        loads = []
        original_load = KineticsLibrary.load

        def tracked_load(self, filename, *args, **kwargs):
            if os.path.realpath(filename) == os.path.join(source, 'reactions.py'):
                loads.append(filename)
            return original_load(self, filename, *args, **kwargs)

        with patch.object(KineticsLibrary, 'load', tracked_load):
            restarted_from_seed.execute()
            # Startup loads restart_edge but does not currently admit it. Exercise
            # its public admission path against the actual first job's seed output.
            edge = restarted_from_seed.database.kinetics.libraries['restart_edge'].get_library_reactions()
            assert any(getattr(rxn, 'library', None) == 'tiny-external' for rxn in edge)
            restarted_from_seed.reaction_model.add_reaction_library_to_edge('restart_edge')
        assert len(loads) == 1
        assert len(restarted_from_seed.database.kinetics.libraries['tiny-external'].entries) == 3
        assert sum(label == 'tiny-external' for label, kind in
                   restarted_from_seed.database.kinetics.library_order) == 1
    else:
        restarted_from_seed.execute()
    assert restarted_from_seed.database.kinetics.external_library_labels[source] == 'tiny-external'

    shutil.move(str(library), str(tmp_path / 'moved-external'))
    missing = tmp_path / 'missing'
    missing.mkdir()
    with pytest.raises(IOError, match='External kinetics library .* recorded by restart seed is missing'):
        RMG(input_file=str(seed_only_input), output_directory=str(missing)).execute()


@pytest.mark.functional
class TestMain:
    @classmethod
    def setup_class(cls):
        """A function that is run ONCE before all unit tests in this class."""
        cls.testDir = os.path.join(originalPath, "..", "test", "rmgpy", "test_data", "mainTest")
        cls.outputDir = "output"
        cls.databaseDirectory = settings["database.directory"]

        # Read database content through symlinks, but keep export destinations
        # local so saveSeedToDatabase exercises the real writer safely.
        cls.seedExport = tempfile.TemporaryDirectory(prefix='rmg-main-seeds-')
        database_root = cls.seedExport.name
        for name in os.listdir(cls.databaseDirectory):
            if name != 'kinetics':
                os.symlink(os.path.join(cls.databaseDirectory, name), os.path.join(database_root, name))
        kinetics_root = os.path.join(database_root, 'kinetics')
        os.mkdir(kinetics_root)
        for name in os.listdir(os.path.join(cls.databaseDirectory, 'kinetics')):
            if name != 'libraries':
                os.symlink(os.path.join(cls.databaseDirectory, 'kinetics', name), os.path.join(kinetics_root, name))
        libraries_root = os.path.join(kinetics_root, 'libraries')
        os.mkdir(libraries_root)
        for name in os.listdir(os.path.join(cls.databaseDirectory, 'kinetics', 'libraries')):
            if name not in ('testSeed', 'testSeed_edge'):
                os.symlink(os.path.join(cls.databaseDirectory, 'kinetics', 'libraries', name),
                           os.path.join(libraries_root, name))
        cls.seedKinetics = os.path.join(libraries_root, 'testSeed')
        cls.seedKineticsEdge = os.path.join(libraries_root, 'testSeed_edge')
        assert not os.path.exists(cls.seedKinetics)
        assert not os.path.exists(cls.seedKineticsEdge)

        output_path = os.path.join(cls.testDir, cls.outputDir)
        if os.path.exists(output_path):
            shutil.rmtree(output_path)
        os.mkdir(output_path)

        cls.rmg = RMG(input_file=os.path.join(cls.testDir, 'input.py'), output_directory=output_path)
        with patch.dict(settings, {'database.directory': database_root}):
            cls.rmg.execute()

    @classmethod
    def teardown_class(cls):
        """A function that is run ONCE after all unit tests in this class."""
        # Reset module level database
        import rmgpy.data.rmg

        rmgpy.data.rmg.database = None

        # Remove output directory
        shutil.rmtree(os.path.join(cls.testDir, cls.outputDir))

        cls.seedExport.cleanup()

    def test_rmg_execute(self):
        """Test that RMG.execute completed successfully."""
        assert isinstance(self.rmg.database, RMGDatabase)
        assert self.rmg.done

    def test_rmg_increases_reactions(self):
        """Test that RMG.execute increases reactions and species."""
        assert len(self.rmg.reaction_model.core.reactions) > 0
        assert len(self.rmg.reaction_model.core.species) > 1
        assert len(self.rmg.reaction_model.edge.reactions) > 0
        assert len(self.rmg.reaction_model.edge.species) > 0

    def test_rmg_seed_mechanism_creation(self):
        """Test that the expected seed mechanisms are created in output directory."""
        seed_dir = os.path.join(self.testDir, self.outputDir, "seed")
        assert os.path.exists

        assert os.path.exists(os.path.join(seed_dir, "seed"))  # kinetics library folder made

        assert os.path.exists(os.path.join(seed_dir, "seed", "dictionary.txt"))  # dictionary file made
        assert os.path.exists(os.path.join(seed_dir, "seed", "reactions.py"))  # reactions file made

    def test_rmg_seed_edge_mechanism_creation(self):
        """Test that the expected seed mechanisms are created in output directory."""
        seed_dir = os.path.join(self.testDir, self.outputDir, "seed")
        assert os.path.exists

        assert os.path.exists(os.path.join(seed_dir, "seed_edge"))  # kinetics library folder made

        assert os.path.exists(os.path.join(seed_dir, "seed_edge", "dictionary.txt"))  # dictionary file made
        assert os.path.exists(os.path.join(seed_dir, "seed_edge", "reactions.py"))  # reactions file made

    def test_rmg_seed_library_creation(self):
        """Test that seed mechanisms are created in the correct database locations."""
        assert self.rmg.save_seed_to_database
        assert os.path.isfile(os.path.join(self.seedKinetics, 'reactions.py'))
        assert os.path.isfile(os.path.join(self.seedKinetics, 'dictionary.txt'))

    def test_rmg_seed_edge_library_creation(self):
        """Test that edge seed mechanisms are created in the correct database locations."""
        assert self.rmg.save_seed_to_database
        assert os.path.isfile(os.path.join(self.seedKineticsEdge, 'reactions.py'))
        assert os.path.isfile(os.path.join(self.seedKineticsEdge, 'dictionary.txt'))

    def test_rmg_rms_mechanism_files_creation(self):
        """Test that rms mechanisms are created in the correct location."""
        assert os.path.exists(os.path.join(self.rmg.output_directory,"rms"))
        assert len(os.listdir(os.path.join(self.rmg.output_directory,"rms"))) != 0
        
    def test_rmg_seed_works(self):
        """Test that the created seed libraries work.

        Note: Since this test modifies the class level RMG instance,
        it can cause other tests to fail if run out of order."""
        # Load the seed libraries into the database
        self.rmg.database.load(
            path=self.databaseDirectory,
            thermo_libraries=[],
            reaction_libraries=[self.seedKinetics, self.seedKineticsEdge],
            seed_mechanisms=[self.seedKinetics, self.seedKineticsEdge],
            kinetics_families="default",
            kinetics_depositories=[],
            depository=False,
        )

        self.rmg.reaction_model = CoreEdgeReactionModel()
        self.rmg.reaction_model.add_reaction_library_to_edge("testSeed")  # try adding seed as library
        assert len(self.rmg.reaction_model.edge.species) > 0
        assert len(self.rmg.reaction_model.edge.reactions) > 0

        self.rmg.reaction_model = CoreEdgeReactionModel()
        self.rmg.reaction_model.add_seed_mechanism_to_core("testSeed")  # try adding seed as seed mech
        assert len(self.rmg.reaction_model.core.species) > 0
        assert len(self.rmg.reaction_model.core.reactions) > 0

        self.rmg.reaction_model = CoreEdgeReactionModel()
        self.rmg.reaction_model.add_reaction_library_to_edge("testSeed_edge")  # try adding seed as library
        assert len(self.rmg.reaction_model.edge.species) > 0
        assert len(self.rmg.reaction_model.edge.reactions) > 0

        self.rmg.reaction_model = CoreEdgeReactionModel()
        self.rmg.reaction_model.add_seed_mechanism_to_core("testSeed_edge")  # try adding seed as seed mech
        assert len(self.rmg.reaction_model.core.species) > 0
        assert len(self.rmg.reaction_model.core.reactions) > 0

    def test_rmg_memory(self):
        """
        test that RMG Memory objects function properly
        """
        for rxnsys in self.rmg.reaction_systems:
            Rmem = RMG_Memory(rxnsys, None)
            Rmem.generate_cond()
            Rmem.get_cond()
            Rmem.add_t_conv_N(1.0, 0.2, 2)
            Rmem.generate_cond()
            Rmem.get_cond()

    def test_make_cantera_input_file_from_ck(self):
        """
        This test ensures that a usable Cantera input file is created via the Chemkin to Cantera conversion.
        """
        import cantera as ct

        cantera_files = os.path.join(self.rmg.output_directory, "cantera_from_ck")
        files = os.listdir(cantera_files)
        for f in files:
            if ".yaml" in f:
                try:
                    ct.Solution(os.path.join(cantera_files, f))
                except:
                    assert False, "The output Cantera file is not loadable in Cantera."
    
    def test_make_cantera_input_file_directly_1(self):
        """
        This tests to ensure that a usable Cantera input file is created via direct yaml writer 1.
        """
        import cantera as ct

        cantera_files = os.path.join(self.rmg.output_directory, "cantera1")
        files = os.listdir(cantera_files)
        for f in files:
            if ".yaml" in f:
                try:
                    ct.Solution(os.path.join(cantera_files, f))
                except:
                    assert False, "The output Cantera file is not loadable in Cantera."

    def test_make_cantera_input_file_directly_2(self):
        """
        This tests to ensure that a usable Cantera input file is created via direct yaml writer 2.
        """
        import cantera as ct

        cantera_files = os.path.join(self.rmg.output_directory, "cantera2")
        files = os.listdir(cantera_files)
        for f in files:
            if ".yaml" in f:
                try:
                    ct.Solution(os.path.join(cantera_files, f))
                except:
                    assert False, "The output Cantera file is not loadable in Cantera."

    def test_cantera_input_files_match_chemkin_later(self):
        """
        Copy the Cantera YAML files (generated directly by RMG and converted from Chemkin)
        to the test data directory so that yaml_cantera1Test can compare them.
        """
        # Copy RMG-generated YAML 1 to test data directory
        cantera_dir = os.path.join(self.rmg.output_directory, "cantera1")
        rmg_yaml_path = os.path.join(cantera_dir, 'chem_annotated.yaml')
        assert os.path.exists(rmg_yaml_path), f"RMG-generated Cantera YAML file {rmg_yaml_path} not found"
        test_data_cantera_target = os.path.join(self.testDir, '..', 'yaml_writer_data', 'cantera1', 'from_main_test.yaml')
        os.makedirs(os.path.dirname(test_data_cantera_target), exist_ok=True)
        shutil.copy(rmg_yaml_path, test_data_cantera_target)

        # Copy RMG-generated YAML 2 to test data directory
        cantera_dir = os.path.join(self.rmg.output_directory, "cantera2")
        rmg_yaml_path = os.path.join(cantera_dir, 'chem_annotated.yaml')
        assert os.path.exists(rmg_yaml_path), f"RMG-generated Cantera YAML file {rmg_yaml_path} not found"
        test_data_cantera_target = os.path.join(self.testDir, '..', 'yaml_writer_data', 'cantera2', 'from_main_test.yaml')
        os.makedirs(os.path.dirname(test_data_cantera_target), exist_ok=True)
        shutil.copy(rmg_yaml_path, test_data_cantera_target)

        # Copy chemkin-converted YAML to test data directory
        cantera_from_ck_dir = os.path.join(
            self.rmg.output_directory, "cantera_from_ck"
        )
        ck_yaml_path = os.path.join(cantera_from_ck_dir, "chem_annotated.yaml")
        assert os.path.exists(ck_yaml_path), f"Chemkin-converted YAML file {ck_yaml_path} not found"
        test_data_chemkin_target = os.path.join(self.testDir, '..', 'yaml_writer_data', 'ck2yaml', 'from_main_test.yaml')
        os.makedirs(os.path.dirname(test_data_chemkin_target), exist_ok=True)
        shutil.copy(ck_yaml_path, test_data_chemkin_target)


@pytest.mark.functional
class TestRestartWithFilters:
    @classmethod
    def setup_class(cls):
        """A function that is run ONCE before all unit tests in this class."""
        cls.testDir = os.path.join(originalPath, "..", "test", "rmgpy", "test_data", "restartTest")
        cls.outputDir = os.path.join(cls.testDir, "output_w_filters")
        cls.databaseDirectory = settings["database.directory"]

        try:
            os.mkdir(cls.outputDir)
        except FileExistsError:
            pass  # output directory already exists
        initialize_log(logging.INFO, os.path.join(cls.outputDir, "RMG.log"))

        cls.rmg = RMG(
            input_file=os.path.join(cls.testDir, "restart_w_filters.py"),
            output_directory=os.path.join(cls.outputDir),
        )

    def test_restart_with_filters(self):
        """
        Test that the RMG restart job with filters included completed without problems
        """
        self.rmg.execute()
        with open(os.path.join(self.outputDir, "RMG.log"), "r") as f:
            assert "MODEL GENERATION COMPLETED" in f.read()

    @classmethod
    def teardown_class(cls):
        """A function that is run ONCE after all unit tests in this class."""
        # Reset module level database
        import rmgpy.data.rmg

        rmgpy.data.rmg.database = None

        # Remove output directory
        shutil.rmtree(cls.outputDir)


@pytest.mark.functional
class TestRestartNoFilters:
    @classmethod
    def setup_class(cls):
        """A function that is run ONCE before all unit tests in this class."""
        cls.testDir = os.path.join(originalPath, "..", "test", "rmgpy", "test_data", "restartTest")
        cls.outputDir = os.path.join(cls.testDir, "output_no_filters")
        cls.databaseDirectory = settings["database.directory"]

        os.mkdir(cls.outputDir)
        initialize_log(logging.INFO, os.path.join(cls.outputDir, "RMG.log"))

        cls.rmg = RMG(
            input_file=os.path.join(cls.testDir, "restart_no_filters.py"),
            output_directory=os.path.join(cls.outputDir),
        )

    def test_restart_no_filters(self):
        """
        Test that the RMG restart job with no filters included completed without problems
        """
        self.rmg.execute()
        with open(os.path.join(self.outputDir, "RMG.log"), "r") as f:
            assert "MODEL GENERATION COMPLETED" in f.read()

    @classmethod
    def teardown_class(cls):
        """A function that is run ONCE after all unit tests in this class."""
        # Reset module level database
        import rmgpy.data.rmg

        rmgpy.data.rmg.database = None

        # Remove output directory
        shutil.rmtree(cls.outputDir)


@pytest.mark.functional
class TestMainFunctions:
    @classmethod
    def setup_class(cls):
        """A function that is run ONCE before all unit tests in this class."""
        cls.testDir = os.path.join(originalPath, "..", "test", "rmgpy", "test_data", "mainTest")
        cls.outputDir = os.path.join(cls.testDir, "output")
        cls.databaseDirectory = settings["database.directory"]

        os.makedirs(os.path.join(cls.testDir, cls.outputDir), exist_ok=True)

        cls.max_iter = 10

        cls.rmg = RMG(
            input_file=os.path.join(cls.testDir, "superminimal_input.py"),
            output_directory=cls.outputDir,
        )

        cls.rmg.execute(max_iterations=cls.max_iter)

    def test_save_seed_modulus(self):
        """
        Test that saveSeedModulus argument from superminimal_input.py saved the correct number of seeds
        """
        path = os.path.join(self.outputDir, "previous_seeds")
        num_dir_actual = sum(os.path.isdir(os.path.join(path, i)) for i in os.listdir(path))
        num_dir_expected = self.max_iter // 2 + 1  # +1 is for saving iteration 0
        assert num_dir_actual == num_dir_expected

    def test_max_iter(self):
        """
        Test the command line argument of -i
        """
        df = pd.read_excel(os.path.join(self.outputDir, "statistics.xls"))
        num_rows = df.shape[0]

        num_iter_actual = num_rows
        num_iter_expected = self.max_iter + 2 # +2 is for saving iteration 0, and the final after the loop ends.
        assert num_iter_actual == num_iter_expected

    @classmethod
    def teardown_class(cls):
        """A function that is run ONCE after all unit tests in this class."""
        # Reset module level database
        import rmgpy.data.rmg

        rmgpy.data.rmg.database = None

        # Remove output directory
        shutil.rmtree(cls.outputDir)


class TestProfiling:
    @classmethod
    def setup_class(cls):
        """A function that is run ONCE before all unit tests in this class."""
        # Making the profile graph requires a display. See if one is available first
        cls.display_found = False

        try:
            cls.display_found = bool(os.environ["DISPLAY"])
        except KeyError:  # This means that no display was found
            pass
        cls.test_dir = os.path.join(originalPath, "..", "test", "rmgpy", "test_data", "mainTest")

    @patch("rmgpy.rmg.main.logging")
    def test_make_profile_graph(self, mock_logging):
        """
        Test that `make_profile_graph` function behaves properly given the current display state
        """
        profile_file = os.path.join(self.test_dir, "RMG.profile")
        make_profile_graph(profile_file)
        if self.display_found:  # Check that the profile graph was made
            assert os.path.exists(os.path.join(self.test_dir, "RMG.profile.dot.pdf"))
        else:  # We can't test making a profile graph on this system, but at least test that this was recognized
            mock_logging.warning.assert_called_with(
                "Could not find a display, which is required in order to generate "
                "the profile graph. This "
                "is likely due to this job being run on a remote server without performing X11 forwarding "
                "or running the job through a job manager like SLURM.\n\n The graph can be generated later "
                "by running with the postprocessing flag `rmg.py -P input.py` from any directory/computer "
                "where both the input file and RMG.profile file are located and a display is available.\n\n"
                "Note that if the postprocessing flag is specified, this will force the graph generation "
                "regardless of if a display was found, which could cause this program to crash or freeze."
            )

    @classmethod
    def teardown_class(cls):
        """A function that is run ONCE after all unit tests in this class."""

        if cls.display_found:  # Remove output PDF
            os.remove(os.path.join(cls.test_dir, "RMG.profile.dot.pdf"))
            os.remove(os.path.join(cls.test_dir, "RMG.profile.dot"))
            os.remove(os.path.join(cls.test_dir, "RMG.profile.dot.ps2"))


class TestCanteraOutputConversion:
    """
    Tests if we can convert Chemkin files to Cantera files without crashing.
    (Or raising an exception for bad files.)
    """
    def setup_class(self):
        self.chemkin_files = {
            """ELEMENTS
	H
	D /2.014/
	T /3.016/
	C
	CI /13.003/
	O
	OI /18.000/
	N

END

SPECIES
    ethane(1)       
    CH3(4)          
END

THERM ALL
   300.000  1000.000  5000.000

ethane(1)               H 6  C 2            G100.000   5000.000  954.52        1
 4.58987205E+00 1.41507042E-02-4.75958084E-06 8.60284590E-10-6.21708569E-14    2
-1.27217823E+04-3.61762003E+00 3.78032308E+00-3.24248354E-03 5.52375224E-05    3
-6.38573917E-08 2.28633835E-11-1.16203404E+04 5.21037799E+00                   4

CH3(4)                  H 3  C 1            G100.000   5000.000  1337.62       1
 3.54144859E+00 4.76788187E-03-1.82149144E-06 3.28878182E-10-2.22546856E-14    2
 1.62239622E+04 1.66040083E+00 3.91546822E+00 1.84153688E-03 3.48743616E-06    3
-3.32749553E-09 8.49963443E-13 1.62856393E+04 3.51739246E-01                   4

END



REACTIONS    KCAL/MOLE   MOLES

CH3(4)+CH3(4)=ethane(1)                             8.260e+17 -1.400    1.000    

END
""": True,
            """ELEMENTS
	CI /13.003/
	O
	OI /18.000/
	N

END

SPECIES
    ethane(1)       
    CH3(4)          
END

THERM ALL
   300.000  1000.000  5000.000

ethane(1)               H 6  C  2            G100.000   5000.000  954.52        1
 4.58987205E+00 1.41507042E-02-4.75958084E-06 8.60284590E-10-6.21708569E-14    2
-1.27217823E+04-3.61762003E+00 3.78032308E+00-3.24248354E-03 5.52375224E-05    3
-6.38573917E-08 2.28633835E-11-1.16203404E+04 5.21037799E+00                   4

CH3(4)                  H 3  C 1            G100.000   5000.000  1337.62       1
 3.54144859E+00 4.76788187E-03-1.82149144E-06 3.28878182E-10-2.22546856E-14    2
 1.62239622E+04 1.66040083E+00 3.91546822E+00 1.84153688E-03 3.48743616E-06    3
-3.32749553E-09 8.49963443E-13 1.62856393E+04 3.51739246E-01                   4

END



REACTIONS    KCAL/MOLE   MOLES

CH3(4)+CH3(4)=ethane(1)                             8.260e+17 -1.400    1.000    

END
""": False,
            """ELEMENTS
	H
	D /2.014/
	T /3.016/
	C
	CI /13.003/
	O
	OI /18.000/
	N

END

SPECIES
    ethane(1)       
    CH3(4)          
END

THERM ALL
   300.000  1000.000  5000.000

ethane(1)               H 6  C 2            G100.000   5000.000  954.52        1
 4.58987205E+00 1.41507042E-02-4.75958084E-06 8.60284590E-10-6.21708569E-14    2
-1.27217823E+04-3.61762003E+00 3.78032308E+00-3.24248354E-03 5.52375224E-05    3
-6.38573917E-08 2.28633835E-11-1.16203404E+04 5.21037799E+00                   4

END

REACTIONS    KCAL/MOLE   MOLES

CH3(4)+CH3(4)=ethane(1)                             8.260e+17 -1.400    1.000    

END
""": False,
        }
        self.rmg = RMG()
        self.dir_name = "temp_dir_for_testing"
        self.rmg.output_directory = os.path.join(originalPath, "..", "test", "rmgpy", "test_data", self.dir_name)

        self.tran_dat = """
! Species         Shape    LJ-depth  LJ-diam   DiplMom   Polzblty  RotRelaxNum Data     
! Name            Index    epsilon/k_B sigma     mu        alpha     Zrot      Source   
ethane(1)           2     252.301     4.302     0.000     0.000     1.500    ! GRI-Mech
CH3(4)              2     144.001     3.800     0.000     0.000     0.000    ! GRI-Mech
        """

    def teardown_class(self):
        os.chdir(originalPath)
        # try to remove the tree. If testChemkinToCanteraConversion properly
        # ran, the files should already be removed.
        try:
            shutil.rmtree(self.dir_name)
        except OSError:
            pass
        # go back to the main RMG-Py directory
        os.chdir("..")

    def test_chemkin_to_cantera_conversion(self):
        """
        Tests that good and bad chemkin files raise proper exceptions
        """

        from cantera.ck2yaml import InputError

        for ck_input, works in self.chemkin_files.items():
            os.chdir(originalPath)
            os.mkdir(self.dir_name)
            os.chdir(self.dir_name)

            f = open("chem001.inp", "w")
            f.write(ck_input)
            f.close()

            f = open("tran.dat", "w")
            f.write(self.tran_dat)
            f.close()

            if works:
                self.rmg.generate_cantera_files_from_chemkin(os.path.join(os.getcwd(), "chem001.inp"))
            else:
                with pytest.raises(InputError):
                    self.rmg.generate_cantera_files_from_chemkin(os.path.join(os.getcwd(), "chem001.inp"))

            # clean up
            os.chdir(originalPath)
            shutil.rmtree(self.dir_name)
