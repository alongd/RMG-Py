from cantera import ck2yaml
import pytest


from compare_yaml_outputs import CompareYaml
from rmgpy.chemkin import load_chemkin_file
from rmgpy.rmg.model import ReactionModel
from rmgpy.yaml_cantera2 import save_cantera_model


@pytest.fixture(scope="module")
def compare_manager(tmp_path_factory):
    """Compare the two YAML writers from committed Chemkin inputs."""
    data = "rmgpy/tools/data/various_kinetics"
    chemkin = f"{data}/chem_annotated.inp"
    dictionary = f"{data}/species_dictionary.txt"
    transport = f"{data}/tran.dat"
    output = tmp_path_factory.mktemp("yaml-writer-comparison")
    converted = output / "ck2yaml.yaml"
    direct = output / "cantera2.yaml"

    species, reactions = load_chemkin_file(
        chemkin, dictionary, transport_path=transport, use_chemkin_names=True
    )
    save_cantera_model(
        ReactionModel(species=species, reactions=reactions), str(direct)
    )
    ck2yaml.Parser().convert_mech(
        chemkin,
        transport_file=transport,
        out_name=str(converted),
        quiet=True,
        permissive=True,
    )
    return CompareYaml(str(converted), str(direct))

def test_compare_number_of_species(compare_manager):
    assert compare_manager.compare_species_count()


def test_compare_species_names(compare_manager):
    assert compare_manager.compare_species_names()


def test_compare_species_count_per_phase(compare_manager):
    assert compare_manager.compare_species_count_per_phase()


def test_compare_reactions(compare_manager):
    assert compare_manager.compare_reactions()
