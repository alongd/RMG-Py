"""Reference-state provider interface and named failure contracts."""

import pytest

from rmgpy.kmc.reference_thermo import GasPhaseRMGReferenceThermo, ThermoUnavailable
from rmgpy.molecule.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species
from rmgpy.thermo import ThermoData


def _real_recombination():
    return Reaction(
        reactants=[Species(molecule=[Molecule(smiles="[CH3]")]) for _ in range(2)],
        products=[Species(molecule=[Molecule(smiles="CC")])],
    )


def test_gas_reference_provider_labels_reference_and_caches_species_lookup():
    class Database:
        def __init__(self):
            self.calls = []

        def get_thermo_data(self, species):
            self.calls.append(species.molecule[0].to_smiles())
            return ThermoData(
                Tdata=([300, 400, 600, 800, 1000], "K"),
                Cpdata=([30, 30, 30, 30, 30], "J/(mol*K)"),
                H298=(0, "kJ/mol"),
                S298=(100, "J/(mol*K)"),
            )

    database = Database()
    provider = GasPhaseRMGReferenceThermo(database, "test-only-database-commit")
    reaction = _real_recombination()
    for species in reaction.reactants + reaction.products:
        species.thermo = ThermoData(
            Tdata=([300, 400, 600, 800, 1000], "K"),
            Cpdata=([30, 30, 30, 30, 30], "J/(mol*K)"),
            H298=(1000, "kJ/mol"),
            S298=(100, "J/(mol*K)"),
        )
    temperatures = (600.0, 700.0, 800.0)
    actual = provider.equilibrium_constants(reaction, temperatures)
    assert actual == [
        reaction.get_equilibrium_constant(temperature, type="Kc")
        for temperature in temperatures
    ]
    assert provider.provenance["reference_thermo"] == "RMG gas-phase Kc"
    assert provider.provenance["rmg_database_sha"] == "test-only-database-commit"
    assert database.calls == ["[CH3]", "CC"]
    assert all(
        species.thermo.H298.value_si == 0
        for species in reaction.reactants + reaction.products
    )
    assert provider.equilibrium_constants(_real_recombination(), temperatures) == actual
    assert database.calls == ["[CH3]", "CC"]
    reordered = _real_recombination()
    reordered.products[0].molecule[0].atoms.reverse()
    assert provider.equilibrium_constants(reordered, temperatures) == actual
    assert database.calls == ["[CH3]", "CC"]


def test_missing_reference_thermo_names_the_unavailable_species():
    provider = GasPhaseRMGReferenceThermo(None, None)
    with pytest.raises(ThermoUnavailable) as failure:
        provider.equilibrium_constants(_real_recombination(), (700.0,))
    assert "thermo unavailable for species [CH3]: no thermo database" == str(
        failure.value
    )


def test_failed_rmg_lookup_preserves_named_species_and_evidence():
    class Database:
        def get_thermo_data(self, species):
            raise ValueError("induced lookup failure")

    provider = GasPhaseRMGReferenceThermo(Database(), "test-only-database-commit")
    with pytest.raises(ThermoUnavailable) as failure:
        provider.equilibrium_constants(_real_recombination(), (700.0,))
    assert (
        str(failure.value)
        == "thermo unavailable for species [CH3]: induced lookup failure"
    )
