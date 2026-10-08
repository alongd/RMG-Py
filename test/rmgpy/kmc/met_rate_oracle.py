"""RMG-direct activation-rate oracle for the expanded termination inventory."""

from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.molecule.molecule import Molecule
from rmgpy.reaction import Reaction
from rmgpy.species import Species
from rmgpy.kmc.barrier_e0 import FixedBBarrierE0Provider
from rmgpy.kmc.reference_thermo import SharedThermoAssignment
from rmgpy.kinetics.arrhenius import ArrheniusBM, ArrheniusEP


def independent_termination_rates(artifact, database, temperature):
    by_id = {record["event_id"]: record for record in artifact["records"]}
    species_cache = {}
    estimated_rates = {}
    rates = {}
    provenance = artifact.get("provenance", {})
    database_sha = provenance.get("rmg_database_sha")
    assignment = SharedThermoAssignment(database.thermo, database_sha)
    provider_data = provenance.get("barrier_e0_provider")
    barrier_provider = (
        FixedBBarrierE0Provider(provider_data["B_K"])
        if provider_data and provider_data.get("enabled")
        else None
    )

    def species(graph):
        if graph not in species_cache:
            result = Species(molecule=[Molecule().from_adjacency_list(graph)])
            result.generate_resonance_structures()
            species_cache[graph] = result
        return species_cache[graph]

    for record in artifact["records"]:
        reactants = [species(graph) for graph in record["reactant_graphs"]]
        products = [species(graph) for graph in record["product_graphs"]]
        radical_delta = sum(
            participant.molecule[0].get_radical_count() for participant in products
        ) - sum(
            participant.molecule[0].get_radical_count() for participant in reactants
        )
        if (
            len(reactants) != 2
            or radical_delta >= 0
            or record["family"] not in {"R_Recombination", "Disproportionation"}
            or record["status"] == "refused"
        ):
            continue
        if record["inventory_class"] == "R1:J_ring":
            kind = record["junction_ops"][0]["junction_kind"]
            rates[record["event_id"]] = (
                1.76793e10 if kind == "J_para" else 3.53586e10 / 2
            ) * temperature**-1.00291
            continue
        direct = (
            record
            if record["rate_source"]["kind"] == "RMG family estimate"
            else by_id[record["reverse_of"]]
        )
        if direct["event_id"] not in estimated_rates:
            direct_reactants = [
                species(graph).copy(deep=True) for graph in direct["reactant_graphs"]
            ]
            direct_products = [
                species(graph).copy(deep=True) for graph in direct["product_graphs"]
            ]
            is_forward = direct["rate_source"]["family_template_direction"] == "forward"
            reaction = TemplateReaction(
                reactants=direct_reactants if is_forward else direct_products,
                products=direct_products if is_forward else direct_reactants,
                family=direct["family"],
                template=direct["template"].split(";") if direct["template"] else [],
                degeneracy=direct["raw_path_degeneracy"],
                is_forward=True,
            )
            kinetics, _, _, estimated_forward = database.kinetics.families[
                direct["family"]
            ].get_kinetics(
                reaction,
                template_labels=reaction.template,
                degeneracy=reaction.degeneracy,
                return_all_kinetics=False,
            )
            assert bool(estimated_forward) == is_forward
            source_reaction = TemplateReaction(
                reactants=direct_reactants if is_forward else direct_products,
                products=direct_products if is_forward else direct_reactants,
                family=direct["family"],
                template=direct["template"].split(";") if direct["template"] else [],
                degeneracy=direct["raw_path_degeneracy"],
                is_forward=True,
            )
            assignment.assign_reaction(source_reaction)
            evaluated = kinetics
            if isinstance(kinetics, (ArrheniusBM, ArrheniusEP)):
                source_reaction.kinetics = kinetics
                if barrier_provider is not None:
                    source_reaction, _ = barrier_provider.prepare_reaction(source_reaction)
                source_reaction.fix_barrier_height()
                evaluated = source_reaction.kinetics
            estimated_rates[direct["event_id"]] = evaluated.get_rate_coefficient(
                temperature
            )
        rate = estimated_rates[direct["event_id"]]
        if direct is not record:
            reference = Reaction(
                reactants=[species(graph) for graph in direct["reactant_graphs"]],
                products=[species(graph) for graph in direct["product_graphs"]],
            )
            assignment.assign_reaction(reference)
            rate /= reference.get_equilibrium_constant(temperature, type="Kc")
        rates[record["event_id"]] = rate * record["ssa_multiplier"]
    return rates
