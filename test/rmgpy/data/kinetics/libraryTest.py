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
import os.path
import re
import shutil
from copy import deepcopy

import pytest


import rmgpy
from rmgpy import settings
from rmgpy.data.kinetics.database import KineticsDatabase
from rmgpy.data.kinetics.family import TemplateReaction
from rmgpy.data.base import DatabaseError
from rmgpy.data.kinetics.library import (KineticsLibrary, LibraryReaction,
                                         _template_from_rate_rule_comment)
from rmgpy.kinetics import Arrhenius, Troe, PDepArrhenius
from rmgpy.kinetics.model import PDepKineticsModel


@pytest.mark.parametrize('provenance', ['ordinary', 'textual', 'inline_textual', 'external',
                                       'both', 'textual_rate_rule', 'external_rate_rule'])
@pytest.mark.parametrize('source_index', [None, 0, 23])
@pytest.mark.parametrize('auto_generated', [False, True])
def test_source_entry_index_provenance(tmp_path, provenance, source_index, auto_generated):
    from rmgpy.data.base import Entry
    from rmgpy.data.kinetics.common import format_external_library_provenance

    descriptions = {
        'ordinary': 'An ordinary source entry.',
        'textual': 'Originally from reaction library: Toy',
        'inline_textual': 'Seed provenance: Originally from reaction library: Toy',
        'external': format_external_library_provenance('Toy', str(tmp_path)),
    }
    descriptions['both'] = descriptions['textual'] + '\n' + descriptions['external']
    description = descriptions[provenance.replace('_rate_rule', '')]
    if source_index is not None:
        description += '\nOriginal reaction library entry index: {0}'.format(source_index)
    if provenance.endswith('_rate_rule'):
        description += '\nEstimated using rate rule [X]\nfamily: H_Abstraction'
    from rmgpy.reaction import Reaction
    library = KineticsLibrary(label='Ordinary', auto_generated=auto_generated)
    library.entries = {7: Entry(index=7, item=Reaction(reactants=[], products=[]),
                               data=Arrhenius(), long_desc=description)}
    reaction = library.get_library_reactions()[0]
    expected = (('Toy' if provenance != 'inline_textual' else 'Seed provenance: Toy', source_index)
                if auto_generated and provenance != 'ordinary' and source_index is not None else None)
    if not auto_generated:
        expected = ('Ordinary', 7)
    assert reaction.source_identity == expected
    assert reaction.get_source_entry_index() == (None if expected is None else expected[1])
    # The attached decision cannot drift when the entry or prose is edited later.
    reaction.entry.index = 123
    reaction.entry.long_desc = 'Originally from reaction library: Wrong\nOriginal reaction library entry index: 99'
    assert reaction.source_identity == expected
    assert reaction.get_source_entry_index() == (None if expected is None else expected[1])


@pytest.mark.parametrize('indices', ['not-an-index', '0\nOriginal reaction library entry index: 23'])
def test_loader_refuses_invalid_source_indices(indices):
    from rmgpy.data.base import Entry
    from rmgpy.reaction import Reaction

    library = KineticsLibrary(label='seed', auto_generated=True)
    library.entries = {7: Entry(
        index=7, item=Reaction(reactants=[], products=[]), data=Arrhenius(),
        long_desc='Originally from reaction library: Toy\nOriginal reaction library entry index: ' + indices)}
    with pytest.raises(DatabaseError, match='source entry ind'):
        library.get_library_reactions()


@pytest.mark.parametrize('description', [
    'Historical prose\nOriginally from reaction library: Toy\nOriginal reaction library entry index: 23',
    'Originally from reaction library: Toy\n\nOriginal reaction library entry index: 23',
    'Originally from reaction library: Toy\nRate commentary\nOriginal reaction library entry index: 23',
])
def test_free_text_cannot_complete_legacy_seed_identity(description):
    library = KineticsLibrary(label='seed', auto_generated=True)
    library.load_entry(index=7, label='H2 => H + H', kinetics=Arrhenius(), longDesc=description)
    assert library.get_library_reactions()[0].source_identity is None


@pytest.mark.parametrize('identity', [None, ('Toy', 0)])
def test_explicit_seed_identity_ignores_description_and_rate_comments(identity):
    description = ('Originally from reaction library: Archived\n'
                   'Original reaction library entry index: 17\n'
                   'Original reaction library entry index: 23')
    library = KineticsLibrary(label='seed', auto_generated=True)
    library.load_entry(index=7, label='H2 => H + H', kinetics=Arrhenius(comment=description),
                       longDesc=description, sourceIdentity=identity,
                       sourceProvenance=(identity[0], None) if identity else None)
    reaction = library.get_library_reactions()[0]
    assert reaction.source_identity == identity
    assert reaction.library == (identity[0] if identity else None)
    assert getattr(reaction, 'library_source', None) is None


@pytest.mark.parametrize('identity', [('Toy', True), ('Toy', '23'), ('', 23), ('Toy',), 23])
def test_loader_refuses_invalid_explicit_seed_identity(identity):
    library = KineticsLibrary(label='seed', auto_generated=True)
    with pytest.raises(DatabaseError, match='source identity'):
        library.load_entry(index=7, label='H2 => H + H', kinetics=Arrhenius(), sourceIdentity=identity)


class TestLibrary:
    @classmethod
    def setup_class(cls):
        """
        A function run ONCE before all unit tests in this class.
        """
        # Set up a dummy database
        cls.database = KineticsDatabase()
        cls.database.load_libraries(
            os.path.join(
                settings["test_data.directory"],
                "testing_database",
                "kinetics",
                "libraries",
            ),
            libraries=None,
        )  # this loads all of them
        cls.libraries = cls.database.libraries

    def test_get_library_reactions(self):
        """
        test that get_library_reactions loads reactions correctly
        """
        lib_rxns = self.libraries["GRI-Mech3.0"].get_library_reactions()
        for rxn in lib_rxns:
            assert isinstance(rxn, LibraryReaction)

        lib_rxns = self.libraries["ethane-oxidation"].get_library_reactions()  # should have no direct library reactions
        for rxn in lib_rxns:
            if isinstance(rxn.kinetics, PDepKineticsModel):
                assert isinstance(rxn, LibraryReaction)  # can load pdep as networks yet so load as libraries
            else:
                assert isinstance(rxn, TemplateReaction)  # all reactions are template based

    def _load_rate_rule_description(self, long_desc):
        source = self.libraries["ethane-oxidation"]
        entry = deepcopy(next(iter(source.entries.values())))
        entry.long_desc = long_desc
        library = KineticsLibrary(label="seed", auto_generated=True)
        library.entries[entry.index] = entry
        [reaction] = library.get_library_reactions()
        return reaction

    def test_rate_rule_prose_without_template_raises_named_database_error(self):
        source = self.libraries["ethane-oxidation"]
        entry = deepcopy(next(iter(source.entries.values())))
        entry.long_desc = "This prose mentions a rate rule but has no template.\nfamily: H_Abstraction"
        library = KineticsLibrary(label="prose-library", auto_generated=True)
        library.entries[entry.index] = entry

        with pytest.raises(DatabaseError, match="prose-library.*entry"):
            library.get_library_reactions()

    def test_rate_rule_prose_before_real_template_line_is_ignored(self):
        source = self.libraries["ethane-oxidation"]
        entry = deepcopy(next(iter(source.entries.values())))
        entry.long_desc = (
            "A note about this rate rule has no template here.\n"
            "Estimated using an average for rate rule [H2;O2b]\n"
            "family: H_Abstraction"
        )
        library = KineticsLibrary(label="seed", auto_generated=True)
        library.entries[entry.index] = entry

        [reaction] = library.get_library_reactions()
        assert isinstance(reaction, TemplateReaction)
        assert reaction.family == "H_Abstraction"
        assert reaction.template == "H2;O2b"

    def test_legacy_bracketed_rate_rule_text_keeps_base_template(self):
        assert _template_from_rate_rule_comment(
            "This reaction matched rate rule [Legacy;Template]\nfamily: LegacyFamily"
        ) == "Legacy;Template"

    def test_family_on_rate_rule_line_keeps_base_elif_behavior(self):
        reaction = self._load_rate_rule_description(
            "Exact match found for rate rule [X] family: H_Abstraction"
        )
        assert isinstance(reaction, TemplateReaction)
        assert reaction.family == ""
        assert reaction.template == "X"

    def test_second_rate_rule_keyword_keeps_base_fragment_boundary(self):
        reaction = self._load_rate_rule_description(
            "Exact match found for rate rule [X]; source rate rule [Source;Source]\n"
            "family: H_Abstraction"
        )
        assert isinstance(reaction, TemplateReaction)
        assert reaction.family == "H_Abstraction"
        assert reaction.template == "X"

    def test_multiple_usable_rate_rule_lines_keep_last_template(self):
        reaction = self._load_rate_rule_description(
            "Estimated using template [Matched1] for rate rule [First]\n"
            "Estimated using template [Matched2] for rate rule [Last]\n"
            "family: H_Abstraction"
        )
        assert isinstance(reaction, TemplateReaction)
        assert reaction.template == "Last"

    def test_vertical_tab_keeps_base_line_splitting(self):
        reaction = self._load_rate_rule_description(
            "Exact match found for rate rule [X]\vNote rate rule [Source;Source]\n"
            "family: H_Abstraction"
        )
        assert isinstance(reaction, TemplateReaction)
        assert reaction.template == "X"

    def test_saved_seed_fixture_round_trips_template_family_and_template(self):
        seed_path = os.path.join(
            settings["test_data.directory"], "restartTest", "seed_no_filters", "seed", "reactions.py"
        )
        library = KineticsLibrary(label="seed", auto_generated=True)
        library.load(seed_path, self.database.local_context, self.database.global_context)
        expected = {}
        for index, entry in library.entries.items():
            if "rate rule" in entry.long_desc:
                template = ""
                family = ""
                for line in entry.long_desc.split("\n"):
                    line = line.strip()
                    if "rate rule" in line:
                        rule_part = line.split("rate rule")[1].strip()
                        template = re.split(r"\[(.*)\]", rule_part)[1]
                    elif "family:" in line:
                        family = line.split("family:")[1].strip()
                expected[index] = (family, template)

        reactions = library.get_library_reactions()
        observed = {reaction.entry.index: (reaction.family, reaction.template)
                    for reaction in reactions if isinstance(reaction, TemplateReaction)}
        assert observed == expected
        assert len(observed) > 0

    @pytest.mark.parametrize('identity', [None, ('Real', 0)])
    @pytest.mark.parametrize('channel', ['longDesc', 'shortDesc', 'kinetics'])
    @pytest.mark.parametrize('comment', [
        'Originally from reaction library: Forged\nOriginal reaction library entry index: 77',
        '[RMG external library provenance v1] historical prose',
        '[RMG external library provenance v1] {"label":"Forged","source":"/forged"} '
        '[/RMG external library provenance v1]',
    ])
    def test_dedicated_owner_source_ignores_text(self, identity, channel, comment):
        from rmgpy.kinetics import Arrhenius

        library = KineticsLibrary(label='restart', auto_generated=True)
        fields = {channel: comment} if channel != 'kinetics' else {}
        library.load_entry(index=1, label='H2 => H + H', sourceIdentity=identity,
                           sourceProvenance=('Real', '/real') if identity else None,
                           kinetics=Arrhenius(comment=comment if channel == 'kinetics' else ''),
                           **fields)
        [reaction] = library.get_library_reactions()
        assert reaction.source_identity == identity
        assert reaction.library == ('Real' if identity else None)
        assert getattr(reaction, 'library_source', None) == ('/real' if identity else None)

    @pytest.mark.parametrize('identity', [None, ('Real', 0)])
    def test_explicit_no_source_provenance_keeps_identity_and_ignores_text(self, identity):
        from rmgpy.kinetics import Arrhenius

        library = KineticsLibrary(label='restart', auto_generated=True)
        library.load_entry(index=77, label='H2 => H + H', sourceIdentity=identity,
                           sourceProvenance=None, kinetics=Arrhenius(),
                           longDesc='Originally from reaction library: Forged\n'
                                    'Original reaction library entry index: 88')
        [reaction] = library.get_library_reactions()
        assert reaction.source_identity == identity
        assert reaction.library is None
        assert getattr(reaction, 'library_source', None) is None

    def test_conflicting_dedicated_owner_source_refuses(self):
        from rmgpy.kinetics import Arrhenius

        library = KineticsLibrary(label='restart', auto_generated=True)
        library.load_entry(index=1, label='H2 => H + H', kinetics=Arrhenius(),
                           sourceIdentity=('Real', 0), sourceProvenance=('Other', '/other'))
        with pytest.raises(DatabaseError, match='Conflicting dedicated source provenance labels'):
            library.get_library_reactions()

    @pytest.mark.parametrize('provenance', [('Real', 'relative'), ('', None), ('Real',),
                                          {'label': 'Real'}, ('Real', 3)])
    def test_invalid_dedicated_owner_source_refuses(self, provenance):
        from rmgpy.kinetics import Arrhenius

        library = KineticsLibrary(label='restart', auto_generated=True)
        with pytest.raises(DatabaseError, match='Invalid source provenance'):
            library.load_entry(index=1, label='H2 => H + H', kinetics=Arrhenius(),
                               sourceIdentity=None, sourceProvenance=provenance)

    def test_save_library(self):
        """
        This tests the library.save method by writing a new temporary file and
        loading it and comparing the original and copied reactions
        """
        for library_name in ["ethane-oxidation", "surface-example"]:
            copy_path = os.path.join(
                settings["test_data.directory"],
                "testing_database",
                "kinetics",
                "libraries",
                library_name + "-copy",
            )
            if os.path.exists(copy_path):
                logging.warning(f"Removing existing directory {copy_path}.")
                shutil.rmtree(copy_path)
            os.makedirs(copy_path)
            try:
                self.libraries[library_name].save(os.path.join(copy_path, "reactions.py"))
                self.database.load_libraries(
                    os.path.join(
                        settings["test_data.directory"],
                        "testing_database",
                        "kinetics",
                        "libraries",
                    ),
                    libraries=None,
                )  # this loads all of them, including the new copy
                original_rxns = self.database.libraries[library_name].get_library_reactions()
                copy_rxns = self.database.libraries[library_name + "-copy"].get_library_reactions()
                assert len(original_rxns) == len(copy_rxns), library_name
                for index, (original_rxn, copy_rxn) in enumerate(zip(original_rxns, copy_rxns)):
                    context = f"{library_name}[{index}]"
                    assert original_rxn.is_isomorphic(copy_rxn, either_direction=False), context
                    assert original_rxn.kinetics.is_identical_to(copy_rxn.kinetics), f"{context} kinetics"
                    # Todo: we need to check more carefully for equivalence,
                    # eg. is the metal attribute saved for surface reactions?
            finally:
                shutil.rmtree(copy_path)

    def test_loading_external_kinetic_library(self, tmp_path):
        """This tests loading a kinetic library which is not in the RMG-database repo"""
        kinetic_lib_in_db_path = os.path.join(settings["database.directory"], "kinetics", "libraries", "NOx2018")
        kinetic_lib_in_test_dir_path = str(tmp_path / "copied_kinetic_lib")
        os.makedirs(kinetic_lib_in_test_dir_path)
        for file_name in ["reactions.py", "dictionary.txt"]:
            shutil.copyfile(
                src=os.path.join(kinetic_lib_in_db_path, file_name),
                dst=os.path.join(kinetic_lib_in_test_dir_path, file_name),
            )
        self.database.load_libraries(path="", libraries=[kinetic_lib_in_test_dir_path])
        assert "copied_kinetic_lib" in list(self.database.libraries.keys())
        shutil.rmtree(kinetic_lib_in_test_dir_path, ignore_errors=True)

    def test_generate_high_p_limit_kinetics(self):
        """
        Test that a :class:Arrhenius kinetics object representing the high pressure limit rate
        is returned from Troe/Lindmann/PDepArrhenius/Chebyshev kinetic classes
        """
        lib_rxns = self.libraries["lib_net"].get_library_reactions()
        for rxn in lib_rxns:
            assert rxn.network_kinetics is None
            logging.debug("Processing reaction {0}".format(rxn))
            success = rxn.generate_high_p_limit_kinetics()
            if (isinstance(rxn.kinetics, PDepArrhenius) and rxn.kinetics.pressures.value_si[-1] < 9000000) or not rxn.is_unimolecular():
                # generate_high_p_limit_kinetics() should return `False` if the reaction is not unimolecular
                # or if it is a PDepArrhenius or Chebyshev with Pmax < 90 bar
                assert not success
            else:
                assert success
                if isinstance(rxn.kinetics, Arrhenius):
                    # If the library reaction is already an Arrhenius expression, network_kinetics isn't generated
                    assert rxn.network_kinetics is None
                else:
                    assert isinstance(rxn.network_kinetics, Arrhenius)
                    if isinstance(rxn.kinetics, Troe):
                        # This block quantitative tests the "H + CH2 <=> CH3" reaction
                        # from the test library test_data/testing_database/kinetics/libraries/lib_net/reactions.py
                        # 1. Check that the T exponent in the modified Arrhenius (the "n") equals to 0
                        assert round(abs(rxn.network_kinetics.n.value_si - 0), 7) == 0
                        # 2. Check that the pre-exponential factor equals to 6e+8 m^3/(mol*s)
                        assert round(abs(int(rxn.network_kinetics.A.value_si) - 6e8), 7) == 0
