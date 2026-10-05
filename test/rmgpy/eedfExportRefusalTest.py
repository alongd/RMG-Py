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

from types import SimpleNamespace

import pytest

from rmgpy.chemkin import ChemkinWriter, save_chemkin_file, write_kinetics_entry
from rmgpy.exceptions import EEDFExportError, MechanismWriterError
from rmgpy.kinetics import EEDFChannel
from rmgpy.reaction import Reaction
from rmgpy.rmg.model import ReactionModel
from rmgpy.species import Species
from rmgpy.thermo import NASA, NASAPolynomial
from rmgpy.yaml_cantera1 import CanteraWriter1, write_cantera
from rmgpy.yaml_cantera2 import CanteraWriter2, save_cantera_model


def _thermo():
    coefficients = [1.0, 0.0, 0.0, 0.0, 0.0, -100.0, 1.0]
    return NASA(
        polynomials=[
            NASAPolynomial(coeffs=coefficients, Tmin=(200, "K"), Tmax=(1000, "K")),
            NASAPolynomial(coeffs=coefficients, Tmin=(1000, "K"), Tmax=(6000, "K")),
        ],
        Tmin=(200, "K"),
        Tmax=(6000, "K"),
    )


def _eedf_model():
    h2 = Species(label="H2", index=1).from_smiles("[H][H]")
    h = Species(label="H", index=2).from_smiles("[H]")
    h2.thermo = _thermo()
    h.thermo = _thermo()
    reaction = Reaction(
        index=1,
        reactants=[h2],
        products=[h, h],
        reversible=False,
        kinetics=EEDFChannel("H2 -> H + H", "hydrogen-v1", "ine"),
    )
    return ReactionModel(species=[h2, h], reactions=[reaction]), reaction


def _assert_qualification_refusal(call):
    with pytest.raises(EEDFExportError, match="EEDF export .* requires qualification") as error:
        call()
    assert isinstance(error.value, MechanismWriterError)


def test_chemkin_reaction_serializer_refuses_eedf_kinetics():
    model, reaction = _eedf_model()

    _assert_qualification_refusal(
        lambda: write_kinetics_entry(reaction, model.species, verbose=False)
    )


def test_chemkin_file_export_refuses_eedf_kinetics_without_landing_file(tmp_path):
    model, _ = _eedf_model()
    path = tmp_path / "chem.inp"

    _assert_qualification_refusal(
        lambda: save_chemkin_file(path, model.species, model.reactions, verbose=False)
    )

    assert not path.exists()


def test_cantera_writer1_refuses_eedf_kinetics_without_landing_file(tmp_path):
    model, _ = _eedf_model()
    path = tmp_path / "cantera1.yaml"

    _assert_qualification_refusal(
        lambda: write_cantera(
            model.species,
            model.reactions,
            elements_in_use=model.get_elements(),
            path=path,
        )
    )

    assert not path.exists()


def test_cantera_writer2_refuses_eedf_kinetics_without_landing_file(tmp_path):
    model, _ = _eedf_model()
    path = tmp_path / "cantera2.yaml"

    _assert_qualification_refusal(lambda: save_cantera_model(model, path))

    assert not path.exists()


@pytest.mark.parametrize('writer_class', [ChemkinWriter, CanteraWriter1, CanteraWriter2])
def test_writer_listener_refuses_eedf_reactor_mode_without_markers(
        tmp_path, writer_class):
    output = tmp_path / writer_class.__name__
    output.mkdir()
    job = SimpleNamespace(
        output_directory=str(output),
        reaction_systems=[SimpleNamespace(eedf_mode=True)],
        reaction_model=SimpleNamespace(),
    )
    writer = writer_class(str(output))

    _assert_qualification_refusal(lambda: writer.update(job))

    assert not any(path.is_file() for path in output.rglob('*'))
