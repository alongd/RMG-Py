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

"""Input-file and restart seams for explicit EN wall and geometry selection."""
from pathlib import Path

import numpy as np
import pytest
import rmgpy
from rmgpy.exceptions import InputError, PlasmaStateError
from rmgpy.rmg.main import RMG
from rmgpy.rmg.input import read_input_file, save_input_file


@pytest.fixture(autouse=True)
def tested_tree():
    assert Path.cwd() in Path(rmgpy.__file__).resolve().parents


PREAMBLE = """
database(thermoLibraries=[],reactionLibraries=[],seedMechanisms=[],kineticsFamilies="none")
species(label='e-',structure=adjacencyList('1 e u1 p0 c-1'))
species(label='Ar',structure=adjacencyList('1 Ar u0 p4 c0'))
species(label='Arp',structure=adjacencyList('multiplicity 2\\n1 Ar u1 p3 c+1'))
species(label='Cl-',structure=adjacencyList('1 Cl u0 p4 c-1'))
"""


def input_text(model='confinedAnion',arm='fullFrequency',extra=''):
    return PREAMBLE + """
plasmaReactor(temperature=(298.15,'K'),pressure=(5,'torr'),electronTemperature=(35000,'K'),
    initialMoleFractions={'Ar':0.999998,'Arp':1.e-6,'e-':1.e-6,'Cl-':0.},
    chamberGeometry={'shape':'cylinder','radius':(5,'cm'),'length':(30,'cm')},
    ionReducedMobilities={'Arp':(1.535e-4,'m^2/(V*s)')},
    electronegativeWallModel=%r,electronegativeWallGeometry=%r,
    anionReducedMobilities={'Cl-':(1.5e-4,'m^2/(V*s)')},
    terminationTime=(1.e-5,'s'),%s)
simulator(atol=1.e-16,rtol=1.e-8)
model(toleranceMoveToCore=0.1,toleranceInterruptSimulation=0.1)
""" % (model,arm,extra)


@pytest.mark.parametrize('model',['confinedAnion','electropositiveBracket'])
@pytest.mark.parametrize('arm',['fullFrequency','radialOnly'])
def test_parse_save_and_reread_preserve_explicit_closure_and_geometry(tmp_path,model,arm):
    original=tmp_path/'input.py';original.write_text(input_text(model,arm))
    first=RMG();read_input_file(str(original),first)
    before=first.reaction_systems[0]
    assert before.electronegative_wall_model == model
    assert before.electronegative_wall_geometry == arm
    assert before.wall_diffusion_components == ((2.405/.05)**2,(np.pi/.30)**2)
    assert before.wall_chamber_geometry == dict(shape='cylinder', radius=.05, length=.30)
    saved=tmp_path/'saved.py';save_input_file(str(saved),first)
    second=RMG();read_input_file(str(saved),second)
    after=second.reaction_systems[0]
    assert after.electronegative_wall_model == model
    assert after.electronegative_wall_geometry == arm
    assert after.wall_diffusion_components == before.wall_diffusion_components
    assert after.wall_chamber_geometry == before.wall_chamber_geometry
    assert after.anion_reduced_mobilities == before.anion_reduced_mobilities


def test_full_frequency_is_the_input_default(tmp_path):
    text=input_text().replace(",electronegativeWallGeometry='fullFrequency'",'')
    path=tmp_path/'input.py';path.write_text(text)
    job=RMG();read_input_file(str(path),job)
    assert job.reaction_systems[0].electronegative_wall_geometry == 'fullFrequency'


@pytest.mark.parametrize('model,arm',[(None,'fullFrequency'),('arbitrary','fullFrequency'),('confinedAnion','arbitrary')])
def test_no_inferred_closure_or_arbitrary_selector(tmp_path,model,arm):
    path=tmp_path/'input.py';path.write_text(input_text(model,arm))
    with pytest.raises(PlasmaStateError):
        read_input_file(str(path),RMG())


def test_wrong_units_for_explicit_components_refuse(tmp_path):
    extra="wallDiffusionComponents={'radial':(1.,'m'),'axial':(1.,'m^-2')},"
    path=tmp_path/'input.py';path.write_text(input_text(extra=extra))
    with pytest.raises(InputError,match='dimensions'):
        read_input_file(str(path),RMG())
