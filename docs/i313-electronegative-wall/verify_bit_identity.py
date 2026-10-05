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

"""Dump binary argon reference results for base/tip cmp.

Uses the established synthetic argon energy fixture, with species-specific
map transport. The baseline lacks EN inputs and refuses even zero-population
core anions. For zero-anion cases its pure-argon model carries an inert zero neutral
placeholder so both solvers have the same packed dimension. Compare the
complete numerical state and flux arrays; the inactive population and its
wall flux must remain exactly zero. This does not relax the base refusal.
"""
import argparse
import inspect
import pickle
import sys
from pathlib import Path

import numpy as np
import rmgpy
import rmgpy.data.rmg as data_module

ROOT = Path.cwd().resolve()
assert ROOT in Path(rmgpy.__file__).resolve().parents, rmgpy.__file__
sys.path.insert(0,str(ROOT/'test/rmgpy/solver'))
import plasmaEnergyBalanceTest as fixture
from rmgpy.species import Species


class MapReferenceReactor(fixture.ToyLibraryReactor):
    closure = None
    arm = 'fullFrequency'
    zero_anion = False

    def __init__(self,*args,**kw):
        scalar = kw.pop('ion_reduced_mobility')
        kw['ion_reduced_mobilities'] = {'Ar+':scalar}
        supported = hasattr(fixture.PlasmaReactor,'electronegative_wall_model')
        self.verify_anion = None
        if supported and self.closure is not None:
            kw.update(electronegative_wall_model=self.closure,
                      electronegative_wall_geometry=self.arm,
                      wall_diffusion_components=((2.405/fixture.RADIUS)**2,(np.pi/fixture.LENGTH)**2))
            if self.zero_anion:
                self.verify_anion=Species(label='Cl-').from_adjacency_list('1 Cl u0 p4 c-1')
                self.verify_anion.thermo=fixture._thermo(0.)
                args=list(args)
                args[2]=dict(args[2]);args[2][self.verify_anion]=0.
                kw['anion_reduced_mobilities']={'Cl-':(1.5e-4,'m^2/(V*s)')}
                kw['thermo_source_assertions']=dict(kw['thermo_source_assertions'], **{'Cl-': 'ion'})
        if not supported and self.zero_anion:
            # Equal-size packed system at the old API: its policy refuses a
            # core anion even at zero, so use a zero inert neutral placeholder.
            # This changes neither the argon state nor its sources, but holds
            # the solver's error-norm dimension fixed for the comparison.
            self.verify_anion=Species(label='Cl').from_adjacency_list('1 Cl u1 p3 c0')
            self.verify_anion.thermo=fixture._thermo(0.)
            args=list(args);args[2]=dict(args[2]);args[2][self.verify_anion]=0.
            kw['wall_single_bath_approximation']=True
            kw['electron_energy_balance']['elastic_collisions']['Cl']={'ignore':'zero inactive reference placeholder'}
        super().__init__(*args,**kw)

    def initialize_model(self,core,rxns,edge,edge_rxns,**kw):
        if self.verify_anion is not None:
            core.append(self.verify_anion)
        return super().initialize_model(core,rxns,edge,edge_rxns,**kw)


def dump_case(output,closure=None,arm='fullFrequency',zero_anion=False):
    data_module.database=None
    MapReferenceReactor.closure=closure
    MapReferenceReactor.arm=arm
    MapReferenceReactor.zero_anion=zero_anion
    reactor,core,rxns=fixture._build(cls=MapReferenceReactor,elastic=False)
    common=[reactor.species_index[s] for s in core if s.label in ('e-','Ar','Ar+')]
    trajectory=[]; wall=[]; frequency=[]; budget=[]
    for time in [0.,1.e-10,1.e-9,1.e-8,1.e-7,1.e-6]:
        if time>0.:
            reactor.advance(time)
        if reactor.verify_anion is not None:
            index=reactor.species_index[reactor.verify_anion]
            assert reactor.y[index] == 0.
            assert reactor.wall_flux[index] == 0.
        state=reactor.y.copy()
        trajectory.append((reactor.t,state.tobytes()))
        wall.append(reactor.wall_flux.tobytes())
        frequency.append(np.float64(reactor.nu_wall_latched).tobytes())
        budget.append(reactor.energy_budget.copy())
    output.mkdir(parents=True,exist_ok=True)
    for name,value in [('trajectory',trajectory),('wall_flux',wall),
                       ('nu_wall_latched',frequency),('energy_budget',budget)]:
        (output/(name+'.bin')).write_bytes(pickle.dumps(value,protocol=4))
    print('DUMP',output)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',required=True)
    args=parser.parse_args();output=Path(args.output)
    dump_case(output/'map_no_anions')
    for closure in ('confinedAnion','o2ReferenceQualifiedUnity'):
        for arm in ('fullFrequency','radialOnly'):
            dump_case(output/(closure+'_'+arm),closure,arm,True)


if __name__ == '__main__':
    main()
