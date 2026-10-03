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

name = 'explicit superelastic pair'
shortDesc = 'Independent irreversible rates for both directions.'
longDesc = ''

# The loader requires opposite-direction entries to be marked as duplicates;
# it only merges duplicate entries in the same direction.
entry(
    index=1,
    label='e- + N2v1 => e- + N2',
    reversible=False,
    duplicate=True,
    kinetics=TwoTemperaturePlasma(A=(1.0e3, 'm^3/(mol*s)'), n=0.0, Ea_g=(0.0, 'J/mol'), Ea_e=(0.0, 'J/mol')),
)
entry(
    index=2,
    label='e- + N2 => e- + N2v1',
    reversible=False,
    duplicate=True,
    kinetics=TwoTemperaturePlasma(A=(1.0e2, 'm^3/(mol*s)'), n=0.0,
                                  Ea_g=(0.0, 'J/mol'), Ea_e=(0.0, 'J/mol')),
)
