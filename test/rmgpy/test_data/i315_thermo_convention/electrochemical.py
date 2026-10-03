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

name = 'convention fixture'
thermoConvention = 'electrochemical'
entry(
    index=1, label='proton', molecule='1 H u0 p0 c+1',
    thermo=NASA(polynomials=[NASAPolynomial(
        coeffs=[2.5, 0, 0, 0, 0, -745.375, 0],
        Tmin=(200, 'K'), Tmax=(6000, 'K'),
    )], Tmin=(200, 'K'), Tmax=(6000, 'K')),
)
entry(
    index=2, label='oxide', molecule='multiplicity 2\n1 O u1 p3 c-1',
    thermo=NASA(polynomials=[NASAPolynomial(
        coeffs=[2.5, 0, 0, 0, 0, -745.375, 0],
        Tmin=(200, 'K'), Tmax=(6000, 'K'),
    )], Tmin=(200, 'K'), Tmax=(6000, 'K')),
)
entry(
    index=3, label='hydrogen', molecule='multiplicity 2\n1 H u1 p0 c0',
    thermo=NASA(polynomials=[NASAPolynomial(
        coeffs=[2.5, 0, 0, 0, 0, -745.375, 0],
        Tmin=(200, 'K'), Tmax=(6000, 'K'),
    )], Tmin=(200, 'K'), Tmax=(6000, 'K')),
)
entry(
    index=4, label='oxygen', molecule='multiplicity 3\n1 O u2 p2 c0',
    thermo=NASA(polynomials=[NASAPolynomial(
        coeffs=[2.5, 0, 0, 0, 0, -745.375, 0],
        Tmin=(200, 'K'), Tmax=(6000, 'K'),
    )], Tmin=(200, 'K'), Tmax=(6000, 'K')),
)
