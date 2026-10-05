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

import pytest

from rmgpy.exceptions import NonEquilibriumReverseRateError
from rmgpy.kinetics import EEDFChannel
from rmgpy.reaction import Reaction
from rmgpy.species import Species


def test_reversible_eedf_channel_refuses_reverse_from_equilibrium():
    ar = Species(label="Ar").from_adjacency_list("1 Ar u0 p4 c0")
    ar_excited = Species(label="Ar*").from_adjacency_list(
        "electronicstate Ar(4s)\n1 Ar u0 p4 c0"
    )
    reaction = Reaction(
        reactants=[ar],
        products=[ar_excited],
        reversible=True,
        kinetics=EEDFChannel("Ar -> Ar*", "argon-v1", "ine"),
    )

    reason = reaction.get_reverse_from_equilibrium_refusal()

    assert "EEDF provider" in reason
    assert "qualified reverse EEDF channel" in reason
    with pytest.raises(NonEquilibriumReverseRateError, match="EEDF provider"):
        reaction.check_reverse_from_equilibrium_supported()
