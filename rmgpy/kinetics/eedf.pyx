# cython: embedsignature=True, cdivision=True

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

"""Kinetics markers whose rates are supplied by an EEDF evaluation context."""

from rmgpy.kinetics.model import KineticsModel


cdef class EEDFChannel(KineticsModel):
    """
    Identify a reaction rate supplied by an EEDF collision channel.

    ``process`` and ``collision_set`` identify the provider-owned channel, and
    ``side`` selects its inferior (``"ine"``) or superior (``"sup"``) rate.
    The marker deliberately carries no standalone rate law.
    """

    def __init__(self, process='', collision_set='', side='ine', Tmin=None, Tmax=None, Pmin=None, Pmax=None,
                 comment=''):
        KineticsModel.__init__(
            self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, comment=comment
        )
        self.process = process
        self.collision_set = collision_set
        self.side = side

    property uses_eedf:
        """Whether this kinetics marker requires an EEDF provider."""
        def __get__(self):
            return True

    property side:
        """The provider rate bound used by this channel: ``"ine"`` or ``"sup"``."""
        def __get__(self):
            return self._side

        def __set__(self, value):
            if value not in ("ine", "sup"):
                raise ValueError("EEDFChannel side must be 'ine' or 'sup'.")
            self._side = value

    def __repr__(self):
        """Return an evaluable representation of this marker."""
        string = "EEDFChannel(process={0!r}, collision_set={1!r}, side={2!r}".format(
            self.process, self.collision_set, self.side
        )
        if self.Tmin is not None:
            string += ', Tmin={0!r}'.format(self.Tmin)
        if self.Tmax is not None:
            string += ', Tmax={0!r}'.format(self.Tmax)
        if self.Pmin is not None:
            string += ', Pmin={0!r}'.format(self.Pmin)
        if self.Pmax is not None:
            string += ', Pmax={0!r}'.format(self.Pmax)
        if self.comment != '':
            string += ', comment="""{0}"""'.format(self.comment)
        return string + ')'

    def __reduce__(self):
        """Return the constructor state used to serialize this marker."""
        return EEDFChannel, (
            self.process, self.collision_set, self.side,
            self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.comment,
        )

    cpdef double get_rate_coefficient(self, double T, double P=0.0) except -1:
        """Refuse standalone evaluation; an EEDF provider must supply this rate."""
        raise NotImplementedError(
            "EEDFChannel has no standalone rate coefficient; evaluate it in an EEDF context."
        )

    cpdef bint is_identical_to(self, KineticsModel other_kinetics) except -2:
        """Return whether another kinetics object identifies the same EEDF channel."""
        if not isinstance(other_kinetics, EEDFChannel):
            return False
        if not KineticsModel.is_identical_to(self, other_kinetics):
            return False
        return (
            self.process == (<EEDFChannel>other_kinetics).process
            and self.collision_set == (<EEDFChannel>other_kinetics).collision_set
            and self.side == (<EEDFChannel>other_kinetics).side
        )
