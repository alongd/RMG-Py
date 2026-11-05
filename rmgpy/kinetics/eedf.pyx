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
    A sourced cascade supplies ``branch_product``, ``branch_yield`` and
    ``yield_table_sha256`` together. These immutable fields select one product
    of a complete family; cascade branches can only use the inferior side.
    The marker deliberately carries no standalone rate law.
    ``role='rate_only'`` supplies forward dissociation chemistry covered by
    the complete N2 energy inventory; it cannot bind an energy-loss entry.
    """

    def __init__(self, process='', collision_set='', side='ine', Tmin=None, Tmax=None, Pmin=None, Pmax=None,
                 comment='', branch_product=None, branch_yield=None, yield_table_sha256=None,
                 role='collision'):
        # Validate before changing state, including explicit reinitialization
        # through native-record reconstruction. Ordinary markers need no source.
        from rmgpy.tools.eedf.ownership import validate_branch_metadata
        from rmgpy.tools.eedf.rate_only import validate_rate_only_marker
        validate_rate_only_marker(role, side, branch_product)
        validate_branch_metadata(branch_product, branch_yield, yield_table_sha256, side)
        KineticsModel.__init__(
            self, Tmin=Tmin, Tmax=Tmax, Pmin=Pmin, Pmax=Pmax, comment=comment
        )
        self.process = process
        self.collision_set = collision_set
        self._branch_product = branch_product
        self._branch_yield = branch_yield
        self._yield_table_sha256 = yield_table_sha256
        self._role = role
        self.side = side

    property role:
        """Collision energy owner or forward chemical-only rate owner."""
        def __get__(self):
            return self._role

    property branch_product:
        """Sourced product state, or None for an ordinary owner."""
        def __get__(self):
            return self._branch_product

    property branch_yield:
        """Sourced product fraction, including an explicit zero."""
        def __get__(self):
            return self._branch_yield

    property yield_table_sha256:
        """Digest of the generated yield source that owns this fraction."""
        def __get__(self):
            return self._yield_table_sha256

    property uses_eedf:
        """Whether this kinetics marker requires an EEDF provider."""
        def __get__(self):
            return True

    property side:
        """The provider rate bound used by this channel: ``"ine"`` or ``"sup"``."""
        def __get__(self):
            return self._side

        def __set__(self, value):
            if self._role == 'rate_only' and value != 'ine':
                raise ValueError('rate-only owner cannot use a superelastic side')
            if value not in ("ine", "sup"):
                raise ValueError("EEDFChannel side must be 'ine' or 'sup'.")
            if self._branch_product is not None and value != 'ine':
                raise ValueError('cascade branch cannot use the superelastic side')
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
        if self.branch_product is not None:
            string += ', branch_product={0!r}, branch_yield={1!r}, yield_table_sha256={2!r}'.format(
                self.branch_product, self.branch_yield, self.yield_table_sha256)
        if self.role != 'collision':
            string += ', role={0!r}'.format(self.role)
        return string + ')'

    def __reduce__(self):
        """Return the constructor state used to serialize this marker."""
        args = (
            self.process, self.collision_set, self.side,
            self.Tmin, self.Tmax, self.Pmin, self.Pmax, self.comment,
        )
        if self.branch_product is not None or self.role != 'collision':
            args += (self.branch_product, self.branch_yield, self.yield_table_sha256)
        if self.role != 'collision':
            args += (self.role,)
        return EEDFChannel, args

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
            and self.role == (<EEDFChannel>other_kinetics).role
            and self.branch_product == (<EEDFChannel>other_kinetics).branch_product
            and self.branch_yield == (<EEDFChannel>other_kinetics).branch_yield
            and self.yield_table_sha256 == (<EEDFChannel>other_kinetics).yield_table_sha256
        )
