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

"""Apply only the approved temporal-growth mixing statement to pinned source."""
import hashlib
import json
from pathlib import Path


def apply_growth_patch(root):
    root = Path(root)
    pin = json.loads((root / 'growth-mixing-provenance.json').read_text())
    patch = (root / 'growth-mixing.patch').read_bytes()
    if hashlib.sha256(patch).hexdigest() != pin['patch_sha256']:
        raise ValueError('growth mixing patch checksum mismatch')
    source = root / 'loki-src/source/ElectronKinetics.cpp'
    original = source.read_bytes()
    if hashlib.sha256(original).hexdigest() != pin['original_source_sha256']:
        raise ValueError('LoKI source does not match the approved patch base')
    old = b'CIEffNew = mixingParameter * CIEffNew + (1 - mixingParameter) * CIEffOld;'
    new = b'CIEffNew = CIEffOld + mixingParameter * (CIEffNew - CIEffOld);'
    start = original.index(b'    // CIEff is <nu_eff>/N')
    position = original.index(old, start)
    updated = original[:position] + new + original[position + len(old):]
    if hashlib.sha256(updated).hexdigest() != pin['patched_source_sha256']:
        raise ValueError('patched LoKI source checksum mismatch')
    source.write_bytes(updated)
    print('LOKI_GROWTH_PATCH_PASS ' + pin['patched_source_sha256'], flush=True)


if __name__ == '__main__':
    apply_growth_patch(Path(__file__).resolve().parent)
