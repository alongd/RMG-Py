# `plasma_tdep` — a two-reaction plasma mechanism, as RMG wrote it

`chem.inp` and `species_dictionary.txt` are verbatim output of a run of
`docs/i123-integration/input.py` (the lithium charge network: electron-impact ionisation
from `PlasmaElectronImpactIonization`, radiative recombination from
`PlasmaRadiativeRecombination`). Nothing in them is hand-authored — in particular the
rate parameters are RMG's own export of the two library rate laws, not values chosen
here.

They exist so `TestPlasmaChemkinRoundTrip` in `test/rmgpy/chemkinTest.py` can exercise the
Chemkin round trip for a `TDEP`-marked plasma rate law without a database, a deck run, or
a network connection. The property under test is that reading this file and writing it
again does not move the rate constant: it used to divide every second-order rate constant
by exactly 10^6 per trip, silently, under an unchanged units header.

If you regenerate these files, regenerate both together — the reader resolves `e-(1)`,
`Li(2)` and `[Lip](3)` out of the dictionary, and a mismatched pair fails as a missing
species rather than as a stale fixture.
