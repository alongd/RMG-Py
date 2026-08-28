# Evidence for the units-pairing repair

Raw output, unedited, from the commands named in each file. All produced in
`/home/alon/Code/RMG-Py-i145-units-laundering` against
`database.directory = ../RMG-database-i123c-audit/input` (the probe prints the resolved
path at its head, so the database each run actually loaded is in the artifact rather than
in this sentence).

`BEFORE` files were produced with the two source hunks reverted and the extensions
rebuilt; `AFTER` files with the repair in place and rebuilt. Nothing else differs between
the two sets.

| file | what it is |
|---|---|
| `deck-chem.inp` | the mechanism the lithium deck generates. Byte-identical before and after the repair — the repair touches the read-back path, not generation. |
| `deck-stdout.tail.log`, `deck-stderr.tail.log` | the deck run's own streams. It exits 1 by design, on a pre-existing Chemkin-to-Cantera translation failure that is a separate ticket. |
| `probe-units-pairing.BEFORE.log` | `probe_units_pairing.py` before: 18 pass, 8 fail, exit 1. Three round trips lose 10^-6, 10^-12, 10^-18. |
| `probe-units-pairing.AFTER.log` | the same probe after: 26 pass, 0 fail, exit 0. Rate constant identical at every trip. |
| `pytest.NEWTESTS.RED.log` | the new tests against the unrepaired build: 4 failed, 3 passed. The failures print the 10^6. |
| `pytest.NEWTESTS.GREEN.log` | the same tests against the repaired build: 7 passed. |
| `merge.BEFORE.stdout.log`, `merge.BEFORE.reactions.inp` | `scripts/mergeModels.py` on the two half-mechanisms before: exit 0, `1.490e+12` and `1.734e+08`. |
| `merge.AFTER.stdout.log`, `merge.AFTER.reactions.inp` | the same merge after: exit 0, `1.490e+18` and `1.734e+14`, i.e. the numbers the deck generated. |

Reproducing, from the repo root:

```bash
conda activate rmg_env
export PYTHONPATH=/home/alon/Code/RMG-Py-i145-units-laundering:$PYTHONPATH
python -c "import rmgpy; print(rmgpy.__file__)"     # must be under this worktree
make build                                          # never bare `make`

mkdir -p /tmp/i145-deck && cp docs/i123-integration/input.py /tmp/i145-deck/
python rmg.py /tmp/i145-deck/input.py \
  > >(tee /tmp/i145-deck/stdout.log) 2> >(tee /tmp/i145-deck/stderr.log >&2)   # exits 1, by design

python docs/i145-units-laundering/probe_units_pairing.py    /tmp/i145-deck    # exit 0
python docs/i145-units-laundering/make_half_mechanisms.py   /tmp/i145-deck
python scripts/mergeModels.py \
  --model1 /tmp/i145-deck/merge-half0/chem.inp /tmp/i145-deck/merge-half0/species_dictionary.txt \
  --model2 /tmp/i145-deck/merge-half1/chem.inp /tmp/i145-deck/merge-half1/species_dictionary.txt
```
