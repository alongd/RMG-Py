"""Red tests (I-213, S0): the MoM producer and the kMC runner on this line
agree on one ``--mom`` contract.

A ``mom.json`` the producer writes (``rmgpy.kmc.mom_producer``, from a
hand-built MoM profile: no RMG run) goes through the runner's own command
line (``rmgpy.kmc.run.main``, what ``python -m rmgpy.kmc.run`` calls): the
header form is accepted and recorded, every quantity the producer exports is
one the runner reports for kMC (no MoM-only row), and a header bound to
another deck or to other conventions is refused.
"""

import csv
import json

from momProducerTest import make_setup, moving_run, produce

from rmgpy.kmc import mom_producer as mp
from rmgpy.kmc import run as kr

END, SNAPSHOTS = 40.0, 9


def produced(tmp_path):
    setup = make_setup(tmp_path)
    out = tmp_path / "mom"
    produce(setup, out)
    moving_run(out, setup["spec"])
    produce(setup, out, mode="convert")
    return setup, out / "mom.json"


def run_cli(setup, mom, out):
    return kr.main([str(setup["deck"]), str(setup["catalog"]), "--temperature", "700", "--end", str(END),
                    "--snapshots", str(SNAPSHOTS), "--seeds", "1", "--max-events", "200",
                    "--out", str(out), "--mom", str(mom)])


def comparison(out):
    with open(out / "comparison.csv", newline="") as f:
        return list(csv.DictReader(line for line in f if not line.startswith("#")))


def test_a_produced_mom_file_round_trips_through_the_runner_command_line(tmp_path):
    setup, mom = produced(tmp_path)
    doc = json.loads(mom.read_text())
    assert run_cli(setup, mom, tmp_path / "run") == 0
    manifest = json.loads((tmp_path / "run" / "run.json").read_text())
    assert manifest["mom"]["header"] == doc["header"]
    assert manifest["mom"]["times"] == len(doc["entries"]) == SNAPSHOTS
    rows = comparison(tmp_path / "run")
    exported = {f"{k}.{name}" if isinstance(v, dict) else k
                for e in doc["entries"] for k, v in e["quantities"].items()
                for name in (v if isinstance(v, dict) else [None])}
    compared = {r["quantity"] for r in rows if r["mom"] not in ("", None)}
    assert exported and exported <= compared
    # every exported quantity has a kMC side: no MoM-only row
    kmc_side = {r["quantity"] for r in rows if r["kmc_mean"] not in ("", None)}
    assert exported - kmc_side == set()


def test_the_whole_ensemble_keys_are_reported_by_both_sides(tmp_path):
    setup, mom = produced(tmp_path)
    assert run_cli(setup, mom, tmp_path / "run") == 0
    whole = {"whole_mu0_mol", "whole_mu1_mol", "whole_mu2_mol", "whole_Mn_g_mol", "whole_Mw_g_mol"}
    assert whole <= set(mp.QUANTITY_SOURCES)
    t0 = [r for r in comparison(tmp_path / "run") if float(r["time_s"]) == 0.0]
    both = {r["quantity"] for r in t0 if r["kmc_mean"] not in ("", None) and r["mom"] not in ("", None)}
    assert whole <= both


def test_a_mom_file_of_another_deck_is_refused_by_the_runner(tmp_path, capsys):
    setup, mom = produced(tmp_path)
    doc = json.loads(mom.read_text())
    doc["header"]["producer"]["deck_sha256"] = "0" * 64
    mom.write_text(json.dumps(doc))
    assert run_cli(setup, mom, tmp_path / "run") == 2
    assert "not the deck being run" in capsys.readouterr().err


def test_a_mom_file_whose_conventions_changed_is_refused_by_the_runner(tmp_path, capsys):
    setup, mom = produced(tmp_path)
    conventions = mom.parent / "conventions.json"
    conventions.write_text(conventions.read_text() + "\n")
    assert run_cli(setup, mom, tmp_path / "run") == 2
    assert "conventions digest" in capsys.readouterr().err
