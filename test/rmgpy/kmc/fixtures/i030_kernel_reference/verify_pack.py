#!/home/alon/anaconda3/envs/rmg_env/bin/python
"""Assemble measured table, or execute the exact pytest sketch from pack.md."""

import argparse
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent
MET_PATH = next(parent / "rmgpy/kmc/met.py" for parent in ROOT.parents
                if (parent / "rmgpy/kmc/met.py").is_file())
RUN_DIRECTORY = Path("/home/alon/runs/i030-met-kernel-reference/fixture-verification")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--assemble", action="store_true")
    parser.add_argument("--output", type=Path, default=RUN_DIRECTORY,
                        help="directory for temporary pytest work and reproduced results")
    args = parser.parse_args()
    pack_path = ROOT / "pack.md"
    pack = pack_path.read_text()
    if args.assemble:
        table = (ROOT / "reference/results/results.md").read_text().rstrip()
        marker = "MEASUREMENT_TABLE_PLACEHOLDER"
        if marker in pack:
            pack = pack.replace(marker, table)
        else:
            pattern = r"<!-- BEGIN REPRODUCED NUMBERS -->.*?<!-- END REPRODUCED NUMBERS -->"
            pack, count = re.subn(pattern, lambda match: table, pack, flags=re.S)
            assert count == 1
        pack_path.write_text(pack)
        print("ASSEMBLED: measured table embedded in pack.md", flush=True)
        return
    match = re.search(r"```python\n(def test_met_kernel_reference\(tmp_path\):.*?)(?:\n```)", pack, re.S)
    if not match:
        raise AssertionError("no executable test sketch found in pack")
    prelude = f'''import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import subprocess
import sys
sys.dont_write_bytecode = True
source = Path({str(MET_PATH)!r})
spec = importlib.util.spec_from_file_location("met_reference_target", source)
met_module = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = met_module
spec.loader.exec_module(met_module)
'''
    env = dict(os.environ, PYTEST_DISABLE_PLUGIN_AUTOLOAD="1", PYTHONDONTWRITEBYTECODE="1",
               MET_KERNEL_REFERENCE_PACK=str(ROOT))
    for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        env[key] = "1"
    args.output = args.output.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="met-kernel-reference-", dir=args.output) as temporary:
        test_path = Path(temporary) / "test_met_kernel_reference.py"
        test_path.write_text(prelude + "\n" + match.group(1) + "\n")
        print("SKETCH VERIFIER: executing the exact proposed pytest body; all twelve ensembles rerun.", flush=True)
        subprocess.run([sys.executable, "-m", "pytest", "-c", "/dev/null", "--rootdir", temporary,
                        "--basetemp", str(Path(temporary) / "pytest-work"),
                        "-p", "no:cacheprovider", "-q", "-s", str(test_path)],
                       cwd=ROOT, env=env, check=True)
        outputs = list((Path(temporary) / "pytest-work").rglob("results.json"))
        if len(outputs) != 1:
            raise AssertionError("expected one complete reference output")
        fresh = json.loads(outputs[0].read_text())
        saved = json.loads((ROOT / "reference/results/results.json").read_text())
        if fresh != saved:
            raise AssertionError("full fresh rerun differs from the frozen research report")
        rendered = outputs[0].with_name("results.md").read_text().rstrip()
        if rendered not in pack:
            raise AssertionError("freshly reproduced table does not match the numbers in pack.md")
        destination = args.output / "reproduced"
        destination.mkdir(exist_ok=True)
        (destination / "results.json").write_text(outputs[0].read_text())
        (destination / "results.md").write_text(rendered + "\n")
        print(f"VERIFIER PASS: all {len(fresh['simulation'])} first-passage runs reproduced exactly; every generated pack number matches.", flush=True)
    print("PACK VERIFIER PASS: reproduced table and proposed pytest assertions passed.", flush=True)


if __name__ == "__main__":
    main()
