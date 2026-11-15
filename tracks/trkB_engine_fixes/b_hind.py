#!/usr/bin/env python3
"""Track trkB_engfix: H1 hindcast FVS-NE chains through engines/fvs/run_hindcast.py (run_chain and main unchanged) with two
patches: the binary is --bin (A.binary) and the default arm's keyword block is the DDM block of --params (fold entry --fold),
one FIXMORT set per cycle of the chain; --params none runs the plain default. Host and cfg YAML are rewritten to the bundle.
usage: b_hind.py --bin FVSne --params P.json|none --fold F --chains-file FILE --restricted DIR --out DIR [--workers N]"""
import argparse, json, os, sys, tempfile
from pathlib import Path
import yaml
sys.path.insert(0, str(Path(__file__).resolve().parent))
import b_lib as BL
sys.path.insert(0, str(BL.PME / "code/engines/fvs"))
import run_hindcast as RH  # noqa: E402

ap = argparse.ArgumentParser(); ap.add_argument("--bin", required=True); ap.add_argument("--params", required=True)
ap.add_argument("--fold", default="all"); ap.add_argument("--chains-file", default=""); ap.add_argument("--restricted", required=True)
ap.add_argument("--out", required=True); ap.add_argument("--workers", type=int, default=4); ap.add_argument("--max-chains", type=int, default=0)
a = ap.parse_args()
P = None if a.params == "none" else json.loads(Path(a.params).read_text())[a.fold]
RH.A.binary = lambda H, arm, v: a.bin
_ak = RH.arm_keywords
RH.arm_keywords = (lambda arm, cal, rows, y0, lens: _ak(arm, cal, rows, y0, lens) if P is None
                   else (_ak(arm, cal, rows, y0, lens) + "\n" + BL.ddm_block(P["a"], P["b"], len(lens))).strip())
host = yaml.safe_load((BL.PME / "code/profiles/cardinal/host.yaml").read_text())
host["fvs"]["bin_dir"] = str(Path(a.bin).parent); host["fvs"]["keyfile_builder_dir"] = str(BL.KF / "calibration/python")
cfg = yaml.safe_load((BL.PME / "code/config/hindcast/perseus_me_hindcast.yaml").read_text())
td = Path(tempfile.mkdtemp()); (td / "host.yaml").write_text(yaml.safe_dump(host))
cfg["engines"]["fvs"]["host"] = str(td / "host.yaml"); (td / "cfg.yaml").write_text(yaml.safe_dump(cfg))
argv = ["run_hindcast.py", "--restricted", a.restricted, "--cfg", str(td / "cfg.yaml"), "--arms", "default", "--out", a.out, "--workers", str(a.workers)]
if a.chains_file:
    argv += ["--chains-file", a.chains_file]
if a.max_chains:
    argv += ["--max-chains", str(a.max_chains)]
sys.argv = argv
RH.main()
