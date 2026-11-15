"""Track trkB_engfix shared helpers. Paths are relative to PME (env PME, default /home/aaron/perseus_me) so the same code
runs on firebreather and in the Cardinal bundle. DDM keyword block (PREREG section 1)."""
import os, json, hashlib
from pathlib import Path
PME = Path(os.environ.get("PME", "/home/aaron/perseus_me"))
KF = Path(os.environ.get("PME_KF", str(Path.home() / "mmframe/fvs_modern_kf")))
DCLS = [(0.0, 3.0), (3.0, 5.0), (5.0, 7.0), (7.0, 9.0), (9.0, 12.0), (12.0, 15.0), (15.0, 999.0)]


def md5(p):
    return hashlib.md5(Path(p).read_bytes()).hexdigest()


def compute_block(date="1"):
    return f"COMPUTE   {date:>10}\nXBSDI = BSDI\nXBSMX = BSDIMAX\nEND\nDATABASE\nCOMPUTDB\nEND"


def ddm_block(a, b, ncyc):
    """FIXMORT replace mode, all species, one record per DBH class and cycle; Event Monitor rate
    1 - (1 + EXP(a_k + b*BSDI/BSDIMAX))**(YEAR - CENDYEAR - 1) = 1 - (1 - p)^(cycle length)"""
    assert len(a) == len(DCLS)
    out = []
    for c in range(1, ncyc + 1):
        for (lo, hi), ak in zip(DCLS, a):
            out.append(f"FIXMORT   {c:>10d}     Parms(0, &\n1-(1+EXP({ak:.6f}+{b:.6f}*BSDI/BSDIMAX))**(YEAR-CENDYEAR-1), &\n{lo:.1f}, {hi:.1f}, 0, 0)")
    return "\n".join(out)
