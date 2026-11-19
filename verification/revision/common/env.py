"""Path setup shared by the revision analyses (imports the paper's estimator code read-only)."""
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REV = HERE.parent                        # verification/revision
RUN = REV.parents[1]                     # paper/
CODE = RUN / "experiment" / "code"
APPV = RUN / "verification" / "application"
for p in (str(CODE), str(APPV), str(HERE)):
    if p not in sys.path:
        sys.path.insert(0, p)
os.environ.setdefault("OMP_NUM_THREADS", "1")


def git_sha():
    import subprocess
    try:
        return subprocess.run(["git", "-C", str(RUN), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    except Exception:  # noqa: BLE001
        return None


def clean(x):
    import numpy as np
    if isinstance(x, dict):
        return {str(k): clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, np.ndarray):
        return clean(x.tolist())
    if isinstance(x, (np.floating, float)):
        x = float(x)
        return x if np.isfinite(x) else None
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return x
