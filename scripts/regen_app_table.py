"""Regenerate the application table with the current generator and the full application run."""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiment" / "code"))
import make_tables as MT  # noqa: E402

src = Path(sys.argv[1])
app = json.loads(src.read_text(encoding="utf-8"))
out, _ = MT.table_application(app)
(ROOT / "tables" / "table_7_application.tex").write_text(out, encoding="utf-8")
print("written", len(out))
