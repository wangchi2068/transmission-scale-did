"""Download the CDC state mask mandate raw dataset (217 MB) and verify its SHA-256 checksum."""
import hashlib
import sys
import urllib.request
from pathlib import Path

DEST = Path(__file__).resolve().parent / "raw" / "cdc_masks" / "rows.csv"
URL = "https://data.cdc.gov/api/views/62d6-pm5i/rows.csv"
EXPECTED_SHA256 = "2354fc3169c01a5fff9297ae3a232fc2c441d4c66e191433b6dc6d86d22d70c0"
EXPECTED_SIZE = 217183558


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(1024 * 1024):
            h.update(chunk)
    return h.hexdigest()


def main():
    DEST.parent.mkdir(parents=True, exist_ok=True)
    if DEST.exists():
        size = DEST.stat().st_size
        print(f"File exists: {DEST} ({size} bytes)")
        if size == EXPECTED_SIZE:
            digest = sha256_file(DEST)
            if digest == EXPECTED_SHA256:
                print("SHA-256 verified successfully.")
                return 0
            print(f"Checksum mismatch: {digest} != {EXPECTED_SHA256}")
    print(f"Downloading from {URL} to {DEST}...")
    urllib.request.urlretrieve(URL, DEST)
    digest = sha256_file(DEST)
    if digest == EXPECTED_SHA256:
        print(f"Download complete and SHA-256 verified: {digest}")
        return 0
    print(f"ERROR: Checksum mismatch after download: {digest}", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main() if len(sys.argv) > 1 and sys.argv[1] == "--run" else 0)
