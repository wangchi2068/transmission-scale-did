"""Targeted regression check for tolerance-based global-root deduplication."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

CODE_DIR = Path(__file__).resolve().parents[2] / "experiment" / "code"
sys.path.insert(0, str(CODE_DIR))

import qdgmm  # noqa: E402


def test_exact_and_near_duplicate_roots_keep_first_estimate():
    first = np.array([0.31, -0.12, 0.50])
    near = first + np.array([4e-9, -3e-9, 2e-9])
    distinct = first + np.array([2e-4, 0.0, 0.0])
    candidates = []

    assert qdgmm._append_unique_root(candidates, first)
    assert not qdgmm._append_unique_root(candidates, first.copy())
    assert not qdgmm._append_unique_root(candidates, near)
    assert qdgmm._append_unique_root(candidates, distinct)
    assert len(candidates) == 2
    assert np.array_equal(candidates[0][1], first)
    assert np.array_equal(candidates[1][1], distinct)


if __name__ == "__main__":
    test_exact_and_near_duplicate_roots_keep_first_estimate()
    print("root dedup check passed: exact/near duplicates collapsed; distinct root retained")
