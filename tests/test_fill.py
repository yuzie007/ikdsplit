"""Tests for `fill`."""

import pathlib
import shutil
import subprocess

import pytest

from ikdsplit.utils import cd


@pytest.mark.parametrize("data_path", ["227", "002"])
def test_fill(data_path, tmp_path) -> None:
    """Test `fill`."""
    src = pathlib.Path(__file__).parent / "testdata" / data_path
    fns = ["atoms_conventional.csv", "cell.dat", "ikdsplit.toml", "FPOSCAR"]
    for _ in fns:
        if (src / _).exists():
            shutil.copy2(src / _, tmp_path)
    with cd(tmp_path):
        assert subprocess.call(["ikdsplit", "fill"]) == 0
