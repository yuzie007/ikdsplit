"""Tests for `run`."""

import pathlib
import shutil
import subprocess

from ikdsplit.utils import cd


def test_run(tmp_path) -> None:
    """Test `run`."""
    src = pathlib.Path(__file__).parent
    fns = ["atoms_conventional.csv", "cell.dat", "ikdsplit.toml", "FPOSCAR"]
    for _ in fns:
        shutil.copy2(src / _, tmp_path)
    with cd(tmp_path):
        assert subprocess.call(["ikdsplit", "run", "-l", "2"]) == 0
