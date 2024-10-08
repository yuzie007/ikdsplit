"""Tests for `run`."""

import pathlib
import shutil
import subprocess

import pytest

from ikdsplit.utils import cd


@pytest.mark.parametrize("data_path", ["002", "194", "227"])
def test_run(data_path: str, tmp_path) -> None:
    """Test `run`."""
    level = 1 if data_path == "194" else 2

    src = pathlib.Path(__file__).parent / "testdata" / data_path
    fns = ["atoms_conventional.csv", "cell.dat", "ikdsplit.toml", "FPOSCAR"]
    for _ in fns:
        if (src / _).exists():
            shutil.copy2(src / _, tmp_path)

    commands = [
        ["ikdsplit", "split", "-r"],
        ["ikdsplit", "fill", "-r", "-l", str(level)],
        ["ikdsplit", "regress", "-r", "-l", str(level)],
        ["ikdsplit", "sort", "-r", "-l", str(level)],
    ]
    with cd(tmp_path):
        for command in commands:
            assert subprocess.call(command) == 0
