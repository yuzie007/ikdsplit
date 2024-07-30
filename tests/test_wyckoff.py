"""Tests for `wyckoff.csv`."""

import os
import pathlib

import pandas as pd
import pytest

import ikdsplit
from ikdsplit.spacegroup import parse_transformation

src = pathlib.Path(ikdsplit.__file__).parent / "database" / "wyckoff"
csv_files = sorted([_ for _ in os.listdir(src) if _.endswith(".csv")])


@pytest.mark.parametrize("fn", csv_files)
def test_wyckoff(fn: str) -> None:
    """Test if Wyckoff positions can be converted into the operator."""
    wyckoffs = pd.read_csv(src / fn)
    for d in wyckoffs.to_dict(orient="records")[::-1]:
        parse_transformation(d["coordinates"])
