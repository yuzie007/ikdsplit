"""Tests for `wycksplit.toml`."""

import os
import pathlib
import tomllib

import ikdsplit


def test_wycksplit() -> None:
    """Test `wycksplit.toml`."""
    src = pathlib.Path(ikdsplit.__file__).parent / "database" / "wycksplit"
    fns = sorted([_ for _ in os.listdir(src) if _.endswith(".toml")])
    for fn in fns:
        print(fn)
        with pathlib.Path.open(src / fn, "rb") as f:
            d = tomllib.load(f)
        print(d)
