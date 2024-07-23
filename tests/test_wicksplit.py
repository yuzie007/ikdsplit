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
        sup_ref, sub_ref = (int(_) for _ in fn.split(".")[0].split("_"))
        with pathlib.Path.open(src / fn, "rb") as f:
            d = tomllib.load(f)
        assert "coset_representatives" in d, fn
        assert d["space_group_number_sup"] == sup_ref, fn
        assert d["space_group_number_sub"] == sub_ref, fn
