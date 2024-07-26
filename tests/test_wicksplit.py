"""Tests for `wycksplit.toml`."""

import os
import pathlib
import tomllib

import pytest
from ase.spacegroup import Spacegroup

import ikdsplit
from ikdsplit.spacegroup import (
    check_equal,
    convert_symmetry_operations,
    get_setting_for_origin_choice_2,
    get_symmetry_operations,
    invert,
)

src = pathlib.Path(ikdsplit.__file__).parent / "database" / "wycksplit"
toml_files = sorted([_ for _ in os.listdir(src) if _.endswith(".toml")])


@pytest.mark.parametrize("fn", toml_files)
def test_wycksplit(fn: str) -> None:
    """Test `wycksplit.toml`."""
    sup_ref, sub_ref = (int(_) for _ in fn.split(".")[0].split("_"))
    with pathlib.Path(src / fn).open("rb") as f:
        d = tomllib.load(f)
    assert "coset_representatives" in d, fn
    assert d["space_group_number_sup"] == sup_ref, fn
    assert d["space_group_number_sub"] == sub_ref, fn


@pytest.mark.parametrize("fn", toml_files)
def test_subgroup(fn: str) -> None:
    """Test subgroup."""
    nsup, nsub = (int(_) for _ in fn.split(".")[0].split("_"))
    with pathlib.Path(src / fn).open("rb") as f:
        d = tomllib.load(f)
    sg_sup = Spacegroup(nsup, setting=get_setting_for_origin_choice_2(nsup))
    symops_sup = get_symmetry_operations(sg_sup)
    sg_sub = Spacegroup(nsub, setting=get_setting_for_origin_choice_2(nsub))
    symops_sub = get_symmetry_operations(sg_sub)
    transformation = invert(d["basis_change"], d["origin_shift"])
    symops_sub = convert_symmetry_operations(symops_sub, transformation)

    is_supergroup = []
    for symop_sub in symops_sub:
        for symop_sup in symops_sup:
            if check_equal(symop_sub, symop_sup, sg_sup.subtrans):
                is_supergroup.append(True)
                break
        else:
            is_supergroup.append(False)
    assert all(is_supergroup)
