"""Tests for `wycksplit.toml`."""

import collections
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
    multiply,
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


def get_symops(nsup: int, nsub: int, d: dict) -> tuple:
    sg_sup = Spacegroup(nsup, setting=get_setting_for_origin_choice_2(nsup))
    symops_sup = get_symmetry_operations(sg_sup)
    sg_sub = Spacegroup(nsub, setting=get_setting_for_origin_choice_2(nsub))
    symops_sub = get_symmetry_operations(sg_sub)
    transformation = invert(d["basis_change"], d["origin_shift"])
    symops_sub = convert_symmetry_operations(symops_sub, transformation)
    return symops_sup, symops_sub, sg_sup.subtrans


@pytest.mark.parametrize("fn", toml_files)
def test_subgroup(fn: str) -> None:
    """Test subgroup."""
    nsup, nsub = (int(_) for _ in fn.split(".")[0].split("_"))
    with pathlib.Path(src / fn).open("rb") as f:
        d = tomllib.load(f)
    symops_sup, symops_sub, subtrans = get_symops(nsup, nsub, d)

    is_supergroup = []
    for symop_sub in symops_sub:
        for symop_sup in symops_sup:
            if check_equal(symop_sub, symop_sup, subtrans):
                is_supergroup.append(True)
                break
        else:
            is_supergroup.append(False)
    assert all(is_supergroup)


@pytest.mark.parametrize("fn", toml_files)
def test_coset_decomposition(fn: str) -> None:
    """Test coset decomposition."""
    nsup, nsub = (int(_) for _ in fn.split(".")[0].split("_"))
    with pathlib.Path(src / fn).open("rb") as f:
        d = tomllib.load(f)
    symops_sup, symops_sub, subtrans = get_symops(nsup, nsub, d)

    is_subgroup = []
    for symop_sup in symops_sup:
        found = False
        for symop_sub in symops_sub:
            for i, _ in enumerate(d["coset_representatives"]):
                coset_representative = _["rotation"], _["translation"]
                tmp = multiply(coset_representative, symop_sub)
                if check_equal(tmp, symop_sup, subtrans):
                    is_subgroup.append(i)
                    found = True
                    break
            if found:
                break
        else:
            is_subgroup.append(-1)

    counts = dict(collections.Counter(is_subgroup))
    print(counts)
    assert -1 not in counts
    assert len(counts) == len(d["coset_representatives"])
    assert all(_ == counts[0] for _ in counts.values()), counts
