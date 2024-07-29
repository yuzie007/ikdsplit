"""Tests for `wycksplit.toml`."""

import collections
import os
import pathlib
import tomllib

import numpy as np
import pytest
from ase.spacegroup import Spacegroup

import ikdsplit
from ikdsplit.spacegroup import (
    check_equal,
    convert_symmetry_operations,
    get_setting_for_origin_choice_2,
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
    symops_sup = sg_sup.get_symop()
    sg_sub = Spacegroup(nsub, setting=get_setting_for_origin_choice_2(nsub))
    symops_sub = sg_sub.get_symop()
    transformation = invert(d["basis_change"], d["origin_shift"])
    symops_sub = convert_symmetry_operations(symops_sub, transformation)
    return symops_sup, symops_sub, sg_sup.subtrans


def reduce_subtranslations(
    subtranslations: np.ndarray,
    coset_representatives: list[dict[str, np.ndarray]],
) -> np.ndarray:
    """Reduce subtranslations according to coset representatives."""
    is_found = []
    rotation = np.eye(3, dtype=int)
    for subtranslation in subtranslations:
        for representative in coset_representatives:
            _ = representative["rotation"], representative["translation"]
            if check_equal(_, (rotation, subtranslation)):
                is_found.append(True)
                break
        else:
            is_found.append(False)
    if all(is_found):
        return np.zeros((1, 3))
    return subtranslations


@pytest.mark.parametrize("fn", toml_files)
def test_subgroup(fn: str) -> None:
    """Test subgroup."""
    nsup, nsub = (int(_) for _ in fn.split(".")[0].split("_"))
    with pathlib.Path(src / fn).open("rb") as f:
        d = tomllib.load(f)
    symops_sup, symops_sub, _ = get_symops(nsup, nsub, d)

    is_supergroup = []
    for symop_sub in symops_sub:
        for symop_sup in symops_sup:
            if check_equal(symop_sub, symop_sup):
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
    symops_sup, symops_sub, subtranslations = get_symops(nsup, nsub, d)

    subtranslations = reduce_subtranslations(
        subtranslations,
        d["coset_representatives"],
    )

    is_subgroup = []
    for symop_sup_orig in symops_sup:
        found = False
        for subtranslation in subtranslations:
            op = (np.eye(3, dtype=int), subtranslation)
            symop_sup = multiply(op, symop_sup_orig)
            for symop_sub in symops_sub:
                for i, _ in enumerate(d["coset_representatives"]):
                    coset_representative = _["rotation"], _["translation"]
                    tmp = multiply(coset_representative, symop_sub)
                    if check_equal(tmp, symop_sup):
                        is_subgroup.append(i)
                        found = True
                        break
                if found:
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
