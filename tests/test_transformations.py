"""Tests for `wycksplit.toml`."""

import numpy as np
import pytest
from ase.spacegroup import Spacegroup

from ikdsplit.spacegroup import (
    check_equal,
    convert_symmetry_operations,
    get_setting_for_origin_choice_2,
    get_transformations,
    invert,
    parse_transformation_csv,
)

transformations = get_transformations().to_dict(orient="records")


def get_symops(
    nsup: int,
    nsub: int,
    change_of_basis: np.ndarray,
    origin_shift: np.ndarray,
) -> tuple:
    sg_sup = Spacegroup(nsup, setting=get_setting_for_origin_choice_2(nsup))
    symops_sup = sg_sup.get_symop()
    sg_sub = Spacegroup(nsub, setting=get_setting_for_origin_choice_2(nsub))
    symops_sub = sg_sub.get_symop()
    transformation = invert(change_of_basis, origin_shift)
    symops_sub = convert_symmetry_operations(symops_sub, transformation)
    return symops_sup, symops_sub


@pytest.mark.parametrize("transformation", transformations)
def test_subgroup(transformation: dict) -> None:
    """Test subgroup."""
    change_of_basis, origin_shift = parse_transformation_csv(transformation)
    symops_sup, symops_sub = get_symops(
        transformation["supergroup"],
        transformation["subgroup"],
        change_of_basis,
        origin_shift,
    )
    is_supergroup = []
    for symop_sub in symops_sub:
        for symop_sup in symops_sup:
            if check_equal(symop_sub, symop_sup):
                is_supergroup.append(True)
                break
        else:
            is_supergroup.append(False)
    assert all(is_supergroup)
