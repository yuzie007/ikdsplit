"""Tests for crystal classes."""

import pytest

from ikdsplit.spacegroup import find_crystal_class


@pytest.mark.parametrize("space_group_number", range(1, 231))
def test_crystal_class(space_group_number: int) -> None:
    """Test `find_crystal_class`."""
    print(space_group_number, find_crystal_class(space_group_number))
