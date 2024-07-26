"""Spacegroup."""

import numpy as np


def get_setting_for_origin_choice_2(space_group_number: int) -> int:
    """Get `setting` for origin choice 2."""
    orthorhombic = [48, 50, 59, 68, 70]
    tetragonal = [85, 86, 88, 125, 126, 129, 130, 133, 134, 137, 138, 141, 142]
    cubic = [201, 203, 222, 224, 227, 228]
    return 2 if space_group_number in orthorhombic + tetragonal + cubic else 1


def invert(
    basis_change: np.ndarray,
    origin_shift: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Invert a general change of the coordinate system.

    Parameters
    ----------
    basis_change : np.ndarray
        P
    origin_shift : np.ndarray
        p

    """
    inv_basis_change = np.linalg.inv(basis_change)
    inv_origin_shift = -1.0 * inv_basis_change @ origin_shift
    return inv_basis_change, inv_origin_shift


def multiply(
    op0: tuple[np.ndarray, np.ndarray],
    op1: tuple[np.ndarray, np.ndarray],
) -> tuple[np.ndarray, np.ndarray]:
    """Multiply two changes of the coordinate system.

    See Sec. (1.2.2.2) in ITA (2016).
    """
    basis_change_0, origin_shift_0 = op0
    basis_change_1, origin_shift_1 = op1
    basis_change = basis_change_0 @ basis_change_1
    origin_shift = basis_change_0 @ origin_shift_1 + origin_shift_0
    return basis_change, origin_shift
