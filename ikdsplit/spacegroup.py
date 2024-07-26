"""Spacegroup."""

import numpy as np


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
