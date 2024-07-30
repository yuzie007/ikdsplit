"""Spacegroup."""

import re

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


def convert_symmetry_operations(
    old: list[tuple[np.ndarray, np.ndarray]],
    transformation: tuple[np.ndarray, np.ndarray],
) -> list[tuple[np.ndarray, np.ndarray]]:
    """Convert symmetry operations by the given transformation."""
    inverse = invert(transformation[0], transformation[1])
    new = [multiply(inverse, multiply(_, transformation)) for _ in old]
    return [(r, t - np.floor(t)) for r, t in new]


def check_equal(
    op0: tuple[np.ndarray, np.ndarray],
    op1: tuple[np.ndarray, np.ndarray],
) -> bool:
    """Check if the symmetry operations are equivalent."""
    d = op0[1] - op1[1]
    d -= np.rint(d)
    return np.allclose(op0[0], op1[0]) and np.allclose(d, 0.0)


def parse_transformation(s: str) -> tuple[np.ndarray, np.ndarray]:
    """Parse transformation string.

    Parameters
    ----------
    s : str
        String like "(x+y+1/2,y,z)".

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Rotation and translation.

    """

    def split_equation(equation: str) -> list[str]:
        """Split equation at "+" and "-"."""
        # https://docs.python.org/3/library/re.html#regular-expression-syntax
        # `(?<!^)` avoids empty string in "+x" -> ["", "+x"]
        return re.split(r"(?<!^)(?=[-+])", equation)

    def parse_coefficient(sr: str) -> int:
        """Parse coefficient preceding "xyz"."""
        if not sr:
            return 1
        if sr == "+":
            return 1
        if sr == "-":
            return -1
        return int(sr)

    def parse_translation(st: str) -> float:
        """Parse translation."""
        frac = st.split("/")
        numerator = int(frac[0])
        denominator = 1 if len(frac) == 1 else int(frac[1])
        return numerator / denominator

    p = s.strip("()").split(",")
    r = np.zeros((3, 3), dtype=int)
    t = np.zeros(3, dtype=float)
    xyz = "xyz"

    for i in range(3):
        for q in split_equation(p[i]):
            if q[-1] in xyz:
                j = xyz.index(q[-1])
                r[i, j] += parse_coefficient(q[:-1])
            else:
                t[i] += parse_translation(q)

    return r, t
