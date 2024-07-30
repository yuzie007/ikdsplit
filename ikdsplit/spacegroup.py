"""Spacegroup."""

import collections
import pathlib
import re

import numpy as np
import pandas as pd
from ase.spacegroup import Spacegroup

import ikdsplit


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


def find_wyckoff(xyz0: np.ndarray, space_group_number: int) -> str:
    """Find Wyckoff position.

    Parameters
    ----------
    xyz0 : np.ndarray
        Position to be checked.
    space_group_number : int
        Space group number.

    Returns
    -------
    str
        Wyckoff letter with multiplicity.

    """
    src = pathlib.Path(ikdsplit.__file__).parent / "database" / "wyckoff"
    wyckoffs = pd.read_csv(src / f"{space_group_number:03d}.csv")

    # get all symmetrically equivalent positions
    # one of them should be equal to the representative Wyckoff coordinates
    xyzs = expand_equivalent_positions(xyz0, space_group_number)

    for d in wyckoffs.to_dict(orient="records")[::-1]:
        # skip if #(symmetrically equivalent points) != multiplicity
        if len(xyzs) != d["multiplicity"]:
            continue
        for xyz in xyzs:
            r, t = parse_transformation(d["coordinates"])
            xyz_refined = r @ xyz + t
            diff = xyz - xyz_refined
            diff -= np.rint(diff)
            if np.allclose(diff, 0.0):
                return str(d["multiplicity"]) + d["wyckoff_letter"]
    raise RuntimeError(xyz0, space_group_number)


def expand_equivalent_positions(
    xyz0: np.ndarray,
    space_group_number: int,
) -> np.ndarray:
    """Expand symmetrically equivalent positions.

    Parameters
    ----------
    xyz0 : np.ndarray
        Original positions.
    space_group_number : int
        Space group number.

    Returns
    -------
    xyzs : np.ndarray
        Symmetrically equivalent positions.

    """
    setting = get_setting_for_origin_choice_2(space_group_number)
    symops = Spacegroup(space_group_number, setting).get_symop()
    xyzs = [xyz0]
    for r, t in symops:
        xyz = r @ xyz0 + t
        diff = xyz[None, :] - np.array(xyzs)
        diff -= np.rint(diff)
        if not np.any(np.all(np.isclose(diff, 0.0), axis=1)):
            xyzs.append(xyz)
    return np.array(xyzs)


def reduce_equivalent_positions(
    xyz0s: np.ndarray,
    space_group_number: int,
) -> np.ndarray:
    """Reduce symmetrically equivalent positions.

    Parameters
    ----------
    xyz0s : np.ndarray
        Positions to be reduced.
    space_group_number : int
        Space group number.

    Returns
    -------
    np.ndarray
        Reduced positions.

    """
    setting = get_setting_for_origin_choice_2(space_group_number)
    symops = Spacegroup(space_group_number, setting).get_symop()
    xyzs = [xyz0s[0]]
    for xyz0 in xyz0s:
        for r, t in symops:
            xyz = r @ xyz0 + t
            diff = xyz[None, :] - np.array(xyzs)
            diff -= np.rint(diff)
            if np.any(np.all(np.isclose(diff, 0.0), axis=1)):
                break
        else:
            xyzs.append(xyz0)
    return np.array(xyzs)


def find_rotation_type(rotation: np.ndarray) -> int:
    """Find rotation type.

    References
    ----------
    - Sec. 1.2.2.4. in International Tables for Crystallography A (2016)
    - Table 1 in Grosse-Kunstleve, Acta Crystallogr A 55, 383 (1999).
    - Table V in Togo et al., arXiv:1808.01590 (2018).

    """
    # https://doi.org/10.1107/S0108767398010186
    tr = np.trace(rotation)
    det = int(round(np.linalg.det(rotation)))
    return {
        (-2, -1): -6,
        (-1, -1): -4,
        (00, -1): -3,
        (+1, -1): -2,
        (-3, -1): -1,
        (+3, +1): +1,
        (-1, +1): +2,
        (00, +1): +3,
        (+1, +1): +4,
        (+2, +1): +6,
    }[(tr, det)]


def find_crystal_class(space_group_number: int) -> str:
    """Find crystal class.

    References
    ----------
    - Table VI in Togo et al., arXiv:1808.01590 (2018).

    """
    sg = Spacegroup(space_group_number)
    rotation_types = [find_rotation_type(r) for r in sg.get_rotations()]
    counter = collections.Counter(rotation_types)
    key = tuple(counter[_] for _ in (-6, -4, -3, -2, -1, 1, 2, 3, 4, 6))
    return {
        (0, 0, 0, 0, 0, 1, 0, 0, 0, 0): "1",
        (0, 0, 0, 0, 1, 1, 0, 0, 0, 0): "-1",
        (0, 0, 0, 0, 0, 1, 1, 0, 0, 0): "2",
        (0, 0, 0, 1, 0, 1, 0, 0, 0, 0): "m",
        (0, 0, 0, 1, 1, 1, 1, 0, 0, 0): "2/m",
        (0, 0, 0, 0, 0, 1, 3, 0, 0, 0): "222",
        (0, 0, 0, 2, 0, 1, 1, 0, 0, 0): "mm2",
        (0, 0, 0, 3, 1, 1, 3, 0, 0, 0): "mmm",
        (0, 0, 0, 0, 0, 1, 1, 0, 2, 0): "4",
        (0, 2, 0, 0, 0, 1, 1, 0, 0, 0): "-4",
        (0, 2, 0, 1, 1, 1, 1, 0, 2, 0): "4/m",
        (0, 0, 0, 0, 0, 1, 5, 0, 2, 0): "422",
        (0, 0, 0, 4, 0, 1, 1, 0, 2, 0): "4mm",
        (0, 2, 0, 2, 0, 1, 3, 0, 0, 0): "-42m",
        (0, 2, 0, 5, 1, 1, 5, 0, 2, 0): "4/mmm",
        (0, 0, 0, 0, 0, 1, 0, 2, 0, 0): "3",
        (0, 0, 2, 0, 1, 1, 0, 2, 0, 0): "-3",
        (0, 0, 0, 0, 0, 1, 3, 2, 0, 0): "32",
        (0, 0, 0, 3, 0, 1, 0, 2, 0, 0): "3m",
        (0, 0, 2, 3, 1, 1, 3, 2, 0, 0): "-3m",
        (0, 0, 0, 0, 0, 1, 1, 2, 0, 2): "6",
        (2, 0, 0, 1, 0, 1, 0, 2, 0, 0): "-6",
        (2, 0, 2, 1, 1, 1, 1, 2, 0, 2): "-6/m",
        (0, 0, 0, 0, 0, 1, 7, 2, 0, 2): "622",
        (0, 0, 0, 6, 0, 1, 1, 2, 0, 2): "6mm",
        (2, 0, 0, 4, 0, 1, 3, 2, 0, 0): "-62m",
        (2, 0, 2, 7, 1, 1, 7, 2, 0, 2): "6/mmm",
        (0, 0, 0, 0, 0, 1, 3, 8, 0, 0): "23",
        (0, 0, 8, 3, 1, 1, 3, 8, 0, 0): "m-3",
        (0, 0, 0, 0, 0, 1, 9, 8, 6, 0): "432",
        (0, 6, 0, 6, 0, 1, 3, 8, 0, 0): "-43m",
        (0, 6, 8, 9, 1, 1, 9, 8, 6, 0): "m-3m",
    }[key]


def find_point_group_order(s: str) -> int:
    """Find order of point group."""
    return {
        "1": 1,
        "-1": 2,
        "2": 2,
        "m": 2,
        "2/m": 4,
        "222": 4,
        "mm2": 4,
        "mmm": 8,
        "4": 4,
        "-4": 4,
        "4/m": 8,
        "422": 8,
        "4mm": 8,
        "-42m": 8,
        "4/mmm": 16,
        "3": 3,
        "-3": 6,
        "32": 6,
        "3m": 6,
        "-3m": 12,
        "6": 6,
        "-6": 6,
        "-6/m": 12,
        "622": 12,
        "6mm": 12,
        "-62m": 12,
        "6/mmm": 24,
        "23": 12,
        "m-3": 24,
        "432": 24,
        "-43m": 24,
        "m-3m": 48,
    }[s]
