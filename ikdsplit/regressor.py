import functools
import os

import ase.io
import numpy as np
import pandas as pd
import yaml
from ase import Atoms
from ase.build import make_supercell

from ikdsplit.utils import format_df


def change_coordinates(
    atoms: Atoms,
    basis_change: np.ndarray,
    origin_shift: np.ndarray,
) -> tuple[Atoms, pd.DataFrame]:
    """Change the coordinate system."""
    atoms.set_scaled_positions(atoms.get_scaled_positions() - origin_shift)
    atoms = make_supercell(atoms, basis_change.T, order="atom-major")
    return atoms


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
    basis_change = basis_change_1 @ basis_change_0
    origin_shift = basis_change_1 @ origin_shift_0 + origin_shift_1
    return basis_change, origin_shift


def cumulate_coordinate_change(
    basis_change_last: np.ndarray,
    origin_shift_last: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Cumulate changes of the coordinate system."""
    ops = []
    fn = ""
    while True:
        fn = "wycksplit.yaml" if not fn else os.path.join("..", fn)
        if not os.path.isfile(fn):
            break
        with open(fn, encoding="utf-8") as f:
            d = yaml.safe_load(f)
        ops.append(invert(d["rotation"], d["translation"]))
    ops.append((basis_change_last, origin_shift_last))
    return functools.reduce(multiply, ops)


def add_arguments(parser):
    parser.add_argument(
        "images",
        nargs="+",
        help="POSCAR files to be updated",
    )
    parser.add_argument(
        "--rotation",
        nargs=9,
        default=[1, 0, 0, 0, 1, 0, 0, 0, 1],
        type=int,
    )
    parser.add_argument(
        "--translation",
        nargs=3,
        default=[0.0, 0.0, 0.0],
        type=float,
    )


def run(args):
    basis_change_last = np.array(args.rotation).reshape(3, 3)
    origin_shift_last = np.array(args.translation)
    basis_change, origin_shift = cumulate_coordinate_change(
        basis_change_last,
        origin_shift_last,
    )
    for fin in args.images:
        print(fin)
        atoms = ase.io.read(fin)
        atoms = change_coordinates(atoms, basis_change, origin_shift)
        atoms.write(f"R{fin}", direct=True)

    df = pd.read_csv("atoms_conventional.csv", skipinitialspace=True)
    df[["x", "y", "z"]] -= origin_shift
    df[["x", "y", "z"]] @= basis_change.T
    df = format_df(df)
    df.to_csv("atoms_regressed.csv", float_format="%24.18f", index=False)
