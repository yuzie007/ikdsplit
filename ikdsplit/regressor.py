import os

import ase.io
import numpy as np
import pandas as pd
import yaml
from ase import Atoms
from ase.build import make_supercell


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


def regress(
    atoms: Atoms,
    basis_change_last: np.ndarray,
    origin_shift_last: np.ndarray,
) -> Atoms:
    """Regress `Atoms` to the original unit cell."""
    atoms = ase.io.read("POSCAR")

    fn = ""
    while True:
        fn = "wycksplit.yaml" if not fn else os.path.join("..", fn)
        if not os.path.isfile(fn):
            break
        with open(fn, encoding="utf-8") as f:
            d = yaml.safe_load(f)
        basis_change, origin_shift = invert(d["rotation"], d["translation"])
        atoms = change_coordinates(atoms, basis_change, origin_shift)

    return change_coordinates(atoms, basis_change_last, origin_shift_last)


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
    basis_change = np.array(args.rotation).reshape(3, 3)
    origin_shift = np.array(args.translation)
    for fin in args.images:
        print(fin)
        atoms = ase.io.read(fin)
        atoms = regress(atoms, basis_change, origin_shift)
        atoms.write(f"S{fin}", direct=True)
