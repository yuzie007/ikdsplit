import os

import ase.io
import numpy as np
import pandas as pd
import yaml
from ase import Atoms
from ase.build import make_supercell


def format_df(df):
    for k in df.columns:
        if df[k].dtype == object:
            df[k] = df[k].str.pad(12)
    return df


def modify_atoms(
    atoms: Atoms,
    df: pd.DataFrame,
    rotation: np.ndarray,
    translation: np.ndarray,
):
    supercell_matrix = np.linalg.inv(rotation).T
    atoms = make_supercell(atoms, supercell_matrix, order="atom-major")
    atoms.set_scaled_positions(atoms.get_scaled_positions() + translation)

    df[["x", "y", "z"]] @= np.array(rotation).T
    df[["x", "y", "z"]] += translation

    return atoms, df


def add_arguments(parser):
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
    atoms = ase.io.read("POSCAR")
    df = pd.read_csv("atoms_conventional.csv", skipinitialspace=True)

    filename = "wycksplit.yaml"
    with open(filename) as f:
        mapping = yaml.safe_load(f)
    atoms, df = modify_atoms(
        atoms,
        df,
        mapping["rotation"],
        mapping["translation"],
    )
    while True:
        filename = os.path.join("..", filename)
        if not os.path.isfile(filename):
            break
        with open(filename) as f:
            mapping = yaml.safe_load(f)
        atoms, df = modify_atoms(
            atoms,
            df,
            mapping["rotation"],
            mapping["translation"],
        )

    rotation = np.array(args.rotation).reshape(3, 3)
    translation = np.array(args.translation)
    atoms, df = modify_atoms(atoms, df, rotation, translation)

    atoms.write("SPOSCAR_regressed", direct=True)

    df = format_df(df)
    df.to_csv("atoms_regressed.csv", float_format="%24.18f", index=False)
