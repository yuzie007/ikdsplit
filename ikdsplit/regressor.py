import itertools
import os

import ase.io
import numpy as np
import pandas as pd
import yaml
from ase.build import make_supercell


def format_df(df):
    for k in df.columns:
        if df[k].dtype == object:
            df[k] = df[k].str.pad(12)
    return df


def modify_atoms(atoms, df, mapping):
    rotation = mapping["rotation"]
    translation = mapping["translation"]

    atoms = make_supercell(atoms, np.linalg.inv(rotation).T, order='atom-major')
    atoms.set_scaled_positions(atoms.get_scaled_positions() + translation)

    df[["x", "y", "z"]] @= np.array(rotation).T
    df[["x", "y", "z"]] += translation

    return atoms, df


def add_arguments(parser):
    pass


def run(args):
    atoms = ase.io.read("POSCAR")
    df = pd.read_csv("atoms_conventional.csv", skipinitialspace=True)

    filename = "wycksplit.yaml"
    with open(filename) as f:
        mapping = yaml.safe_load(f)
    atoms, df = modify_atoms(atoms, df, mapping)
    while True:
        filename = os.path.join("..", filename)
        if not os.path.isfile(filename):
            break
        with open(filename) as f:
            mapping = yaml.safe_load(f)
        atoms, df = modify_atoms(atoms, df, mapping)

    mapping = {}
    mapping['rotation'] = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    mapping['translation'] = [0.125, 0.125, 0.125]  # setting 1 -> 2
    atoms, df = modify_atoms(atoms, df, mapping)

    atoms.write("SPOSCAR_ref", direct=True)

    df = format_df(df)
    df.to_csv("atoms_ref.csv", float_format="%24.18f", index=False)
