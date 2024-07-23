"""Convert `atoms_conventional.csv`."""

import string
import tomllib

import numpy as np
import pandas as pd

from ikdsplit.utils import format_df


def parse_wycksplit(d: dict) -> dict:
    """Parse `wycksplit.toml`."""
    reps = d["coset_representatives"]
    for k0, v0 in d["wycksplit"].items():
        d["wycksplit"][k0] = {k1: reps[v1] for k1, v1 in v0.items()}
    return d


def convert():
    cell = np.loadtxt("../cell.dat")

    with open("wycksplit.toml", "rb") as f:
        mapping = tomllib.load(f)
    mapping = parse_wycksplit(mapping)

    df = pd.read_csv("../atoms_conventional.csv", skipinitialspace=True)

    cell = (cell.T @ mapping["basis_change"]).T

    np.savetxt("cell.dat", cell, fmt="%24.18f")

    symbols = []
    wyckoffs = []
    wyckoffs_orig = []
    basis = []
    for i, s in df.iterrows():
        for wyckoff, d in mapping["wycksplit"][s["wyckoff"]].items():
            symbols.append(s["symbol"])
            wyckoffs.append(wyckoff.rstrip(string.digits))
            wyckoffs_orig.append(s["wyckoff"])

            # representative for split wyckoff sites
            xyz = s[["x", "y", "z"]].to_numpy(float)
            basis.append(d["basis_change"] @ xyz + d["origin_shift"])

    basis = np.array(basis) - mapping["origin_shift"]
    basis = (np.linalg.inv(mapping["basis_change"]) @ basis.T).T
    basis -= np.rint(basis)

    d = {
        "symbol": symbols,
        "wyckoff": wyckoffs,
        "wyckoff_orig": wyckoffs_orig,
        "x": basis[:, 0].tolist(),
        "y": basis[:, 1].tolist(),
        "z": basis[:, 2].tolist(),
    }
    df = pd.DataFrame(d)
    df = format_df(df)
    filename = "atoms_conventional.csv"
    df.to_csv(filename, float_format="%24.18f", index=False)


def add_arguments(parser):
    pass


def run(args):
    convert()
