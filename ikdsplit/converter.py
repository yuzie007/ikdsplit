"""Convert `atoms_conventional.csv` for the subgroup."""

import argparse
import pathlib
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


def convert() -> None:
    """Convert `atoms_conventional.csv` for the subgroup."""
    cell = np.loadtxt("../cell.dat")

    with pathlib.Path("wycksplit.toml").open("rb") as f:
        wycksplit = tomllib.load(f)
    wycksplit = parse_wycksplit(wycksplit)

    df = pd.read_csv("../atoms_conventional.csv", skipinitialspace=True)

    cell = (cell.T @ wycksplit["basis_change"]).T

    np.savetxt("cell.dat", cell, fmt="%24.18f")

    ds = []
    for _, s in df.iterrows():
        for wyckoff, op in wycksplit["wycksplit"][s["wyckoff"]].items():
            d = {}
            d["symbol"] = s["symbol"]
            d["fill"] = s["fill"]
            d["wyckoff"] = wyckoff.rstrip(string.digits)
            d["wyckoff_orig"] = s["wyckoff"]

            # representative for split wyckoff sites
            xyz = s[["x", "y", "z"]].to_numpy(float)
            xyz = op["basis_change"] @ xyz + op["origin_shift"]
            d["x"], d["y"], d["z"] = xyz

            ds.append(d)

    df = pd.DataFrame(ds)

    # change coordinates for subgroup
    xyzs = df[["x", "y", "z"]].to_numpy() - wycksplit["origin_shift"]
    xyzs = (np.linalg.inv(wycksplit["basis_change"]) @ xyzs.T).T
    xyzs -= np.rint(xyzs)
    df[["x", "y", "z"]] = xyzs

    df = format_df(df)
    filename = "atoms_conventional.csv"
    df.to_csv(filename, float_format="%24.18f", index=False)


def add_arguments(parser: argparse.ArgumentParser) -> None:
    pass


def run(args: argparse.Namespace) -> None:
    convert()
