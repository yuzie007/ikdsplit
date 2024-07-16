import string

import numpy as np
import pandas as pd
import yaml


def format_df(df):
    for k in df.columns:
        if df[k].dtype == object:
            df[k] = df[k].str.pad(12)
    return df


def add_arguments(parser):
    pass


def run(args):
    cell = np.loadtxt("../cell.dat")

    with open("wycksplit.yaml") as f:
        mapping = yaml.safe_load(f)

    df = pd.read_csv(f"../atoms_conventional.csv", skipinitialspace=True)

    cell = (mapping["rotation"] @ cell.T).T

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
            xyz = s[["x", "y", "z"]].values.astype(float)
            basis.append(d["rotation"] @ xyz + d["translation"])

    basis = np.array(basis) - mapping["translation"]
    basis = (np.linalg.inv(mapping["rotation"]) @ basis.T).T
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
