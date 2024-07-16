import itertools
import string

import numpy as np
import pandas as pd
import yaml
from ase.spacegroup import crystal, get_spacegroup


def format_df(df):
    for k in df.columns:
        if df[k].dtype == object:
            df[k] = df[k].str.pad(12)
    return df


def make_atoms(df: pd.DataFrame, spacegroup: int, cell: np.ndarray):
    n = len(df)

    # hexagonal cell (consistent with Bilbao) for rhombohedral space groups
    setting = 1 if (3 <= spacegroup <= 15 or 143 <= spacegroup <= 194) else 2

    for included in itertools.product([True, False], repeat=n):
        df_included = df[list(included)]
        symbols = df_included["symbol"].str.strip().values
        basis = df_included[["x", "y", "z"]]
        if len(symbols) == 0:
            return
        try:
            atoms = crystal(
                symbols,
                basis=basis,
                spacegroup=spacegroup,
                cell=cell,
                setting=setting,
            )
        except Exception:
            atoms = crystal(
                symbols,
                basis=basis,
                spacegroup=spacegroup,
                cell=cell,
            )
        spacegroup_actual = get_spacegroup(atoms).todict()["number"]
        print(included, spacegroup_actual)
        if spacegroup_actual == spacegroup:
            break
    else:
        raise RuntimeError
    atoms.write("POSCAR", direct=True)


def add_arguments(parser):
    pass


def run(args):
    cell = np.loadtxt("../cell.dat")

    with open("wycksplit.yaml") as f:
        mapping = yaml.safe_load(f)

    spacegroup_sup = mapping["space_group_number_sup"]
    spacegroup_sub = mapping["space_group_number_sub"]

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

    make_atoms(df, spacegroup_sub, cell)
