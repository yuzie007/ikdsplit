import itertools

import numpy as np
import pandas as pd
import yaml
from ase.spacegroup import crystal, get_spacegroup


def make_atoms(df: pd.DataFrame, spacegroup: int, cell: np.ndarray):
    # hexagonal cell (consistent with Bilbao) for rhombohedral space groups
    setting = 1 if (3 <= spacegroup <= 15 or 143 <= spacegroup <= 194) else 2

    symbols = df["symbol"].str.strip().values
    basis = df[["x", "y", "z"]]
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
    return atoms


def make_images(df: pd.DataFrame, spacegroup: int, cell: np.ndarray):
    n = len(df)

    for included in itertools.product([True, False], repeat=n):
        df_included = df[list(included)]
        atoms = make_atoms(df_included, spacegroup, cell)
        spacegroup_actual = get_spacegroup(atoms).todict()["number"]
        print(included, spacegroup_actual)
        if spacegroup_actual == spacegroup:
            break
    else:
        raise RuntimeError
    atoms.write("POSCAR", direct=True)


def add_arguments(parser):
    parser.add_argument(
        "-s",
        "--symbols",
        nargs="+",
    )


def run(args):
    cell = np.loadtxt("cell.dat")

    with open("wycksplit.yaml") as f:
        mapping = yaml.safe_load(f)

    spacegroup = mapping["space_group_number_sub"]

    df = pd.read_csv("atoms_conventional.csv", skipinitialspace=True)

    if args.symbols:
        df_included = df[df["symbol"].isin(args.symbols)]
        atoms = make_atoms(df_included, spacegroup, cell)
        atoms.write("POSCAR_test", direct=True)
    else:
        make_images(df, spacegroup, cell)
