"""Run recursively."""

import argparse
import copy
import math
import pathlib
import shutil
import tomllib

import pandas as pd

from ikdsplit.converter import convert
from ikdsplit.io import make_default_config, parse_config, write_config
from ikdsplit.spacegroup import invert
from ikdsplit.utils import cd, get_subgroups, print_group


def write_wycksplit_toml_orig(group: int) -> None:
    """Write `wycksplit.toml` for the parent group."""
    s = f"""\
space_group_number_sup = {group}
space_group_number_sub = {group}
basis_change = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
origin_shift = [0.0, 0.0, 0.0]
"""
    with pathlib.Path("wycksplit.toml").open("w", encoding="utf-8") as f:
        f.write(s)


def update_config(superconfig: dict) -> dict:
    """Update `config` for subgroup."""
    config = copy.deepcopy(superconfig)

    with pathlib.Path("wycksplit.toml").open("rb") as f:
        wycksplit = tomllib.load(f)

    config["space_group_number"] = wycksplit["space_group_number_sub"]

    config["cell"] = (config["cell"].T @ wycksplit["basis_change"]).T

    op = invert(wycksplit["basis_change"], wycksplit["origin_shift"])
    op = {"basis_change": op[0].tolist(), "origin_shift": op[1].tolist()}
    config["regress"]["transformations"].insert(0, op)

    return config


def count_configurations() -> int:
    """Count number of atomic configurations."""
    with pathlib.Path("ikdsplit.toml").open("rb") as f:
        config = tomllib.load(f)
    df = pd.read_csv("atoms_conventional.csv", skipinitialspace=True)
    return math.prod([len(config["fill"][_]) for _ in df["symbol"]])


def recur_prepare(
    superconfig: dict,
    supergroup: int,
    group: int,
    level: int,
    max_level: int,
) -> None:
    """Prepare each subgroup recursively."""
    print_group(group, level, end=" ")

    src = pathlib.Path(__file__).parent / "database"
    dn = pathlib.Path(f"{group:03d}")
    dn.mkdir(parents=True, exist_ok=True)
    with cd(dn):
        if supergroup:
            fn = src / "wycksplit" / f"{supergroup:03d}_{group:03d}.toml"
            shutil.copy(fn, "wycksplit.toml")
            convert()
        else:
            write_wycksplit_toml_orig(group)
            shutil.copy2("../atoms_conventional.csv", ".")

        config = update_config(superconfig)

        write_config(config)

        print(f"({count_configurations()} configurations)")

        subgroups = get_subgroups(group)
        for subgroup in subgroups:
            recur_prepare(config, group, subgroup, level + 1, max_level)


def start(max_level: int = 1) -> None:
    """Start calculations."""
    config = make_default_config()
    config.update(parse_config())

    print("prepare ...")
    recur_prepare(config, None, config["space_group_number"], 0, max_level)
    print()


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Add arguments."""
    parser.add_argument(
        "-l",
        "--level",
        default=1,
        type=int,
        help="level up to which maximal subgroups are checked",
    )


def run(args: argparse.Namespace) -> None:
    """Run."""
    start(args.level)
