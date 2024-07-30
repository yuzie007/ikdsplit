"""Run recursively."""

import argparse
import copy
import pathlib

import pandas as pd

from ikdsplit.converter import add_wyckoff, convert
from ikdsplit.io import (
    fetch_transformation,
    make_default_config,
    parse_config,
    write_config,
)
from ikdsplit.spacegroup import (
    find_crystal_class,
    find_point_group_order,
    invert,
)
from ikdsplit.utils import (
    cd,
    count_configurations,
    format_df,
    get_subgroups,
    print_group,
)


def update_config(superconfig: dict, spg_sup: int, spg_sub: int) -> dict:
    """Update `config` for subgroup."""
    config = copy.deepcopy(superconfig)

    config["space_group_number"] = spg_sub

    change_of_basis, origin_shift = fetch_transformation(spg_sup, spg_sub)

    config["cell"] = (config["cell"].T @ change_of_basis).T

    op = invert(change_of_basis, origin_shift)
    op = {"basis_change": op[0].tolist(), "origin_shift": op[1].tolist()}
    config["regress"]["transformations"].insert(0, op)

    return config


def recur_prepare(
    superconfig: dict,
    supergroup: int,
    group: int,
    level: int,
    max_level: int,
) -> None:
    """Prepare each subgroup recursively."""
    print_group(group, level, end=" ")

    dn = pathlib.Path(f"{group:03d}")
    dn.mkdir(parents=True, exist_ok=True)
    with cd(dn):
        if supergroup:
            convert(supergroup, group)
        else:
            fn = "../atoms_conventional.csv"
            df = pd.read_csv(fn, skipinitialspace=True)
            df = add_wyckoff(df, group)
            df = format_df(df)
            fn = "atoms_conventional.csv"
            df.to_csv(fn, float_format="%21.15f", index=False)

        config = update_config(superconfig, supergroup, group)

        write_config(config)

        order = find_point_group_order(find_crystal_class(group))
        print(f"(order: {order})", end=" ")
        ncs = count_configurations()
        print(f"({ncs} configurations)")

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
