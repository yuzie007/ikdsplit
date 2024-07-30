"""Run recursively."""

import argparse
import pathlib

from ikdsplit.filler import fill
from ikdsplit.io import make_default_config, parse_config
from ikdsplit.regressor import regress
from ikdsplit.sorter import sort_all
from ikdsplit.spacegroup import find_crystal_class, find_point_group_order
from ikdsplit.utils import cd, count_configurations, get_subgroups, print_group


def recur_run(
    config: dict,
    supergroup: int,
    group: int,
    level: int,
    max_level: int,
    min_order: int,
    max_configurations: int,
) -> None:
    """Run each subgroup recursively."""
    print_group(group, level, end=" ")
    order = find_point_group_order(find_crystal_class(group))
    print(f"(order: {order})", end=" ")
    ncs = count_configurations()
    print(f"({ncs} configurations)", end=" ")
    if order < min_order or ncs > max_configurations:
        print("... skipped")
        return
    else:
        print()

    dn = pathlib.Path(f"{group:03d}")
    dn.mkdir(parents=True, exist_ok=True)
    with cd(dn):
        fill()

        regress()

        d = config["sort"]
        sort_all(d["reference"])

        if 0 < max_level <= level:
            return

        subgroups = get_subgroups(group)
        for subgroup in subgroups:
            recur_run(
                config,
                group,
                subgroup,
                level + 1,
                max_level,
                min_order,
                max_configurations,
            )


def start(
    max_level: int = 1,
    min_order: int = 4,
    max_configurations: int = 2**12,  # 4096
) -> None:
    """Start calculations."""
    config = make_default_config()
    config.update(parse_config())

    print("run ...")
    recur_run(
        config,
        None,
        config["space_group_number"],
        0,
        max_level,
        min_order,
        max_configurations,
    )
    print()


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Add arguments."""
    parser.add_argument(
        "-l",
        "--level",
        default=0,
        type=int,
        help=(
            "level up to which maximal subgroups are checked "
            "(default (0): all levels are checked)"
        ),
    )
    parser.add_argument(
        "-o",
        "--order",
        default=4,
        type=int,
        help="minimum order of space group to be checked",
    )
    parser.add_argument(
        "-c",
        "--configurations",
        default=2**12,  # 4096
        type=int,
        help="maximum configurations to be checked",
    )


def run(args: argparse.Namespace) -> None:
    """Run."""
    start(args.level, args.order, args.configurations)
