"""Run recursively."""

import argparse
import pathlib

from ikdsplit.filler import fill
from ikdsplit.io import make_default_config, parse_config
from ikdsplit.regressor import regress
from ikdsplit.sorter import sort_all
from ikdsplit.spacegroup import (
    find_crystal_class,
    find_point_group_order,
    get_subgroups,
)
from ikdsplit.utils import cd, count_configurations, print_group


def recur_run(
    config: dict,
    supergroup: int,
    group: int,
    level: int,
    criteria: dict[str, int],
) -> None:
    """Run each subgroup recursively."""
    print_group(group, level, end=" ")

    dn = pathlib.Path(f"{group:03d}")
    dn.mkdir(parents=True, exist_ok=True)
    with cd(dn):
        order = find_point_group_order(find_crystal_class(group))
        print(f"(order: {order})", end=" ")
        ncs = count_configurations()
        print(f"({ncs} configurations)", end=" ")
        if order < criteria["min_order"]:
            print("... skipped")
            return
        if ncs > criteria["max_configurations"]:
            print("... skipped")
            return
        print()

        fill()

        regress()

        d = config["sort"]
        sort_all(d["reference"])

        if 0 < criteria["max_level"] <= level:
            return

        subgroups = get_subgroups(group)
        for subgroup in subgroups:
            recur_run(config, group, subgroup, level + 1, criteria)


def get_default_criteria() -> dict[str, int]:
    """Get default criteria."""
    return {
        "max_level": 1,
        "min_order": 4,
        "max_configurations": 2**12,  # 4096
    }


def start(criteria: dict[str, int]) -> None:
    """Start calculations."""
    config = make_default_config()
    config.update(parse_config())

    criteria = get_default_criteria() | criteria

    print("run ...")
    recur_run(config, None, config["space_group_number"], 0, criteria)
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
    criteria = {
        "max_level": args.level,
        "min_order": args.order,
        "max_configurations": args.configurations,
    }
    start(criteria)


def main() -> None:
    """Run as a script."""
    formatter_class = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=formatter_class)
    add_arguments(parser)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
