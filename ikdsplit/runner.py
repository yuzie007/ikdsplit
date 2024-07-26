"""Run recursively."""

import argparse
import pathlib

from ikdsplit.filler import fill
from ikdsplit.io import make_default_config, parse_config
from ikdsplit.preparer import recur_prepare
from ikdsplit.regressor import regress
from ikdsplit.sorter import sort_all
from ikdsplit.utils import cd, get_subgroups, print_group


def recur_run(
    config: dict,
    supergroup: int,
    group: int,
    level: int,
    max_level: int,
) -> None:
    """Run each subgroup recursively."""
    print_group(group, level)

    dn = pathlib.Path(f"{group:03d}")
    dn.mkdir(parents=True, exist_ok=True)
    with cd(dn):
        fill()

        regress()

        d = config["sort"]
        sort_all(d["reference"])

        if 0 <= max_level <= level:
            return

        subgroups = get_subgroups(group)
        for subgroup in subgroups:
            recur_run(config, group, subgroup, level + 1, max_level)


def start(max_level: int = 1) -> None:
    """Start calculations."""
    config = make_default_config()
    config.update(parse_config())

    print("prepare ...")
    recur_prepare(config, None, config["space_group_number"], 0, max_level)
    print()

    print("run ...")
    recur_run(config, None, config["space_group_number"], 0, max_level)
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
