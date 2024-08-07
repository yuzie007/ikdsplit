"""Run recursively."""

import pathlib
import typing

from ikdsplit.io import make_default_config, parse_config
from ikdsplit.spacegroup import get_subgroups
from ikdsplit.splitter import check_criteria, get_default_criteria
from ikdsplit.utils import cd, print_group


def run_recursive(
    fun: typing.Callable,
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
        if not check_criteria(group, criteria):
            return

        fun()

        if 0 < criteria["max_level"] <= level:
            return

        subgroups = get_subgroups(group)
        for subgroup in subgroups:
            run_recursive(fun, config, group, subgroup, level + 1, criteria)


def start_recursive(fun: typing.Callable, criteria: dict[str, int]) -> None:
    """Start calculations."""
    config = make_default_config()
    config.update(parse_config())

    criteria = get_default_criteria() | criteria

    print("run ...")
    run_recursive(fun, config, None, config["space_group_number"], 0, criteria)
    print()
