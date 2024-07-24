"""Run recursively."""

import argparse
import pathlib
import shutil
import tomllib

from ikdsplit.converter import convert
from ikdsplit.filler import fill
from ikdsplit.regressor import regress
from ikdsplit.sorter import sort_all
from ikdsplit.utils import cd, get_subgroups


def make_default_config() -> dict:
    """Make default `config`."""
    config = {}
    config["regress"] = {"transformations": []}
    config["sort"] = {"reference": None}
    return config


def parse_config(config: dict) -> dict:
    """Parse `ikdsplit.toml`."""
    if config["sort"]["reference"] is not None:
        path = pathlib.Path(config["sort"]["reference"])
        config["sort"]["reference"] = path.resolve()
    return config


def write_config(config: dict) -> None:
    """Write `ikdsplit.toml` for subgroup."""
    space_group_number = config["space_group_number"]
    with pathlib.Path("ikdsplit.toml").open("w", encoding="utf-8") as f:
        f.write(f"space_group_number = {space_group_number:d}\n")
        f.write("\n")

        f.write("[fill]\n")
        for k, v in config["fill"].items():
            r = repr(v)
            f.write(f"{k} = {r}\n")


def print_group(group: int, level: int) -> None:
    """Print group with indent."""
    s = ""
    if level > 1:
        s += "   " * (level - 1)
    if level > 0:
        s += "-> "
    s += f"{group:03d}"
    print(s)


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


def recur(
    config: dict,
    supergroup: int,
    group: int,
    level: int,
    max_level: int,
) -> None:
    """Run each subgroup recursively."""
    print_group(group, level)

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
            shutil.copy2("../cell.dat", ".")

        write_config(config | {"space_group_number": group})

        fill(group, config["fill"])

        d = config["regress"]
        regress(d["transformations"])

        d = config["sort"]
        sort_all(d["reference"])

        if 0 <= max_level <= level:
            return

        subgroups = get_subgroups(group)
        for subgroup in subgroups:
            recur(config, group, subgroup, level + 1, max_level)


def start(max_level: int = 1) -> None:
    """Start calculations."""
    config = make_default_config()
    with pathlib.Path("ikdsplit.toml").open("rb") as f:
        config.update(tomllib.load(f))
    config = parse_config(config)
    recur(config, None, config["space_group_number"], 0, max_level)


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
