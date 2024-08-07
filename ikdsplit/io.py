"""IO."""

import io
import pathlib
import tomllib

import numpy as np
from ase.geometry import cellpar_to_cell


def make_default_config() -> dict:
    """Make default `config`."""
    config = {}
    config["regress"] = {"transformations": []}
    config["sort"] = {"reference": None}
    return config


def _parse_cell(cell: np.ndarray) -> np.ndarray:
    """Parse cell."""
    if isinstance(cell, int | float):
        return cell * np.eye(3)
    cell = np.array(cell)
    if cell.shape != (3, 3):
        return cellpar_to_cell(cell)
    return cell


def parse_config() -> dict:
    """Parse `ikdsplit.toml`."""
    with pathlib.Path("ikdsplit.toml").open("rb") as f:
        config = tomllib.load(f)

    config["cell"] = _parse_cell(config["cell"])

    if (
        "sort" in config
        and "reference" in config["sort"]
        and config["sort"]["reference"] is not None
    ):
        path = pathlib.Path(config["sort"]["reference"])
        config["sort"]["reference"] = path.resolve().absolute().as_posix()

    return config


def _write_cell(cell: np.ndarray, f: io.FileIO) -> None:
    """Write cell."""
    _ = "24.18f"
    f.write("cell = [\n")
    f.write(f"  [{cell[0, 0]:{_}},{cell[0, 1]:{_}},{cell[0, 2]:{_}}],\n")
    f.write(f"  [{cell[1, 0]:{_}},{cell[1, 1]:{_}},{cell[1, 2]:{_}}],\n")
    f.write(f"  [{cell[2, 0]:{_}},{cell[2, 1]:{_}},{cell[2, 2]:{_}}],\n")
    f.write("]\n")
    f.write("\n")


def write_config(config: dict) -> None:
    """Write `ikdsplit.toml` for subgroup."""
    space_group_number = config["space_group_number"]
    cell = config["cell"]
    with pathlib.Path("ikdsplit.toml").open("w", encoding="utf-8") as f:
        f.write(f"space_group_number = {space_group_number:d}\n")
        f.write("\n")

        _write_cell(cell, f)

        f.write("[fill]\n")
        for k, v in config["fill"].items():
            r = repr(v)
            f.write(f"{k} = {r}\n")
        f.write("\n")

        for transformation in config["regress"]["transformations"]:
            f.write("[[regress.transformations]]\n")
            for k, v in transformation.items():
                r = repr(v)
                f.write(f"{k} = {r}\n")
        f.write("\n")

        f.write("[sort]\n")
        for k, v in config["sort"].items():
            if v is None:
                continue
            r = repr(v)
            f.write(f"{k} = {r}\n")
