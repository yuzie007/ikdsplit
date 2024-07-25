"""IO."""

import pathlib


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
        f.write("\n")

        for transformation in config["regress"]["transformations"]:
            f.write("[[regress.transformations]]\n")
            for k, v in transformation.items():
                r = repr(v)
                f.write(f"{k} = {r}\n")
