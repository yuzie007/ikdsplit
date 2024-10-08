"""ikdsplit."""

import argparse

import ikdsplit.converter
import ikdsplit.filler
import ikdsplit.regressor
import ikdsplit.sorter
import ikdsplit.splitter


def main() -> None:
    """Command."""
    formatter_class = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=formatter_class)
    subparsers = parser.add_subparsers(dest="command")

    commands = {
        "split": ikdsplit.splitter,
        "convert": ikdsplit.converter,
        "fill": ikdsplit.filler,
        "regress": ikdsplit.regressor,
        "sort": ikdsplit.sorter,
    }
    for key, value in commands.items():
        value.add_arguments(subparsers.add_parser(key))

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        return

    commands[args.command].run(args)
