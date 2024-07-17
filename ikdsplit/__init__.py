"""Module ikdtools"""

import argparse

import ikdsplit.atoms
import ikdsplit.converter
import ikdsplit.regressor
import ikdsplit.sorter


def main():
    """main"""
    formatter_class = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=formatter_class)
    subparsers = parser.add_subparsers(dest="command")

    commands = {
        "symmetry": ikdsplit.atoms,
        "convert": ikdsplit.converter,
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
