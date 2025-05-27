import os
import click
import ast
from pathlib import Path
from monty.serialization import loadfn

from dftutils.cli import *
from dftutils.core.utils import *
from dftutils.core.matching import *

class PythonLiteralOption(click.Option):
    def type_cast_value(self, ctx, value):
        if value is None:
            return None
        try:
            return ast.literal_eval(value)
        except:
            raise click.BadParameter(value)


@dftutils.command(
    name="match",
    context_settings=CONTEXT_SETTINGS,
    no_args_is_help=True,
    cls=CommandWithConfigFile("config"),
)
@click.option(
    "--a",
    "-a",
    help="Path to first structure.",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "--b",
    "-b",
    help="Path to second structure.",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "--type",
    "-t",
    help="Type of matching 'all' matches all sites, 'select' matches sites defined in parameter 'selection'.",
    required=False,
    default='all',
    show_default=True,
    type=click.STRING,
)
@click.option(
    "--selection",
    "-sl",
    help="Selection of sites defined for 'select' type matching, either a list of integer sites or a list of elements.",
    cls=PythonLiteralOption,
    default=None,
    required=False,
    show_default=True,
)
@click.option(
    "--sort-first",
    "-s",
    help="Sort structures before matching.",
    required=False,
    default=False,
    is_flag=True,
    show_default=True,
)
@click.option(
    "--config",
    "-conf",
    help="Config file for advanced settings.",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    show_default=True,
)
def match(a, b, type, selection, sort_first, config):
    user_settings = loadfn(config) if config is not None else {}
    func_args = list(locals().keys())

    if user_settings:
        valid_args = [
            "a",
            "b",
            "type",
            "selection",
            "sort_first",
            "config,"
        ]
        for key in func_args:
            if key in user_settings:
                user_settings.pop(key, None)

        for key in list(user_settings.keys()):
            # remove non-sense keys from user_settings
            if key not in valid_args:
                user_settings.pop(key)

    match_config = {
        "type": type,
        "selection": selection,
        "sort_first": sort_first,
    }
    
    matcher = StructureMatcher.from_paths_and_config(a, b, match_config)
    matched1, matched2 = matcher.match()
    base, ext = os.path.splitext(a)
    matched1.to(base + "_matched" + ext)
    base, ext = os.path.splitext(b)
    matched2.to(base + "_matched" + ext)