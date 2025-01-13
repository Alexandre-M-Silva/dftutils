import os
import click
from monty.serialization import loadfn

from dftutils.cli import *
from dftutils.core.utils import *

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
    "--config",
    "-conf",
    help="Config file for advanced settings.",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    show_default=True,
)
def match(a, b, config):
    user_settings = loadfn(config) if config is not None else {}
    func_args = list(locals().keys())

    if user_settings:
        valid_args = [
            "a",
            "b",
            "config,"
        ]
        for key in func_args:
            if key in user_settings:
                user_settings.pop(key, None)

        for key in list(user_settings.keys()):
            # remove non-sense keys from user_settings
            if key not in valid_args:
                user_settings.pop(key)

    if isinstance(a, Structure) and isinstance(b, Structure):
        match_indices_from_structs(a, b)
    else:
        match_indices_from_paths(a, b)
