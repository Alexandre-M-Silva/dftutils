import os
import click
from monty.serialization import loadfn

from dftutils.cli import *
from dftutils.core.neb import *
from dftutils.core.utils import *

@dftutils.command(
    name="neb",
    context_settings=CONTEXT_SETTINGS,
    no_args_is_help=True,
    cls=CommandWithConfigFile("config"),
)
@click.option(
    "--path",
    "-p",
    help="Path to directory containing the desired pathway.",
    required=False,
    default=None,
    type=click.Path(exists=True, dir_okay=True),
)
@click.option(
    "--initial",
    "-i",
    help="Path to initial structure containing the desired pathway.",
    required=False,
    default=None,
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "--final",
    "-f",
    help="Path to final structure containing the desired pathway.",
    required=False,
    default=None,
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "--number",
    "-n",
    help="Number of images to generate.",
    required=False,
    type=int,
)
@click.option(
    "--output-path",
    "-o",
    help="Path to output the neb directories.",
    required=False,
    default=None,
    type=click.Path(exists=True, dir_okay=True),
)
@click.option(
    "--config",
    "-conf",
    help="Config file for advanced settings.",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    show_default=True,
)
def neb(path, initial, final, number, output_path, config):
    user_settings = loadfn(config) if config is not None else {}
    func_args = list(locals().keys())

    if user_settings:
        valid_args = [
            "path",
            "strain",
            "config,"
        ]
        for key in func_args:
            if key in user_settings:
                user_settings.pop(key, None)

        for key in list(user_settings.keys()):
            # remove non-sense keys from user_settings
            if key not in valid_args:   
                user_settings.pop(key)

    if not path is None:
        neb = Neb.from_path(path, number)
        neb.to_path(output_path)
    elif path is None and not initial is None and not final is None:
        neb = Neb.from_initial_and_final(initial, final, number)
        neb.to_path(output_path)
    else:
        raise Exception(
            f"Path or initial and final images were not defined."
        )

