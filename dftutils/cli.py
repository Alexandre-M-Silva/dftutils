import os
import sys
import warnings
import click
from monty.serialization import dumpfn, loadfn

from dftutils.polarization import PolarizationPlotter
from dftutils.utils import match_structure_indices

def CommandWithConfigFile(config_file_param_name,):  # can also set CLI options using config file
    """
    Set CLI options using config file.

    Args:
        config_file_param_name (:obj:`str`):
            name of config file
    """

    class CustomCommandClass(click.Command):
        def invoke(self, ctx):
            config_file = ctx.params[config_file_param_name]
            if config_file is not None:
                config_data = loadfn(config_file)
                for param, _val in ctx.params.items():
                    if (
                        ctx.get_parameter_source(param) == click.core.ParameterSource.DEFAULT
                        and param in config_data
                    ):
                        ctx.params[param] = config_data[param]
            return super().invoke(ctx)

    return CustomCommandClass

# CLI Commands:
CONTEXT_SETTINGS = {"help_option_names": ["-h", "--help"]}

@click.group("dftutils", context_settings=CONTEXT_SETTINGS, no_args_is_help=True)
def dftutils():
    """DftUtils: Utilities for VASP DFT calculation workflows."""

@dftutils.command(
    name="polarization",
    context_settings=CONTEXT_SETTINGS,
    no_args_is_help=True,
    cls=CommandWithConfigFile("config"),
)
@click.option(
    "--path",
    "-p",
    help="Path to polarization calculation.",
    required=True,
    type=click.Path(exists=True, dir_okay=True),
)
@click.option(
    "--save",
    "-s",
    help="Highlight midpoint switch pathway.",
    required=False,
    default=False,
    is_flag=True,
    show_default=True,
)
@click.option(
    "--branch",
    "-b",
    help="Highlight midpoint switch branch.",
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
def polarization(path, save, branch, config):
    user_settings = loadfn(config) if config is not None else {}
    func_args = list(locals().keys())

    if user_settings:
        valid_args = [
            "path",
            "save",
            "branch",
            "config",
        ]
        for key in func_args:
            if key in user_settings:
                user_settings.pop(key, None)

        for key in list(user_settings.keys()):
            # remove non-sense keys from user_settings
            if key not in valid_args:
                user_settings.pop(key)

    pol = PolarizationPlotter(path)
    pol.plot(branch, save)

@dftutils.command(
    name="match",
    context_settings=CONTEXT_SETTINGS,
    no_args_is_help=True,
    cls=CommandWithConfigFile("config"),
)
@click.option(
    "--structure",
    "-s",
    help="Path to structure you want to have matching indices.",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "--match",
    "-m",
    help="Path to structure you want to match indices to.",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "--out",
    "-o",
    help="Path of output file.",
    required=False,
    default=None,
    type=click.Path(exists=False, dir_okay=False),
)
@click.option(
    "--config",
    "-conf",
    help="Config file for advanced settings.",
    default=None,
    type=click.Path(exists=True, dir_okay=False),
    show_default=True,
)
def match(structure, match, out, config):
    user_settings = loadfn(config) if config is not None else {}
    func_args = list(locals().keys())

    if user_settings:
        valid_args = [
            "structure",
            "match",
            "out",
            "config,"
        ]
        for key in func_args:
            if key in user_settings:
                user_settings.pop(key, None)

        for key in list(user_settings.keys()):
            # remove non-sense keys from user_settings
            if key not in valid_args:
                user_settings.pop(key)

    output_path = os.path.join(structure, "_matched")
    if not out is None:
        output_path = out
            
    match_structure_indices(structure, match, out)
