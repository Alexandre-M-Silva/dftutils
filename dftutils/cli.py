import os
import sys
import warnings
import click
from monty.serialization import dumpfn, loadfn
import numpy as np

from pymatgen.core.structure import Structure

from dftutils.polarization import PolarizationPlotter
from dftutils.utils import match_structure_indices

from dftutils.strain import apply_strain
from dftutils.strain import format_strain_directory

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

    match_structure_indices(a, b)

@dftutils.command(
    name="strain",
    context_settings=CONTEXT_SETTINGS,
    no_args_is_help=True,
    cls=CommandWithConfigFile("config"),
)
@click.option(
    "--path",
    "-p",
    help="Path to structure.",
    required=True,
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "--strain",
    "-s",
    help="Strain components.",
    required=False,
    type=click.Tuple([float, float, float]),
)
@click.option(
    "--scan",
    "-sc",
    help="If the generation should be a scan.",
    default=False,
    is_flag=True,
    show_default=True,
)
@click.option(
    "--min-strain",
    "-mins",
    help="Min strain components.",
    type=click.Tuple([float, float, float]),
    default=None,
)
@click.option(
    "--max-strain",
    "-maxs",
    help="Max strain components.",
    type=click.Tuple([float, float, float]),
    default=None,
)
@click.option(
    "--number",
    "-n",
    help="Number of strain value to scan.",
    type=int,
    default=0,
)
@click.option(
    "--scan-output",
    "-so",
    help="Output path for scan folders.",
    default=None,
    type=click.Path(exists=True, dir_okay=True),
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
def strain(path, strain, scan, min_strain, max_strain, number, scan_output, config):
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

    structure = Structure.from_file(path)
    if not scan:
        new_structure = apply_strain(structure, strain)
        new_structure.to(path + "_strained", fmt="poscar")
    else:
        t = np.linspace(0, 1, number)
        strains = np.array(min_strain) + t*(np.array(max_strain) - np.array(min_strain))
        structures = [apply_strain(structure, s) for s in strains]
        folders = [os.path.join(scan_output, format_strain_directory(s)) for s in strain]
        
        # Create strain scan folders
        for sf in folders:
            if not os.path.exists(sf):
                os.makedirs(sf)
        
        # Create directories and export structures
        for s, f in zip(structures, folders):
            if not os.path.exists(f):
                os.makedirs(f)
            s.to(os.path.join(f, "POSCAR"), fmt="poscar")
        

