import os
import sys
import warnings
import click
from monty.serialization import dumpfn, loadfn
import numpy as np

from pymatgen.core.structure import Structure

from dftutils.polarization import *
from dftutils.utils import *
from dftutils.strain import *
from dftutils.neb import *

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
    no_args_is_help=False,
    cls=CommandWithConfigFile("config"),
)
@click.option(
    "--outcar",
    "-o",
    help="Path to OUTCAR.",
    required=False,
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "--poscar",
    "-p",
    help="Path to POSCAR/CONTCAR.",
    required=False,
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
def polarization(outcar, poscar, config):
    user_settings = loadfn(config) if config is not None else {}
    func_args = list(locals().keys())

    if user_settings:
        valid_args = [
            "outcar",
            "poscar",
            "config",
        ]
        for key in func_args:
            if key in user_settings:
                user_settings.pop(key, None)

        for key in list(user_settings.keys()):
            # remove non-sense keys from user_settings
            if key not in valid_args:
                user_settings.pop(key)

    if outcar is None:
        outcar = "OUTCAR"
    if poscar is None:
        poscar = ["CONTCAR", "POSCAR"]

    p, q = polarization_from_outcar_structure(outcar, poscar)
    for j in range(5, -6, -1):
        print(f"{p[0]+j*q[0]:10.2f}\t{p[1]+j*q[1]:10.2f}\t{p[2]+j*q[2]:10.2f}")

@dftutils.command(
    name="polarization-scatter",
    context_settings=CONTEXT_SETTINGS,
    no_args_is_help=True,
    cls=CommandWithConfigFile("config"),
)
@click.option(
    "--path",
    "-p",
    help="Path to polarization calculations root.",
    required=True,
    type=click.Path(exists=True, dir_okay=True),
)
@click.option(
    "--axis",
    "-a",
    help="Select scatter axis.",
    required=False,
    type=int,
    default=None,
)
@click.option(
    "--raw",
    "-r",
    help="Create figures of Px, Py and Pz scatter.",
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
def polarization_scatter(path, axis, raw, config):
    user_settings = loadfn(config) if config is not None else {}
    func_args = list(locals().keys())

    if user_settings:
        valid_args = [
            "path",
            "axis",
            "raw",
            "config",
        ]
        for key in func_args:
            if key in user_settings:
                user_settings.pop(key, None)

        for key in list(user_settings.keys()):
            # remove non-sense keys from user_settings
            if key not in valid_args:
                user_settings.pop(key)

    use_matplotlib_style()

    filenames = ['Px', 'Py', 'Pz'] 

    pol = Polarization(path)
    if axis is not None:
        for i, fn in enumerate(filenames):
            pol.plot(os.path.join(path, f'{fn}.png'), axis=i, raw=raw)
    else:
        pol.plot(os.path.join(path, f'{filenames[axis]}.png'), axis=axis, raw=raw)
        
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
        ts = np.linspace(0, 1, int(number))
        strains = [np.array(min_strain) + t*(np.array(max_strain) - np.array(min_strain)) for t in ts]
        structures = [apply_strain(structure, s) for s in strains]
        folders = [os.path.join(scan_output, format_strain_directory(s)) for s in strains]
        
        # Create directories and export structures
        for s, f in zip(structures, folders):
            if not os.path.exists(f):
                os.makedirs(f)
            s.to(os.path.join(f, "POSCAR"), fmt="poscar")
        
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

