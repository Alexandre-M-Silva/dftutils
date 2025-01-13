import os
import click
from monty.serialization import loadfn

from dftutils.cli import *
from dftutils.core.utils import *
from dftutils.core.strain import *

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
        