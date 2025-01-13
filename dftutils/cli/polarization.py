import os
import click
from monty.serialization import loadfn

from dftutils.cli import *
from dftutils.core.polarization import *
from dftutils.core.utils import use_matplotlib_style

@dftutils.command(
    name="polarization",
    context_settings=CONTEXT_SETTINGS,
    no_args_is_help=False,
    cls=CommandWithConfigFile("config"),
)
@click.option(
    "--scatter",
    "-s",
    help="Builds a scatter for the --path that points to a directory cointaining numerically numbered subdirectories with polarization calculations.",
    required=False,
    default=False,
    is_flag=True,
    show_default=True,
)
@click.option(
    "--path",
    "-p",
    help="Path that points to a directory cointaining numerically numbered subdirectories with polarization calculations.",
    required=False,
    type=click.Path(exists=True, dir_okay=True),
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
    "-c",
    help="Path to POSCAR/CONTCAR.",
    required=False,
    type=click.Path(exists=True, dir_okay=False),
)
@click.option(
    "--axis",
    "-a",
    help="Scatter axis.",
    required=False,
    type=int,
    default=None,
)
@click.option(
    "--raw",
    "-r",
    help="Create figures of Px, Py and Pz scatter without determining the switching pathway.",
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
def polarization(scatter, path, outcar, poscar, axis, raw, config):
    user_settings = loadfn(config) if config is not None else {}
    func_args = list(locals().keys())

    if user_settings:
        valid_args = [
            "scatter",
            "path",
            "outcar",
            "poscar",
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

    if scatter:
        if path is None:
            path = '.'

        use_matplotlib_style()

        filenames = ['Px', 'Py', 'Pz'] 

        pol = Polarization(path)
        if axis is None:
            for i, fn in enumerate(filenames):
                pol.plot(os.path.join(path, f'{fn}.png'), axis=i, raw=raw)
        else:
            pol.plot(os.path.join(path, f'{filenames[axis]}.png'), axis=axis, raw=raw)
    else:
        if path is not None:
            outcar = os.path.join(path, "OUTCAR")
            poscar = [os.path.join(path, "CONTCAR"),
                        os.path.join(path, "POSCAR")]
        else:
            if outcar is None:
                outcar = "OUTCAR"
            if poscar is None:
                poscar = ["CONTCAR", "POSCAR"]

        p, q = polarization_from_outcar_structure(outcar, poscar)
        if axis is None:
            for j in range(5, -6, -1):
                print(f"{p[0]+j*q[0]:10.2f}\t{p[1]+j*q[1]:10.2f}\t{p[2]+j*q[2]:10.2f}")
        else:
            for j in range(5, -6, -1):
                print(f"{p[axis]+j*q[axis]:10.2f}")