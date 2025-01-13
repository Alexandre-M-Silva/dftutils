import click
from monty.serialization import loadfn

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
