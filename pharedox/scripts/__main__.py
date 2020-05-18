import click
from pathlib import Path
import logging


@click.group()
@click.option("--debug/--no-debug", default=False)
def cli(debug):
    """Useful scripts for analyzing ratiometric microscopy data"""
    if debug:
        loglevel = logging.DEBUG
    else:
        loglevel = logging.INFO

    logging.basicConfig(
        format="%(asctime)s %(levelname)s:%(message)s",
        level=loglevel,
        datefmt="%I:%M:%S",
    )

    click.echo("Debug mode is %s" % ("ON" if debug else "OFF"))


@cli.command()
@click.argument("directory", type=click.Path(exists=True))
@click.option("--gui/--command-line", default=False)
def analyze(directory, gui):
    """Analyze an experiment"""
    from pharedox import experiment

    exp = experiment.Experiment(Path(directory))
    if gui:
        from pharedox.gui.gui_napari import App

        app = App(experiment=exp)
        app.run()
    else:
        exp.full_pipeline()


@cli.command()
@click.option("--output", "-o", default=None, type=click.Path(exists=False))
@click.argument("filepath", type=click.Path(exists=True))
def split_nc(filepath):
    raise NotImplementedError
