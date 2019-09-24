from pathlib import Path

import click
import experiment


@click.command()
@click.option("--count", default=1, help="Number of greetings.")
@click.option("--name", prompt="Your name", help="The person to greet.")
def hello(count, name):
    """Simple program that greets NAME for a total of COUNT times."""
    for _ in range(count):
        click.echo("Hello %s!" % name)


def pipeline(experiment_dir: str, imaging_scheme: str, register: bool):
    experiment_dir = Path(experiment_dir)
    experiment.PairExperiment(
        experiment_dir=experiment_dir,
        imaging_scheme=imaging_scheme,
        should_register=register,
    ).full_pipeline()
