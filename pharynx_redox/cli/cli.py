"""
command-line interface for running an analysis
"""

from pathlib import Path
from pharynx_redox import experiment


def analyze_experiment(experiment_dir: str):
    experiment_dir = Path(experiment_dir)
    exp = experiment.Experiment(experiment_dir).full_pipeline()
