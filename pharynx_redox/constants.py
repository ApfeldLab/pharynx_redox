import matlab

trimmed_regions = {
    "pm3": [0.07, 0.28],
    "pm4": [0.33, 0.45],
    "pm5": [0.53, 0.70],
    "pm6": [0.80, 0.86],
    "pm7": [0.88, 0.96],
}

untrimmed_regions = {
    "pm3": [0.18, 0.33],
    "pm4": [0.38, 0.46],
    "pm5": [0.52, 0.65],
    "pm6": [0.70, 0.75],
    "pm7": [0.76, 0.82],
}

untrimmed_regions_with_medial = {
    "pm3": [0.18, 0.33],
    "pm4": [0.38, 0.46],
    "pm5": [0.52, 0.65],
    "pm6": [0.70, 0.75],
    "pm7": [0.76, 0.82],
    "medial_axis": [0.18, 0.82],
}

trimmed_regions_with_medial = {
    "pm3": [0.07, 0.28],
    "pm4": [0.33, 0.45],
    "pm5": [0.53, 0.70],
    "pm6": [0.80, 0.86],
    "pm7": [0.88, 0.96],
    "medial_axis": [0.07, 0.96],
}

opt_trimmed_regions = {}

opt_untrimmed_regions = {
    "pm3": [0.28166667, 0.29833333],
    "pm4": [0.405, 0.42166667],
    "pm5": [0.52666667, 0.54333333],
    "pm6": [0.75333333, 0.77],
    "pm7": [0.81, 0.82666667],
}

# matlab_engine = matlab.engine.start_matlab()
