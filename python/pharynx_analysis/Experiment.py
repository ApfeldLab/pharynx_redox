from pharynx_analysis import pharynx_io as pio


class Experiment:

    def __init__(self, raw_image_path: str):
        self.raw_image_path = raw_image_path


class PairExperiment(Experiment):

    def __init__(self, raw_image_path: str, imaging_scheme: str, strain_map: [str]):
        super().__init__(raw_image_path)
        self.imaging_scheme = imaging_scheme
        self.strain_map = strain_map
        self.raw_image_data = pio.load_images(self.raw_image_path, imaging_scheme, strain_map)


class TimeSeriesExperiment(Experiment):

    def __init__(self, raw_image_path: str):
        super().__init__(raw_image_path)
