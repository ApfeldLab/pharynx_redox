import numpy as np

from pharynx_analysis import image_processing as ip
from pharynx_analysis import pharynx_io as pio


class Experiment:
    """
    The abstract Experiment class

    """
    regions = {
        'pm7': [.07, .13],
        'pm6': [.15, .21],
        'pm5': [.25, .47],
        'pm4': [.53, .64],
        'pm3': [.68, .93],
    }

    # R -> OxD parameters
    r_min = 0.852
    r_max = 6.65
    instrument_factor = 0.171

    # OxD -> E parameters
    midpoint_potential = -265
    z = 2
    temperature = 22

    def __init__(self, raw_image_path: str):
        self.raw_image_path = raw_image_path


def get_non_tl(data_array):
    return data_array.where(data_array.wavelength != 'TL', drop=True)


class PairExperiment(Experiment):
    """
    This is the paired ratio experiment

    Attributes
    ----------
    raw_image_path
    imaging_scheme
    strain_map
    midline_smoothing
    """
    midline_map = {
        '410_1': '410_1',
        '470_1': '470_1',
        '410_2': '410_2',
        '470_2': '470_2',
    }

    wavelengths = list(midline_map.keys())

    image_display_order = [
        '410_1', '470_1', 'r1',
        '410_2', '470_2', 'r2'
    ]

    trimmed_profile_length = 100
    n_midline_pts = 200
    seg_threshold = 2000

    def __init__(self, raw_image_path: str, imaging_scheme: str, strain_map: [str], midline_smoothing=1e8):
        super().__init__(raw_image_path)
        self.imaging_scheme = imaging_scheme
        self.strain_map = strain_map
        self.raw_image_data = pio.load_images(self.raw_image_path, imaging_scheme, strain_map)
        self.seg_stack = ip.segment_pharynxes(self.raw_image_data)

        self.rot_fl = []
        self.rot_seg = []

        self.rot_fl, self.rot_seg = ip.center_and_rotate_pharynxes(self.raw_image_data, self.seg_stack)

        # TODO still need to align PA

        self.midlines = ip.calculate_midlines(self.rot_seg)

        # Measure under midlines
        # TODO broken when cropped

        step = self.raw_image_data.x.size // 6
        self.midline_xs = np.linspace(step, self.raw_image_data.x.size - step, self.n_midline_pts)
        self.raw_intensity_data = ip.measure_under_midlines(
            self.rot_fl, self.midlines, (step, self.raw_image_data.x.size - step), n_points=self.n_midline_pts
        )

        self.raw_intensity_data = ip.align_pa(self.raw_intensity_data)

        # Trim
        self.trimmed_intensity_data = ip.trim_profiles(
            self.raw_intensity_data, self.seg_threshold, self.trimmed_profile_length
        )

        # Calculate Redox Measurements
        self.r = self.trimmed_intensity_data.sel(wavelength='410') / self.trimmed_intensity_data.sel(wavelength='470')
        self.oxd = ip.r_to_oxd(self.r)
        self.e = ip.oxd_to_redox_potential(self.oxd)

    def flip_at(self, idx):
        np.fliplr(self.rot_fl[:, idx])
        np.fliplr(self.rot_seg[:, idx])
        np.fliplr(self.raw_intensity_data[:, idx])

    def get_scaled_region_boundaries(self):
        return {
            region: np.int_(self.trimmed_profile_length * np.asarray(self.regions[region]))
            for region in self.regions.keys()
        }


class TimeSeriesExperiment(Experiment):
    def __init__(self, raw_image_path: str):
        super().__init__(raw_image_path)


if __name__ == '__main__':
    img_path = "/Users/sean/code/wormAnalysis/data/paired_ratio_movement_data_sean/2017_02_22-HD233_SAY47/2017_02_22-HD233_SAY47.tif"
    strain_map_path = "/Users/sean/code/wormAnalysis/data/paired_ratio_movement_data_sean/2017_02_22-HD233_SAY47/indexer.csv"
    strains = pio.load_strain_map_from_disk(strain_map_path)
    ex = PairExperiment(img_path, "TL/470/410/470/410", strains)
