import numpy as np
import xarray as xr

from pharynx_analysis import image_processing as ip
from pharynx_analysis import pharynx_io as pio


class Experiment:

    def __init__(self, raw_image_path: str):
        self.raw_image_path = raw_image_path


def get_non_tl(data_array):
    return data_array.where(data_array.wavelength != 'TL', drop=True)


class PairExperiment(Experiment):

    def __init__(self, raw_image_path: str, imaging_scheme: str, strain_map: [str]):
        super().__init__(raw_image_path)
        self.imaging_scheme = imaging_scheme
        self.strain_map = strain_map
        self.raw_image_data = pio.load_images(self.raw_image_path, imaging_scheme, strain_map)
        self.fl_images = self.raw_image_data.where(self.raw_image_data.wavelength != 'TL', drop=True)

        self.trimmed_profile_length = 100
        self.n_midline_pts = 200
        self.seg_threshold = 2000
        self.seg_stack = ip.segment_pharynxes(get_non_tl(self.raw_image_data))

        self.midline_map = {
            '410_1': '410_1',
            '470_1': '470_1',
            '410_2': '410_2',
            '470_2': '470_2',
        }

        self.regions = {
            'pm7': [.07, .13],
            'pm6': [.15, .21],
            'pm5': [.25, .47],
            'pm4': [.53, .64],
            'pm3': [.68, .93],
        }

        self.wavelengths = list(self.midline_map.keys())

        self.rot_fl = []
        self.rot_seg = []

        # Do the rotations
        for wvl in self.midline_map.keys():
            rot_fl, rot_seg = ip.center_and_rotate_pharynxes(get_non_tl(self.raw_image_data).sel(wavelength=wvl),
                                                             self.seg_stack.sel(wavelength=wvl))
            # crop_height=30, crop_width=70)
            self.rot_fl.append(rot_fl)
            self.rot_seg.append(rot_seg)

        self.rot_fl = xr.concat(self.rot_fl, dim='wavelength')
        self.rot_seg = xr.concat(self.rot_seg, dim='wavelength')

        # Calculate midlines
        self.midlines = {}
        for wvl in self.midline_map.keys():
            self.midlines[wvl] = ip.calculate_midlines(self.rot_seg.sel(wavelength=wvl))

        # Measure under midlines
        # TODO broken when cropped
        self.raw_intensity_data = []
        for img_wvl, mid_wvl in self.midline_map.items():
            img_stk = self.rot_fl.sel(wavelength=img_wvl)

            i_data = []
            for i in range(self.raw_image_data.shape[0]):
                midline = self.midlines[mid_wvl][i]
                xs = xr.DataArray(np.linspace(40, 120, self.n_midline_pts), dims='z')
                ys = xr.DataArray(midline(xs), dims='z')

                i_data.append(img_stk[i].interp(x=xs, y=ys))

            self.raw_intensity_data.append(xr.concat(i_data, dim='strain'))
        self.raw_intensity_data = xr.concat(self.raw_intensity_data, dim='wavelength')

        # Trim
        self.trimmed_intensity_data = ip.trim_profiles(self.raw_intensity_data, self.seg_threshold,
                                                       self.trimmed_profile_length)

        # Calculate Redox Measurements
        self.redox = self.trimmed_intensity_data
        # TODO calculate redox

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
