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

        self.midline_map = {
            '410_1': '410_1',
            '470_1': '470_1',
            '410_2': '410_2',
            '470_2': '470_2',
        }

        self.n_midline_pts = 500
        self.seg_threshold = 2000
        self.seg_stack = ip.segment_pharynxes(get_non_tl(self.raw_image_data))

        self.rot_fl = []
        self.rot_seg = []

        # Do the rotations
        for wvl in self.midline_map.keys():
            rot_fl, rot_seg = ip.center_and_rotate_pharynxes(get_non_tl(self.raw_image_data).sel(wavelength=wvl),
                                                             self.seg_stack.sel(wavelength=wvl))
            self.rot_fl.append(rot_fl)
            self.rot_seg.append(rot_seg)

        self.rot_fl = xr.concat(self.rot_fl, dim='wavelength')
        self.rot_seg = xr.concat(self.rot_seg, dim='wavelength')

        # Calculate midlines
        self.midlines = {}
        for wvl in self.midline_map.keys():
            self.midlines[wvl] = ip.calculate_midlines(self.rot_seg.sel(wavelength=wvl))

        # Measure under midlines
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

    def flip_at(self, idx):
        np.fliplr(self.rot_fl[:, idx])
        np.fliplr(self.rot_seg[:, idx])
        np.fliplr(self.raw_intensity_data[:, idx])


class TimeSeriesExperiment(Experiment):
    def __init__(self, raw_image_path: str):
        super().__init__(raw_image_path)
