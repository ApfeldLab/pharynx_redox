from skfda import FDataGrid
from skfda.preprocessing.registration import elastic_registration


def register_profiles(raw_profile_data, lam=0.05):
    reg_profile_data = raw_profile_data.copy()
    for pair in raw_profile_data.pair.data:
        data410 = raw_profile_data.sel(pair=pair, wavelength='410')
        data470 = raw_profile_data.sel(pair=pair, wavelength='470')

        fd410 = FDataGrid(data410)
        fd470 = FDataGrid(data470)

        reg410 = elastic_registration(fd410, fd470, lam=lam)

        # the data matrix has an extra dimension in the last position... so we just index into that dimension directly
        # to obtain the desired shape
        reg_profile_data.loc[dict(pair=pair, wavelength='410')] = reg410.data_matrix[..., 0]

    return reg_profile_data
