import numpy as np
import xarray as xr
from scipy.sparse import csr_matrix

def save_impact_data_to_NetCDF(impact_object, filename, include_eai_exp=True, include_imp_mat=True):
    """
    Save eai_exp and/or imp_mat data from an impact object to a NetCDF file.

    Parameters:
        impact_object: The impact object containing the data.
        filename: The filename for the NetCDF file.
        include_eai_exp: Flag to include eai_exp data in the file.
        include_imp_mat: Flag to include imp_mat data in the file.
    """
    coords = {}
    data_vars = {}

    # Process eai_exp data if requested
    if include_eai_exp and hasattr(impact_object, 'eai_exp'):
        unique_lats, lat_inverse = np.unique(impact_object.coord_exp[:, 0], return_inverse=True)
        unique_lons, lon_inverse = np.unique(impact_object.coord_exp[:, 1], return_inverse=True)

        eai_exp_reshaped = np.full((len(unique_lats), len(unique_lons)), np.nan)
        eai_exp_reshaped[lat_inverse, lon_inverse] = impact_object.eai_exp

        coords.update({'latitude': unique_lats, 'longitude': unique_lons})
        data_vars.update({'eai_exp': (('latitude', 'longitude'), eai_exp_reshaped)})

    # Process imp_mat data if requested
    if include_imp_mat and hasattr(impact_object, 'imp_mat'):
        dense_imp_mat = impact_object.imp_mat.toarray() if isinstance(impact_object.imp_mat, csr_matrix) else impact_object.imp_mat
        num_time_steps, num_spatial_points = dense_imp_mat.shape

        imp_mat_3d = np.full((num_time_steps, len(unique_lats), len(unique_lons)), np.nan)
        for time_step in range(num_time_steps):
            for spatial_point in range(num_spatial_points):
                value = dense_imp_mat[time_step, spatial_point]
                lat = impact_object.coord_exp[spatial_point, 0]
                lon = impact_object.coord_exp[spatial_point, 1]
                lat_idx = np.where(unique_lats == lat)[0][0]
                lon_idx = np.where(unique_lons == lon)[0][0]
                imp_mat_3d[time_step, lat_idx, lon_idx] = value

        coords.update({'time': impact_object.date})
        data_vars.update({'imp_mat': (('time', 'latitude', 'longitude'), imp_mat_3d)})

    # Create xarray Dataset
    ds = xr.Dataset(data_vars, coords=coords)

    # Save to NetCDF
    ds.to_netcdf(filename)

# Example usage:
# save_impact_data_to_NetCDF(heat_waves_impact, 'impact_data_F.nc', include_eai_exp=True, include_imp_mat=True)
