"""
This file is part of CLIMADA.

Copyright (C) 2017 ETH Zurich, CLIMADA contributors listed in AUTHORS.

CLIMADA is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free
Software Foundation, version 3.

CLIMADA is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with CLIMADA. If not, see <https://www.gnu.org/licenses/>.
"""

# Standard library imports
import os
from datetime import datetime
import warnings

# Third-party library imports
import numpy as np
import xarray as xr
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# GIS and geospatial libraries
import rasterio
from osgeo import gdal, osr

# Cartopy for geographic plotting
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.feature



### Functions for 02 Warm-up and 03 Scenario ####


class ImpactWriteGeoTIFF:
    """
    A class to write impact data to a GeoTIFF file.

    Attributes
    ----------
    impact_object : Impact
        An object containing impact data, such as coordinates and values, to be written to a GeoTIFF file.

    Methods
    -------
    __init__(impact_object)
        Initialises the ImpactWriteGeoTIFF instance with the provided impact object.
    write_to_geotiff(filename)
        Writes the impact data to a GeoTIFF file at the specified location.
    """

    def __init__(self, impact_object):
        """
        Initializes the ImpactWriteGeoTIFF with an impact object.

        Parameters
        ----------
        impact_object : Impact
            An object containing impact data (e.g., coordinates and impact values) to be saved to a GeoTIFF file.
        """
        self.impact_object = impact_object

    def write_to_geotiff(self, filename):
        """
        Write the impact data to a GeoTIFF file.

        Parameters
        ----------
        filename : str
            Path and name of the file where the GeoTIFF data will be saved.

        Process
        -------
        - Extracts coordinates and impact values from the impact object.
        - Aggregates and averages values for each unique coordinate pair.
        - Constructs a GeoTIFF raster dataset with appropriate geospatial metadata.
        - Writes the data into the GeoTIFF format for spatial visualization and analysis.

        Raises
        ------
        ValueError
            Raised if the provided impact object does not contain the necessary data or if coordinates are invalid.
        """
        # Extract coordinates and impact values from the impact object
        coords = self.impact_object.coord_exp
        eai_exp_values = self.impact_object.eai_exp

        # Create a dictionary to hold summed eai_exp values and counts for each unique (lat, lon) pair
        eai_exp_dict = {}
        count_dict = {}

        # Sum and count eai_exp values for each unique coordinate
        for i in range(len(coords)):
            lat, lon = coords[i]
            key = (lat, lon)
            if key in eai_exp_dict:
                eai_exp_dict[key] += eai_exp_values[i]
                count_dict[key] += 1
            else:
                eai_exp_dict[key] = eai_exp_values[i]
                count_dict[key] = 1

        # Average the eai_exp values for each coordinate
        for key in eai_exp_dict:
            eai_exp_dict[key] /= count_dict[key]

        # Generate arrays for unique latitudes and longitudes
        unique_lats = np.unique(coords[:, 0])
        unique_lons = np.unique(coords[:, 1])

        # Initialize an array for reshaped eai_exp values
        eai_exp_reshaped = np.zeros((len(unique_lats), len(unique_lons)))

        # Reshape and assign eai_exp values to the array
        for i, lat in enumerate(unique_lats):
            for j, lon in enumerate(unique_lons):
                key = (lat, lon)
                if key in eai_exp_dict:
                    eai_exp_reshaped[i, j] = eai_exp_dict[key]

        # Flip the array along the latitude axis for correct geographical representation
        eai_exp_reshaped = np.flipud(eai_exp_reshaped)

        # Create a new GeoTIFF file
        driver = gdal.GetDriverByName("GTiff")
        dataset = driver.Create(
            filename, len(unique_lons), len(unique_lats), 1, gdal.GDT_Float32
        )

        # Set the geotransform for the dataset
        dataset.SetGeoTransform(
            [
                unique_lons.min(),  # Origin longitude
                (unique_lons.max() - unique_lons.min())
                / len(unique_lons),  # Pixel size in longitude
                0,
                unique_lats.max(),  # Origin latitude
                0,
                -(unique_lats.max() - unique_lats.min())
                / len(unique_lats),  # Pixel size in latitude (negative)
            ]
        )

        # Set the projection to WGS84
        srs = osr.SpatialReference()
        srs.ImportFromEPSG(4326)
        dataset.SetProjection(srs.ExportToWkt())

        # Write the eai_exp data to the raster band
        dataset.GetRasterBand(1).WriteArray(eai_exp_reshaped)
        dataset.FlushCache()  # Ensure data is written to disk


class ImpactReaderGeoTIFF:
    """
    A class to read and visualize impact data from a GeoTIFF file.

    Attributes
    ----------
    filename : str
        The full path to the GeoTIFF file containing impact data.

    Methods
    -------
    __init__(filename)
        Initializes the reader with the GeoTIFF file path.
    plot_geotiff(scale='normal')
        Visualizes the GeoTIFF data using cartographic features and optional scaling.
    _read_data(src)
        Reads raster data and spatial bounds from the GeoTIFF file.
    _create_plot()
        Creates a cartographic plot using matplotlib and Cartopy.
    _plot_data(ax, data, bounds, scale)
        Plots the GeoTIFF data on the map with specified scaling.
    _add_features(ax)
        Adds geographical features (e.g., borders, coastlines) to the plot.
    _add_labels_and_ticks(ax, bounds)
        Adds axis labels, titles, and grid ticks for better interpretability.
    """

    def __init__(self, filename):
        """
        Initializes the ImpactReaderGeoTIFF with the provided file path.

        Parameters
        ----------
        filename : str
            Path to the GeoTIFF file containing impact data.
        """
        self.filename = filename

    def plot_geotiff(self, scale="normal"):
        """
        Visualizes the impact data from the GeoTIFF file.

        Parameters
        ----------
        scale : str, optional
            The scaling to apply to the data ('normal' for linear or 'log' for logarithmic scaling). Defaults to 'normal'.

        Process
        -------
        - Reads data and spatial bounds from the GeoTIFF file.
        - Sets up a map projection using Cartopy.
        - Visualizes data with optional linear or logarithmic scaling.
        - Adds features like borders, coastlines, and gridlines for geographic context.

        Raises
        ------
        ValueError
            If the scale is not 'normal' or 'log'.
        """
        # Open the GeoTIFF file and read data
        with rasterio.open(self.filename) as src:
            data, bounds = self._read_data(src)
            fig, ax = self._create_plot()
            self._plot_data(ax, data, bounds, scale)
            self._add_features(ax)
            self._add_labels_and_ticks(ax, bounds)
            plt.show()

    def _read_data(self, src):
        """
        Reads data and spatial bounds from the GeoTIFF file.

        Parameters
        ----------
        src : rasterio.io.DatasetReader
            The open GeoTIFF file.

        Returns
        -------
        tuple
            A tuple containing:
            - data (numpy.ndarray): The raster data.
            - bounds (rasterio.coords.BoundingBox): The spatial bounds of the data.
        """
        data = src.read(1)
        bounds = src.bounds
        return data, bounds

    def _create_plot(self):
        """
        Creates a matplotlib plot with Cartopy for geographic visualization.

        Returns
        -------
        tuple
            A tuple containing:
            - fig (matplotlib.figure.Figure): The figure object.
            - ax (matplotlib.axes._axes.Axes): The axes object.
        """
        return plt.subplots(
            figsize=(10, 6), subplot_kw={"projection": ccrs.PlateCarree()}
        )

    def _plot_data(self, ax, data, bounds, scale):
        """
        Plots the GeoTIFF data on the map with specified scaling.

        Parameters
        ----------
        ax : matplotlib.axes._axes.Axes
            The axes object to plot on.
        data : numpy.ndarray
            The raster data to visualize.
        bounds : rasterio.coords.BoundingBox
            The spatial bounds of the data.
        scale : str
            The scaling to apply ('normal' for linear or 'log' for logarithmic scaling).

        Raises
        ------
        ValueError
            If the scale is not 'normal' or 'log'.
        """
        # Set the normalization and label based on the scale type
        if scale == "log":
            norm = LogNorm()
            label = "eai_exp (Log Scale)"
        else:
            norm = None
            label = "eai_exp"

        # Create a meshgrid for X and Y coordinates
        ny, nx = data.shape
        x = np.linspace(bounds.left, bounds.right, nx)
        y = np.linspace(bounds.top, bounds.bottom, ny)
        X, Y = np.meshgrid(x, y)

        # Plot the data using pcolormesh
        sc = ax.pcolormesh(X, Y, data, cmap="RdYlBu_r", norm=norm)
        cbar = plt.colorbar(sc, ax=ax, pad=0.05)  # Add a colorbar
        cbar.set_label(label)  # Label the colorbar

    def _add_features(self, ax):
        """
        Adds geographical features to the plot for better context.

        Parameters
        ----------
        ax : matplotlib.axes._axes.Axes
            The axes object to add features to.
        """
        ax.coastlines()  # Add coastlines
        ax.gridlines()  # Add gridlines
        ax.add_feature(cfeature.BORDERS, linestyle=":")  # Add country borders

    def _add_labels_and_ticks(self, ax, bounds):
        """
        Adds axis labels, titles, and ticks to the plot for better interpretability.

        Parameters
        ----------
        ax : matplotlib.axes._axes.Axes
            The axes object to add labels and ticks to.
        bounds : rasterio.coords.BoundingBox
            The spatial bounds of the data.
        """
        ax.set_xlabel("Longitude")
        ax.set_ylabel("Latitude")
        ax.set_title("Expected Annual Impact")

        # Set ticks for longitude and latitude
        num_ticks = 10
        x_ticks = np.linspace(bounds.left, bounds.right, num_ticks)
        y_ticks = np.linspace(bounds.bottom, bounds.top, num_ticks)

        ax.set_xticks(x_ticks)
        ax.set_yticks(y_ticks)
        ax.set_xticklabels([f"{x:.2f}" for x in x_ticks])
        ax.set_yticklabels([f"{y:.2f}" for y in y_ticks])


### Functions for 02 Warm-up and 03 Scenario ####
# save and read impact as netccdf


def parse_date(date_str):
    """
    Parses a date string using multiple formats to ensure compatibility.

    Parameters
    ----------
    date_str : str
        The date string to parse (e.g., '2024-12-08' or '2024-12-08T14:30:00').

    Returns
    -------
    str
        The parsed date in the format 'dd-mm-yyyy'.

    Raises
    ------
    ValueError
        If the date string does not match any of the expected formats.
    """
    for fmt in ("%Y-%m-%dT%H:%M:%S", "%Y-%m-%d"):
        try:
            return datetime.strptime(date_str, fmt).strftime("%d-%m-%Y")
        except ValueError:
            pass
    raise ValueError(f"Date format for '{date_str}' is not recognized")


def save_impact_data_to_NetCDF(
    impact_object,
    filename,
    include_eai_exp=True,
    include_imp_mat=True,
    log_scale_imp=False,
    log_scale_eai=False,
    time_attribute="event_name",
):
    """
    Saves impact data to a NetCDF file with optional logarithmic scaling.

    Parameters
    ----------
    impact_object : object
        An object containing environmental data (`eai_exp`) and impact matrix data (`imp_mat`).
    filename : str
        The file path and name for the output NetCDF file.
    include_eai_exp : bool, optional
        If True, includes environmental data (`eai_exp`) in the output. Default is True.
    include_imp_mat : bool, optional
        If True, includes impact matrix data (`imp_mat`) in the output. Default is True.
    log_scale_imp : bool, optional
        If True, applies logarithmic scaling to the impact matrix data. Default is False.
    log_scale_eai : bool, optional
        If True, applies logarithmic scaling to the environmental data. Default is False.
    time_attribute : str, optional
        The attribute in the `impact_object` containing time information for labeling the time dimension. Default is 'event_name'.

    Outputs
    -------
    NetCDF file
        A file at the specified location containing the processed data.

    Raises
    ------
    ValueError
        If the `time_attribute` is not found in the `impact_object`.
    TypeError
        If the time data in the `impact_object` is empty or None.

    Notes
    -----
    - If `log_scale_imp` or `log_scale_eai` is True, data is transformed using `log(value + 1)`.
    - Coordinate arrays (`latitude` and `longitude`) are derived from `impact_object.coord_exp`.
    """
    coords = {}
    data_vars = {}

    # Validate presence of time attribute
    if (
        not hasattr(impact_object, time_attribute)
        or getattr(impact_object, time_attribute) is None
    ):
        raise ValueError(
            f"Attribute '{time_attribute}' not found in impact_object or is None"
        )

    time_data = getattr(impact_object, time_attribute)
    if not time_data:
        raise TypeError(f"Time data in '{time_attribute}' is empty or None")

    # Process and transform eai_exp data
    if include_eai_exp and hasattr(impact_object, "eai_exp"):
        unique_lats, lat_inverse = np.unique(
            impact_object.coord_exp[:, 0], return_inverse=True
        )
        unique_lons, lon_inverse = np.unique(
            impact_object.coord_exp[:, 1], return_inverse=True
        )
        eai_exp_reshaped = np.full((len(unique_lats), len(unique_lons)), np.nan)
        for idx, value in enumerate(impact_object.eai_exp):
            transformed_value = np.log(value + 1) if log_scale_eai else value
            eai_exp_reshaped[lat_inverse[idx], lon_inverse[idx]] = transformed_value
        coords.update(
            {
                "latitude": ("latitude", unique_lats, {"units": "degrees_north"}),
                "longitude": ("longitude", unique_lons, {"units": "degrees_east"}),
            }
        )
        data_vars.update({"eai_exp": (("latitude", "longitude"), eai_exp_reshaped)})

    # Process and transform imp_mat data
    if include_imp_mat and hasattr(impact_object, "imp_mat"):
        dense_imp_mat = (
            impact_object.imp_mat.toarray()
            if isinstance(impact_object.imp_mat, csr_matrix)
            else impact_object.imp_mat
        )
        num_time_steps, num_spatial_points = dense_imp_mat.shape

        if "latitude" not in coords:
            unique_lats, lat_inverse = np.unique(
                impact_object.coord_exp[:, 0], return_inverse=True
            )
        if "longitude" not in coords:
            unique_lons, lon_inverse = np.unique(
                impact_object.coord_exp[:, 1], return_inverse=True
            )

        imp_mat_3d = np.full(
            (num_time_steps, len(unique_lats), len(unique_lons)), np.nan
        )
        for time_step in range(num_time_steps):
            for spatial_point in range(num_spatial_points):
                value = dense_imp_mat[time_step, spatial_point]
                transformed_value = np.log(value + 1) if log_scale_imp else value
                lat_idx, lon_idx = (
                    lat_inverse[spatial_point],
                    lon_inverse[spatial_point],
                )
                imp_mat_3d[time_step, lat_idx, lon_idx] = transformed_value

        if isinstance(time_data[0], str):
            time_data = [parse_date(date) for date in time_data]
        coords.update({"time": ("time", time_data)})
        data_vars.update({"imp_mat": (("time", "latitude", "longitude"), imp_mat_3d)})

    ds = xr.Dataset(data_vars, coords=coords)
    ds.attrs["crs"] = "EPSG:4326"

    if "imp_mat" in data_vars:
        ds["imp_mat"].encoding["_FillValue"] = -9999
        ds["imp_mat"] = ds["imp_mat"].astype("float32")

    compression_opts = dict(zlib=True, complevel=5)
    encoding = {var: compression_opts for var in ds.data_vars}
    ds.to_netcdf(filename, mode="w", encoding=encoding)

    absolute_path = os.path.abspath(filename)
    print(f"Data saved to NetCDF file at {absolute_path}")


# Suppress specific runtime warnings from shapely
warnings.filterwarnings(
    "ignore", "invalid value encountered in intersects", RuntimeWarning
)

### NetCDF reader ###
class ImpactReaderNetCDF:
    """
    A class to read and visualize environmental or impact matrix data from a NetCDF file.

    Attributes
    ----------
    filename : str
        The path to the NetCDF file to read.
    """

    def __init__(self, filename):
        """
        Initializes the ImpactReaderNetCDF with the specified NetCDF file.

        Parameters
        ----------
        filename : str
            The full path to the NetCDF file containing impact data.
        """
        self.filename = filename
        self.ds = None

    def read_netcdf(self):
        """
        Reads the data from the specified NetCDF file into an xarray dataset.

        Raises
        ------
        FileNotFoundError
            If the specified NetCDF file cannot be found.
        """
        self.ds = xr.open_dataset(self.filename)

    def visualize(self, data_type="eai_exp", time_step=0, scale="linear"):
        """
        Visualizes data from the NetCDF file with options for scaling and time-step selection.

        Parameters
        ----------
        data_type : str, optional
            The type of data to visualize ('eai_exp' for environmental data or 'imp_mat' for impact matrix). Default is 'eai_exp'.
        time_step : int, optional
            The time step to visualize (only applies to `imp_mat`). Default is 0.
        scale : str, optional
            The scale for visualization ('linear' for normal scaling or 'log' for logarithmic scaling). Default is 'linear'.

        Process
        -------
        - Extracts the specified data type (`eai_exp` or `imp_mat`) from the NetCDF file.
        - If `data_type` is 'imp_mat', selects the specified time step for visualization.
        - Creates a geographic plot using Cartopy and matplotlib.
        - Adds features like coastlines and borders for better context.

        Raises
        ------
        ValueError
            If the `data_type` is not 'eai_exp' or 'imp_mat'.
        IndexError
            If the specified `time_step` is out of range for 'imp_mat'.
        """
        if data_type not in ["eai_exp", "imp_mat"]:
            raise ValueError("data_type must be 'eai_exp' or 'imp_mat'")

        data = (
            self.ds[data_type].isel(time=time_step)
            if data_type == "imp_mat"
            else self.ds[data_type]
        )
        lat = self.ds["latitude"].values
        lon = self.ds["longitude"].values

        fig, ax = plt.subplots(
            figsize=(10, 6), subplot_kw={"projection": ccrs.PlateCarree()}
        )

        if scale == "log":
            norm = LogNorm()  # Avoid log(0) issues by ensuring min is > 0
            label = f"{data_type} (Log Scale)"
        else:
            norm = None
            label = f"{data_type}"

        sc = ax.pcolormesh(
            lon,
            lat,
            data.values,
            cmap="RdYlBu_r",
            transform=ccrs.PlateCarree(),
            norm=norm,
        )

        ax.set_extent(
            [lon.min(), lon.max(), lat.min(), lat.max()], crs=ccrs.PlateCarree()
        )
        ax.coastlines()
        gl = ax.gridlines(draw_labels=True)
        gl.top_labels = False
        gl.right_labels = False
        ax.add_feature(cartopy.feature.BORDERS, linestyle=":")

        plt.colorbar(sc, ax=ax, label=label, pad=0.05)

        # Set title based on data type and conditional inclusion of time step and date for 'imp_mat'
        if data_type == "imp_mat":
            date = str(self.ds["time"].values[time_step])[
                :10
            ]  # Extracting date string in YYYY-MM-DD format
            title = f"Visualization of {data_type} - Time Step {time_step} ({date})"
        else:
            title = f"Visualization of {data_type}"

        plt.title(title)
        plt.show()


### Functions for 04 Uncertainty and Sensitivity ####
# save and read impact as NetCDF of scenario data and create plots
# Ignore specific warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
# warnings.filterwarnings("ignore", category=mpl.MatplotlibDeprecationWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning, module="shapely")


class ImpactReaderNetCDF_Scenario:
    """
    A class designed to read and visualize scenario-based impact data from NetCDF files.

    Attributes
    ----------
    filename : str
        The full path to the NetCDF file containing scenario impact data.
    """

    def __init__(self, filename):
        """
        Initializes the ImpactReaderNetCDF_Scenario with the specified NetCDF file.

        Parameters
        ----------
        filename : str
            The full path to the NetCDF file to be read.
        """
        self.filename = filename
        self.ds = None  # Will hold the xarray dataset after reading the NetCDF file
        self.lat = None  # Will hold the latitude coordinates
        self.lon = None  # Will hold the longitude coordinates

    def read_netcdf(self):
        """
        Reads the specified NetCDF file and extracts latitude and longitude coordinates.

        Raises
        ------
        KeyError
            If latitude or longitude coordinates are not found in the dataset.
        """
        self.ds = xr.open_dataset(self.filename)
        # Attempt to find the latitude and longitude coordinates with common naming conventions
        lat_names = ["lat", "latitude"]
        lon_names = ["lon", "longitude"]
        for name in lat_names:
            if name in self.ds.coords:
                self.lat = self.ds[name].values
                break
        for name in lon_names:
            if name in self.ds.coords:
                self.lon = self.ds[name].values
                break
        if self.lat is None or self.lon is None:
            raise KeyError("Latitude or longitude coordinate not found in dataset.")

    def visualize_eai_exp(self, scale="normal"):
        """
        Visualizes the Expected Annual Impact Exposure (eai_exp) data.

        Parameters
        ----------
        scale : str, optional
            Scaling for visualization ('normal' for linear or 'log' for logarithmic). Default is 'normal'.

        Raises
        ------
        ValueError
            If the `eai_exp` data is not available in the dataset.
        """
        data = self.ds["eai_exp"].values
        self._visualize(data, "eai_exp", scale)

    def visualize_imp_mat(self, scale="normal"):
        """
        Visualizes the Impact Matrix (imp_mat) data from the NetCDF file.

        Parameters
        ----------
        scale : str, optional
            Scaling for visualization ('normal' for linear or 'log' for logarithmic). Default is 'normal'.

        Process
        -------
        - Averages the impact matrix data over time if a time dimension exists.
        - Calls the internal `_visualize` method for plotting.

        Raises
        ------
        ValueError
            If the `imp_mat` data is not available in the dataset.
        """
        if "time" in self.ds["imp_mat"].dims:
            data = self.ds["imp_mat"].isel(time=0).values
        else:
            data = self.ds["imp_mat"].values
        self._visualize(data, "imp_mat", scale)

    def _visualize(self, data, title, scale):
        """
        A helper method for visualizing impact data with linear or logarithmic scaling.

        Parameters
        ----------
        data : numpy.ndarray
            The data array to visualize.
        title : str
            The title for the plot (e.g., 'eai_exp' or 'imp_mat').
        scale : str
            Scaling for visualization ('normal' or 'log').

        Raises
        ------
        ValueError
            If `scale` is not 'normal' or 'log'.
        """
        fig, ax = plt.subplots(
            figsize=(10, 6), subplot_kw={"projection": ccrs.PlateCarree()}
        )
        norm = LogNorm() if scale == "log" else None
        label = f'{title} ({"Log Scale" if scale == "log" else "Linear Scale"})'
        sc = ax.pcolormesh(
            self.lon,
            self.lat,
            data,
            cmap="RdYlBu_r",
            transform=ccrs.PlateCarree(),
            norm=norm,
        )
        ax.set_extent(
            [self.lon.min(), self.lon.max(), self.lat.min(), self.lat.max()],
            crs=ccrs.PlateCarree(),
        )
        ax.coastlines()
        gl = ax.gridlines(draw_labels=True)
        gl.top_labels = False
        gl.right_labels = False
        ax.add_feature(cartopy.feature.BORDERS, linestyle=":")
        plt.colorbar(sc, ax=ax, label=label, pad=0.05)
        plt.title(f"{title} Visualization")
        plt.show()
