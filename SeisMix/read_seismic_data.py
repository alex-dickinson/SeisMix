"""
This module is part of SeisMix.
Copyright 2021 Alex Dickinson.
Licensed under the GNU General Public License 3.0 (see LICENSE file).

read_seismic_data.py
~~~~~~~~~~~
Contains functions for:
- Importing seismic images in SEG-Y format
- Selecting portions of imported seismic images

Required dependencies:
- [`numpy`](http://numpy.org)
- [`segyio`](https://segyio.readthedocs.io/en/1.5.3/index.html) 
"""

import numpy as np
import segyio


def read_segy(filename):
    """
    Convert seismic image from SEG-Y format to numpy objects.

    Parameters
    ----------
    filename : str
            Path to the input SEG-Y file

    Returns
    ----------
    cmp_array : numpy.ndarray
            1D numpy array listing all common midpoint numbers
    cmp_grid : numpy.ndarray
            2D numpy array listing common midpoint number at every point in seismic image. Shape is (np.size(twt_s_array), np.size(cmp_array))
    twt_s_array : numpy.ndarray
            1D numpy array listing all two-way travel times in seconds
    twt_s_grid : numpy.ndarray
            2D numpy array listing two-way travel time (in seconds) at every point in seismic image. Shape is (np.size(twt_s_array), np.size(cmp_array))
    data : numpy.ndarray
            2D numpy array listing seismic amplitude at every point in seismic image. Shape is (np.size(twt_s_array), np.size(cmp_array))
    """
    with segyio.open(filename, ignore_geometry=True) as f:
        # Get basic attributes
        n_traces = f.tracecount
        sampling_ms = segyio.tools.dt(f) / 1000
        sampling_s = segyio.tools.dt(f) / 1e6
        n_samples = f.samples.size
        twt_ms_array = f.samples
        twt_s_array = twt_ms_array / 1e3
        data = f.trace.raw[:]  # Get all data into memory (could cause on big files)
        # Load headers
        bin_headers = f.bin

        # To do - extract cmp_min from segy
        cmp_min = 1
        cmp_array = np.arange(cmp_min, n_traces + 1, 1)
        data = data.T
        cmp_grid, twt_s_grid = np.meshgrid(cmp_array, twt_s_array)

    return (cmp_array, cmp_grid, twt_s_array, twt_s_grid, data)


def select_image_region(
    cmp_array,
    cmp_grid,
    cmp_select_min,
    cmp_select_max,
    twt_s_array,
    twt_s_grid,
    twt_select_s_min,
    twt_select_s_max,
    data,
):
    """
    Select region within seismic image.

    Parameters
    ----------
    cmp_array : numpy.ndarray
            1D numpy array listing all common midpoint numbers
    cmp_grid : numpy.ndarray
            2D numpy array listing common midpoint number at every point in seismic image. Shape is (np.size(twt_s_array), np.size(cmp_array))
    cmp_select_min : int
            Minimum common midpoint value defining edge of region to be selected
    cmp_select_max : int
            Maximum common midpoint value defining edge of region to be selected
    twt_s_array : numpy.ndarray
            1D numpy array listing all two-way travel times in seconds
    twt_s_grid : numpy.ndarray
            2D numpy array listing two-way travel time (in seconds) at every point in seismic image. Shape is (np.size(twt_s_array), np.size(cmp_array))
    twt_select_s_min : int
            Minimum value of two-way travel time (in seconds) defining edge of region to be selected
    twt_select_s_max : int
            Maximum value of two-way travel time (in seconds) defining edge of region to be selected
    data : numpy.ndarray
            2D numpy array listing seismic amplitude at every point in seismic image. Shape is (np.size(twt_s_array), np.size(cmp_array))

    Returns
    ----------
    cmp_array_select : numpy.ndarray
            1D numpy array listing all common midpoint numbers within selected region
    cmp_grid_select : numpy.ndarray
            2D numpy array listing common midpoint number at every point within selected region. Shape is (np.size(twt_s_array_select), np.size(cmp_array_select))
    twt_s_array_select : numpy.ndarray
            1D numpy array listing all two-way travel times (in seconds) within selected region
    twt_s_grid_select : numpy.ndarray
            2D numpy array listing two-way travel time (in seconds) at every point within selected region. Shape is (np.size(twt_s_array_select), np.size(cmp_array_select))
    data_select : numpy.ndarray
            2D numpy array listing seismic amplitude at every point within selected region. Shape is (np.size(twt_s_array_select), np.size(cmp_array_select))
    """
    cmp_select_min_index = np.argmax(cmp_array > cmp_select_min) - 1
    cmp_select_max_index = np.argmin(cmp_array < cmp_select_max) + 1
    twt_select_min_index = np.argmax(twt_s_array > twt_select_s_min) - 1
    twt_select_max_index = np.argmin(twt_s_array < twt_select_s_max) + 1
    data_select = data[
        twt_select_min_index:twt_select_max_index,
        cmp_select_min_index:cmp_select_max_index,
    ]
    cmp_grid_select = cmp_grid[
        twt_select_min_index:twt_select_max_index,
        cmp_select_min_index:cmp_select_max_index,
    ]
    cmp_array_select = cmp_array[cmp_select_min_index:cmp_select_max_index]
    twt_s_grid_select = twt_s_grid[
        twt_select_min_index:twt_select_max_index,
        cmp_select_min_index:cmp_select_max_index,
    ]
    twt_s_array_select = twt_s_array[twt_select_min_index:twt_select_max_index]

    return (
        cmp_array_select,
        cmp_grid_select,
        twt_s_array_select,
        twt_s_grid_select,
        data_select,
    )
