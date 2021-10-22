"""
This module is part of SeisMix.
Copyright 2021 Alex Dickinson.
Licensed under the GNU General Public License 3.0 (see LICENSE file).

miscellaneous_functions.py
~~~~~~~~~~~
Contains functions for:
- Rounding values to appropriate number of significant figures

Required dependencies:
- [`numpy`](http://numpy.org)
"""

import numpy as np


def round_to_sf(x, sigma_x):
    """
    Round values to number of significant figures determined by uncertainty.

    Parameters
    ----------
    x : float
        Value before rounding
    sigma_x : float
        Uncertainty on x before rounding

    Returns
    ----------
    x_round : float
        Value after rounding
    sigma_x_round : float
        Uncertainty on x_round after rounding
    """
    if x == 0 or sigma_x == 0:
        x_round = x
    else:
        sigma_x_power = -int(np.floor(np.log10(abs(sigma_x))))
        sigma_x_first_digit = round(sigma_x, sigma_x_power) * np.power(
            10, float(sigma_x_power)
        )

        if sigma_x_first_digit < 3:
            n = 2
        else:
            n = 1

        sigma_x_round = round(sigma_x, sigma_x_power + (n - 1))
        x_round = round(x, sigma_x_power + (n - 1))

    return (x_round, sigma_x_round)
