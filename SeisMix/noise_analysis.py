"""
This module is part of SeisMix.
Copyright 2021 Alex Dickinson.
Licensed under the GNU General Public License 3.0 (see LICENSE file).

noise_analysis.py
~~~~~~~~~~~
Contains functions for:
- Estimating signal-to-noise ratio for a seismic image
- Computing amplitude and phase spectra for a seismic image
- Computing the direct data transform of a seismic image

Required dependencies:
- [`mtspec`](https://krischer.github.io/mtspec/)
- [`numpy`](http://numpy.org)
- [`scipy`](https://scipy.org)
"""

import numpy as np
from scipy import signal
from scipy import fftpack
from scipy.optimize import fmin
from mtspec import mtspec
import matplotlib.pyplot as plt


def estimate_signal_to_noise(data):
    """
    Estimate signal-to-noise ratio, SNR, for seismic image.

    The SNR for two adjacent seismic traces is defined by

    sqrt( |c| / |a -c|)

    where c is maximum value of the cross-correlation between these traces, and a is the value of the zero-lag autocorrelation of the first trace.

    Parameters
    ----------
    data : numpy.ndarray
        2D numpy array listing seismic amplitude at every point in image

    Returns
    ----------
    snr_median : int
        median value of SNR for all pairs of adjacent traces within image
    snr_mean : int
        mean value of SNR for all pairs of adjacent traces within image
    snr_stdev : int
        standard deviation of SNR for all pairs of adjacent traces within image
    """

    number_cmps = np.size(data[0, :])
    snr_array = np.zeros(number_cmps - 1)

    for i in range(number_cmps - 1):
        autocorr = np.correlate(data[:, i], data[:, i])
        autocorr.shape
        crosscorr = np.correlate(data[:, i], data[:, i + 1])
        snr_array[i] = np.power((abs(crosscorr) / abs(autocorr - crosscorr)), 0.5)
    snr_median = np.median(snr_array)
    snr_mean = np.mean(snr_array)
    snr_stdev = np.std(snr_array)

    return (snr_median, snr_mean, snr_stdev)


def make_amplitude_and_phase_spectra(data, sampling_s, cmp_spacing_m):
    """
    Compute amplitude and phase spectra for seismic image.

    Parameters
    ----------
    data : numpy.ndarray
        2D numpy array listing seismic amplitude at every point in image
    sampling_s : float
        Sampling interval in seconds
    cmp_spacing_m : float
        Distance between adjacent common midpoints in metres

    Returns
    ----------
    f : numpy.ndarray
        1D array of frequency values at which spectrum computed
    kx : numpy.ndarray
        1D array of horizontal-wavenumber values at which spectrum computed
    amplitude_spectrum_sorted : numpy.ndarray
        2D array containing values of amplitude spectrum
    phase_spectrum_sorted : numpy.ndarray
        2D array containing values of phase spectrum
    kx_spectrum : numpy.ndarray
        1D array containing horizontal-wavenumber spectrum (i.e. amplitude spectrum summed over frequency axis)
    f_spectrum : numpy.ndarray
        1D array containing frequency spectrum (i.e. amplitude spectrum summed over horizontal-wavenumber axis)
    """

    ### Make fk spectrum of seismic data
    fk = np.real(fftpack.fft2(data))
    amplitude_spectrum = np.abs(fk)
    phase_spectrum = np.angle(fk)
    number_k = np.size(fk[0, :])
    number_f = np.size(fk[:, 0])

    ### Normalise amplitude spectrum
    amplitude_spectrum = amplitude_spectrum / np.max(amplitude_spectrum)

    ### Sort fk spectrum into conventional order
    amplitude_spectrum_sorted = np.fft.fftshift(amplitude_spectrum)
    phase_spectrum_sorted = np.fft.fftshift(phase_spectrum)

    ### Calculate kx values corresponding to wavenumber axis and f values corresponding to frequency axis
    f = fftpack.fftfreq(number_f, d=float(sampling_s))
    kx = fftpack.fftfreq(number_k, d=float(cmp_spacing_m))
    index = np.argsort(kx)
    kx = kx[index]
    index = np.argsort(f)
    f = f[index]

    ### Compute 1D wavenumber and frequency spectra by summing over opposite axes
    kx_spectrum = np.sum(amplitude_spectrum_sorted, axis=0)
    f_spectrum = np.sum(amplitude_spectrum_sorted, axis=1)

    return (
        f,
        kx,
        amplitude_spectrum_sorted,
        phase_spectrum_sorted,
        kx_spectrum,
        f_spectrum,
    )


def make_direct_data_transform(data, cmp_spacing_m, jres, Kres):
    """
    Compute direct data transform for seismic image (i.e. compute horizontal-wavenumber spectrum at each depth within image and sum all spectra). See Holbrook et al. (2013). Spectra are computed using multitaper Fourier transforms (package mtspcec; see Thomson (1982)).

    Parameters
    ----------
    data : numpy.ndarray
        2D numpy array listing seismic amplitude at every point in image
    cmp_spacing_m : float
        Distance between adjacent common midpoints in metres
    jres : int
        Parameter controlling wavenumber resolution of multitaper Fourier transform
    Kres : int
        Parameter controlling power resolution of multitaper Fourier transform

    Returns
    ----------
    kx : numpy.ndarray
        1D array of horizontal-wavenumber values at which direct data transform computed
    ddt_power_spectrum : numpy.ndarray
        1D array containing power of direct data transform at each horizontal wavenumber
    ddt_power_spectrum_stdev : numpy.ndarray
        1D array containing uncertainty of direct data transform at each horizontal wavenumber
    ddt_slope_spectrum : numpy.ndarray
        1D array containing direct data transform scaled by (2*pi*kx)^2. This transform emphasises the transition from the internal-wavenumber spectral subrange (scaled spectral slope <= 0) to the turbulent spectral subrange (scaled spectral subrange = +1/3)
    ddt_slope_spectrum_stdev : numpy.ndarray
        1D array containing uncertainty of scaled direct data transform at each horizontal wavenumber
    """

    def make_level_spectra(data, cmp_spacing_m, jres, Kres):
        number_levels = np.size(data[:, 0])
        for i in range(number_levels):
            level = data[i, :]
            power_spectrum, kx, jackknife, _, _ = mtspec(
                data=level,
                delta=float(cmp_spacing_m),
                time_bandwidth=jres,
                number_of_tapers=Kres,
                nfft=None,
                statistics=True,
            )
            if i == 0:
                power_spectra_array = np.zeros(
                    [np.size(data[:, 0]), np.size(power_spectrum)]
                )
            power_spectra_array[i, :] = power_spectrum
        ddt_power_spectrum = np.mean(power_spectra_array, axis=0)
        ddt_power_spectrum_stdev = np.std(power_spectra_array, axis=0)
        #         todo - error on DDT

        ### Remove first five points to remove spectral roll-off
        ddt_power_spectrum = ddt_power_spectrum[5::]
        ddt_power_spectrum_stdev = ddt_power_spectrum_stdev[5::]
        kx = kx[5::]
        return (kx, ddt_power_spectrum, ddt_power_spectrum_stdev)

    kx, ddt_power_spectrum, ddt_power_spectrum_stdev = make_level_spectra(
        data, cmp_spacing_m, jres, Kres
    )
    ddt_slope_spectrum = np.power((2.0 * np.pi * kx), 2.0) * ddt_power_spectrum
    ddt_slope_spectrum_stdev = (
        np.power((2.0 * np.pi * kx), 2.0) * ddt_power_spectrum_stdev
    )
    #     TODO: when is correct point to average error?

    return (
        kx,
        ddt_power_spectrum,
        ddt_power_spectrum_stdev,
        ddt_slope_spectrum,
        ddt_slope_spectrum_stdev,
    )
