import numpy as np
from scipy import signal
from scipy import fftpack
from scipy.optimize import fmin
from mtspec import mtspec
import matplotlib.pyplot as plt

# ----------------------------------------------------------------


def estimate_signal_to_noise(tx):
    number_cmps = np.size(tx[0, :])
    snr_array = np.zeros(number_cmps - 1)

    for i in range(number_cmps - 1):
        autocorr = np.correlate(tx[:, i], tx[:, i])
        autocorr.shape
        crosscorr = np.correlate(tx[:, i], tx[:, i + 1])
        snr_array[i] = np.power((abs(crosscorr) / abs(autocorr - crosscorr)), 0.5)
    snr_median = np.median(snr_array)
    snr_mean = np.mean(snr_array)
    snr_stdev = np.std(snr_array)

    return (snr_median, snr_mean, snr_stdev)


# ----------------------------------------------------------------


def make_amplitude_and_phase_spectra(tx, sampling, cmp_spacing):
    number_t = np.size(tx[:, 0])
    number_x = np.size(tx[0, :])

    ### Make fk spectrum of seismic data
    fk = np.real(fftpack.fft2(tx))
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
    f = fftpack.fftfreq(number_f, d=float(sampling))
    kx = fftpack.fftfreq(number_k, d=float(cmp_spacing))
    index = np.argsort(kx)
    kx = kx[index]
    index = np.argsort(f)
    f = f[index]

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


# ----------------------------------------------------------------


def make_direct_data_transform(tx, cmp_spacing, jres, K):
    def make_level_spectra(tx, cmp_spacing, jres, K):
        number_levels = np.size(tx[:, 0])
        for i in range(number_levels):
            level = tx[i, :]
            power_spectrum, kx, jackknife, _, _ = mtspec(
                data=level,
                delta=float(cmp_spacing),
                time_bandwidth=jres,
                number_of_tapers=K,
                nfft=None,
                statistics=True,
            )
            if i == 0:
                power_spectra_array = np.zeros(
                    [np.size(tx[:, 0]), np.size(power_spectrum)]
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
        tx, cmp_spacing, jres, K
    )
    ddt_slope_spectrum = np.power((2.0 * np.pi * kx), 2.0) * ddt_power_spectrum
    ddt_slope_spectrum_stdev = (
        np.power((2.0 * np.pi * kx), 2.0) * ddt_power_spectrum_stdev
    )
    #     Todo: when is correct point to average error?

    return (
        kx,
        ddt_power_spectrum,
        ddt_power_spectrum_stdev,
        ddt_slope_spectrum,
        ddt_slope_spectrum_stdev,
    )
