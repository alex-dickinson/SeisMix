import numpy as np
import scipy.fftpack as fftpackage
from scipy import signal
from scipy import stats

import matplotlib.pyplot as plt

#-----------------------------------------------------------------------------------------------------------------

##def estimate_sigma_logamp(K, sig_level):
##	### Degrees of freedom are estimated as 2K.
##	degfree = 2*K
##	### Use two-sided significance level.
##	p = (1.-sig_level) / 2.
##	Qvp = stats.chi2.ppf(p, degfree)
##	Qvp1 = stats.chi2.ppf((1.-p), degfree)
##	upper_bound = np.log10(degfree / Qvp)
##	lower_bound = np.log10(degfree / Qvp1)
###	print upper_bound
###	print lower_bound
###	x = np.linspace(stats.chi2.ppf(0.01, degfree), stats.chi2.ppf(0.99, degfree), 100)
###	plt.plot(x, stats.chi2.pdf(x, degfree))
###	plt.plot([Qvp, Qvp], [np.min(stats.chi2.pdf(x, degfree)), np.max(stats.chi2.pdf(x, degfree))])
###	plt.plot([Qvp1, Qvp1], [np.min(stats.chi2.pdf(x, degfree)), np.max(stats.chi2.pdf(x, degfree))])
###	plt.show()
###	print "here"
##	### For now, take larger of two error bounds on logamp. This is not correct. TODO think about this. ###### Need to include effect of multiplication by (2pikx)^2.
##	sigma_logamp = np.max(abs(upper_bound), abs(lower_bound))
##	return(sigma_logamp)

#-----------------------------------------------------------------------------------------------------------------

def make_welch_average_spectrum(temp_cmp_spacing_m, average_spectrum_number_points, temp_midpoints_dict, temp_twt_min_s, temp_twt_max_s):
    
    all_welch_spectra_dict = {}
    spectra_counter = 0

    for temp_midpoints in temp_midpoints_dict:      
        temp_cmps = temp_midpoints_dict[temp_midpoints][:,0]
        temp_twts = temp_midpoints_dict[temp_midpoints][:,1]
        temp_twt_mean = np.mean(temp_twts)

        if temp_twt_mean > temp_twt_min_s and temp_twt_mean < temp_twt_max_s:
            detrended_twts = signal.detrend(temp_twts, type='linear')
            fqs = 1. / temp_cmp_spacing_m

            ### TODO - change this depth conversion. Depth conversion using constant value of 750 m s^{-1} - could change this
            detrended_depths = 750. * detrended_twts
            average_spectrum_midpoints_length = average_spectrum_number_points
            
            ### Default window is Hann window
            kx, welch_spectrum = signal.welch(detrended_depths, fs=fqs, nperseg=average_spectrum_midpoints_length, scaling='density')

            ### Remove value at zero wavenumber
            kx = kx[1::]
            welch_spectrum = welch_spectrum[1::]
            
            ### Convert to slope spectrum and take log10 values
            slope_spectrum = np.power((2.*np.pi*kx), 2.) * welch_spectrum
            logkx = np.log10(kx)
            logphi = np.log10(slope_spectrum)
            all_welch_spectra_dict[temp_midpoints] = np.column_stack((kx, slope_spectrum, logkx, logphi))

            if spectra_counter == 0:
                all_detrended_slope_spectra_array = np.column_stack((kx, slope_spectrum))
                all_detrended_logphi_array = np.column_stack((logkx, logphi))
            else:
                all_detrended_slope_spectra_array = np.column_stack((all_detrended_slope_spectra_array, slope_spectrum))
                all_detrended_logphi_array = np.column_stack((all_detrended_logphi_array, logphi))

            spectra_counter = spectra_counter + 1

    ###### Average spectra
    ### Average in linear space and then transform to log-space
    average_slope_spectrum = np.mean(all_detrended_slope_spectra_array[:,1::], axis=1)
    std_slope_spectrum = np.std(all_detrended_slope_spectra_array[:,1::], axis=1)
    log10_mean_linear = np.log10(average_slope_spectrum)
    log10_mean_linear_std = np.log10(np.exp(1)) * std_slope_spectrum / average_slope_spectrum

    ### Average in log-space
    mean_log10 = np.mean(all_detrended_logphi_array[:,1::], axis=1)
    mean_log10_std = np.std(all_detrended_logphi_array[:,1::], axis=1)
    
    avespec_array = np.column_stack((kx, average_slope_spectrum, std_slope_spectrum, logkx, log10_mean_linear, log10_mean_linear_std, mean_log10, mean_log10_std))

    return(avespec_array)



