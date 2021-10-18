import numpy as np
from scipy import optimize

import matplotlib.pyplot as plt


def estimate_spectral_variance(kx, spectrum, kx_min_cpm, kx_max_cpm):
	
	kx_min_arg = np.argmax(kx >= kx_min_cpm)
	kx_max_arg = np.argmin(kx <= kx_max_cpm)
	
	kx_min = kx[kx_min_arg]
	kx_max = kx[kx_max_arg]
	kx_integrand = kx[kx_min_arg:kx_max_arg]
	
	spectrum_integrand = spectrum[kx_min_arg:kx_max_arg]
	
	spectrum_variance = np.trapz(spectrum_integrand, kx_integrand)
	
	return(kx_integrand, spectrum_integrand, spectrum_variance)


def compute_potential_energy(spectral_variance, N_rads):
	
	gpe = 0.5 * np.power(N_rads, 0.5) * spectral_variance
	
	return(gpe)
