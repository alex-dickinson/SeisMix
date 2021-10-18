import numpy as np
from scipy import optimize

import matplotlib.pyplot as plt





### Power spectra
### Fit gm power spectral form using least squares
def gm75_residual_power(vars, x, y, sigma_y):
	
	C = vars[0]
	s = vars[1]
	t = vars[2]
	kstar = vars[3]
	
	kx = np.power(10., x)
	gm75_model = np.log10(C) - np.log10(kstar) - (t / s) * np.log10(1. + np.power(( kx / kstar), s))
	
	### Define residuals weighted by error in logphi
	residual = (y - gm75_model) / sigma_y
	
	return(residual)

def fit_gm75_kx_power_parametrisation_leastsq(logkx, logphi, sigma_logphi, par):
	
#	print par, "par"
	
	C_initial = par[0]
	s_initial = par[1]
	t_initial = par[2]
	kstar_initial = par[3]
	
	kx = np.power(10., logkx)
	gm75_model_initial = np.log10(C_initial) - np.log10(kstar_initial) - (t_initial / s_initial) * np.log10(1. + np.power(( kx / kstar_initial), s_initial))
	
	vars = par
#	
#	plt.plot(logkx, logphi)
#	plt.plot(logkx, gm75_model_initial)
#	plt.show()
	
	best_par = optimize.leastsq(gm75_residual_power, par, args=(logkx, logphi, sigma_logphi))
	
#	print best_par
#	print best_par[0]
#	print best_par[1]
#	
	C_fitted = best_par[0][0]
	s_fitted = best_par[0][1]
	t_fitted = best_par[0][2]
	kstar_fitted = best_par[0][3]
	
#	print best_par, "best_par"
	
#	print C_fitted, "C_fitted"
#	print s_fitted, "s_fitted"
#	print t_fitted, "t_fitted"
#	print kstar_fitted, "kstar_fitted"
	
#	print "leastsquares"
	
	gm75_model_fitted = np.log10(C_fitted) - np.log10(kstar_fitted) - (t_fitted / s_fitted) * np.log10(1. + np.power(( kx / kstar_fitted), s_fitted))
	
#	plt.plot(logkx, logphi, label='logpower')
#	plt.plot(logkx, gm75_model_initial, label='gm initial')
#	plt.plot(logkx, gm75_model_fitted, label='gm_fitted')
#	plt.legend()
#	plt.show()
#	
#	print best_par
	return(C_fitted, s_fitted, t_fitted, kstar_fitted)

def make_fitted_power_spectrum(logkx, C_fitted, s_fitted, t_fitted, kstar_fitted):
	
	kx = np.power(10., logkx)
	gm75_model_fitted = np.log10(C_fitted) - np.log10(kstar_fitted) - (t_fitted / s_fitted) * np.log10(1. + np.power(( kx / kstar_fitted), s_fitted))
	
#	plt.plot(logkx, gm75_model_fitted)
#	plt.show()
	
	return(gm75_model_fitted)



### Slope spectra
### Fit gm spectral form using least squares
def gm75_residual(vars, x, y, sigma_y):
	
	C = vars[0]
	s = vars[1]
	t = vars[2]
	kstar = vars[3]
	
	kx = np.power(10., x)
	gm75_model = np.log10(C) - np.log10(kstar) + 2.*x - (t / s) * np.log10(1. + np.power(( kx / kstar), s))
	
	### Define residuals weighted by error in logphi
	residual = (y - gm75_model) / sigma_y
	
	return(residual)

def fit_gm75_kx_parametrisation_leastsq(logkx, logphi, sigma_logphi, par):
	
#	print par, "par"
	
	C_initial = par[0]
	s_initial = par[1]
	t_initial = par[2]
	kstar_initial = par[3]
	
	kx = np.power(10., logkx)
	gm75_model_initial = np.log10(C_initial) - np.log10(kstar_initial) + 2.*logkx - (t_initial / s_initial) * np.log10(1. + np.power(( kx / kstar_initial), s_initial))
	
	vars = par
#	
#	plt.plot(logkx, logphi)
#	plt.plot(logkx, gm75_model_initial)
#	plt.show()
	
	best_par = optimize.leastsq(gm75_residual, par, args=(logkx, logphi, sigma_logphi))
	
#	print best_par
#	print best_par[0]
#	print best_par[1]
#	
	C_fitted = best_par[0][0]
	s_fitted = best_par[0][1]
	t_fitted = best_par[0][2]
	kstar_fitted = best_par[0][3]
	
#	print best_par, "best_par"
	
#	print C_fitted, "C_fitted"
#	print s_fitted, "s_fitted"
#	print t_fitted, "t_fitted"
#	print kstar_fitted, "kstar_fitted"
	
#	print "leastsquares"
	
	gm75_model_fitted = np.log10(C_fitted) - np.log10(kstar_fitted) + 2.*logkx - (t_fitted / s_fitted) * np.log10(1. + np.power(( kx / kstar_fitted), s_fitted))
	
#	plt.plot(logkx, logphi)
#	plt.plot(logkx, gm75_model_initial)
#	plt.plot(logkx, gm75_model_fitted)
#	plt.show()
	
#	print best_par
	return(C_fitted, s_fitted, t_fitted, kstar_fitted)

def make_fitted_spectrum(logkx, C_fitted, s_fitted, t_fitted, kstar_fitted):
	
	kx = np.power(10., logkx)
	gm75_model_fitted = np.log10(C_fitted) - np.log10(kstar_fitted) + 2.*logkx - (t_fitted / s_fitted) * np.log10(1. + np.power(( kx / kstar_fitted), s_fitted))
	
#	plt.plot(logkx, np.log10(gm75_model_fitted))
#	plt.show()
	
	return(gm75_model_fitted)








