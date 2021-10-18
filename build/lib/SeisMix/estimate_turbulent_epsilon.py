
import numpy as np
from scipy import optimize

from python_mixing_functions import chi_square_analysis

import matplotlib.pyplot as plt


###### Functions for analysis of turbulent subrange ######
def fit_intercept_with_fixed_gradient(xdata, ydata, sigma_y, m, x0):
	### Find best-fit intercept with least squares. Equivalent to finding mean of ydata - ymodel.
#	print sigma_y
#	print np.mean(sigma_y)
	
#	def func(c, x, y):
#		return (y - m*x - c)
#	line = optimize.leastsq(func, x0, args=(xdata, ydata), full_output=1)
#	### Find error in best-fit intercept. Equivalent to finding standard error in the mean.
#	best_c = line[0]
#	print best_c
#	cov_matrix = line[1]
#	v11 = cov_matrix[0]
#	sigma_c = sigma_y*np.sqrt(v11)
#	print "root v11 * sigma_logphi = ", sigma_c
#	
#	### Find parameters using direct expressions:
#	print 'Direct'
#	c = np.mean(ydata[:] - m*xdata[:])
#	print c
#	sigma = sigma_y / np.sqrt(np.size(xdata))
#	print " direct expression sigma c = ", sigma
#	print " mean direct expression sigma c = ", np.mean(sigma)
	
	### Find parameters using direct expressions derived using Numerical Recipes pg. 782
	S = np.sum(1. / np.power(sigma_y, 2.))
	Sx = np.sum(xdata / np.power(sigma_y, 2.))
	Sy = np.sum(ydata / np.power(sigma_y, 2.))
	
	best_c = (Sy - m * Sx) / S
	sigma_c = np.power((1. / S), 0.5)
	print(best_c, sigma_c)
	
	
	
#	ax = plt.subplot(111)
#	ax.plot(xdata, ydata)
#	ax.plot(xdata, ydata + sigma_y)
#	ax.plot(xdata, ydata - sigma_y)
##	ax.plot(turb_fitted_line[:,0], turb_fitted_line[:,1])
##	ax.set_xlim([np.min(xdata)-0.1, np.max(xdata)+0.1])
##	ax.set_ylim([np.min(ydata)-0.1, np.max(ydata)+0.1])
##	plt.plot(noise_logkx, noise_logphi, 'g')
#	plt.show()
	
	
	return(best_c, sigma_c)

### Estimate logepsilon and logkt using the best fitting intercept, c, of the straight line with gradient 1/3.
def epsilon_kt_from_turb_intercept(c, gamma, CT, N_rads):
	logepsilon = 1.5*c - 1.5*np.log10(gamma) - 1.5*np.log10(CT) + 3.*np.log10(N_rads) - 3.5*np.log10(2.) - 2.*np.log10(np.pi)
	logkt = 1.5*c - 0.5*np.log10(gamma) - 1.5*np.log10(CT) + np.log10(N_rads) - 3.5*np.log10(2.) - 2.*np.log10(np.pi)
	return(logepsilon, logkt)

### Estimate error in logepsilon and logkt using error propagation formula.
def epsilon_kt_error_propagation(c2_turb, sigma_c2_turb, N_rads, sigma_N_rads):
	sigma_logepsilon = np.power( ( np.power((1.5*sigma_c2_turb), 2.) + np.power((3.*np.log10(np.exp(1.))*(sigma_N_rads/N_rads)), 2.) ), 0.5)
	sigma_logkt = np.power( ( np.power((1.5*sigma_c2_turb), 2.) + np.power((np.log10(np.exp(1.))*(sigma_N_rads/N_rads)), 2.) ), 0.5)
	return(sigma_logepsilon, sigma_logkt)


def estimate_epsilon_from_turbulent_subrange(logkx, logphi, sigma_logphi, N_rads, sigma_N_rads, gamma, CT):
	### Find best-fitting intercept for fixed gradient, m, of 1/3. x0 is initial guess of intercept.
	
	m = 1./3.
	x0 = 0.0
	c2_turb, sigma_c2_turb = fit_intercept_with_fixed_gradient(logkx, logphi, sigma_logphi, m, x0)
	
	### Do chi-squared analysis for straight line with fitted intercept. Degrees of freedom (dof) = M - 1, where M is number of data points in turbulent subrange.
	dof = np.size(logkx) - 1
	chi_sq_data_turb = chi_square_analysis.chi_squared(logkx, logphi, m, c2_turb, sigma_logphi)
	
	### Find coefficients for chi-squared well, i.e. equation chi_squared = a_well*c^2 + b_well*c + d_well, where c is fitted gradient. Valid only for when only the intercept (and not the gradient) is being fitted.
	a_well, b_well, d_well = chi_square_analysis.chi_squared_well_direct_expression(logkx, logphi, m, c2_turb, sigma_logphi)
	### Find values of c at which chi_squared is n times its minimum value. This function is valid only for when only the intercept (and not the gradient) is being fitted.
	n = 2.
	plus_c2, minus_c2 = chi_square_analysis.chi_squared_acceptance_levels(logkx, logphi, m, c2_turb, sigma_logphi, chi_sq_data_turb, n)
	
	### Find value of chi-squared for specified significance level, sig. sig = 0.95 is significance level of 95%.
	sig_level_turb = 0.67
	chi_sq_sig_turb = chi_square_analysis.chi_squared_significance(sig_level_turb, dof)
	
	### Estimate logepsilon and logkt from the best-fitted intercept, c. Estimate their standard errors using propagation formula.
	logepsilon, logkt = epsilon_kt_from_turb_intercept(c2_turb, gamma, CT, N_rads)
	sigma_logepsilon, sigma_logkt = epsilon_kt_error_propagation(c2_turb, sigma_c2_turb, N_rads, sigma_N_rads)
	
	print(logepsilon, sigma_logepsilon, "logepsilon, sigma_logepsilon")
	print(logkt, sigma_logkt, "logkt, sigma_logkt")
	
#	plt.show()
	
	return(c2_turb, sigma_c2_turb, logepsilon, sigma_logepsilon, logkt, sigma_logkt)



















