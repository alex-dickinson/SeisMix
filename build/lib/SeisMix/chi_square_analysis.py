import numpy as np
from scipy import stats


###### Functions for chi-squared analysis. ######
### Chi-squared function for straight line fitting with constant error sigma_y in ydata 
def chi_squared(xdata, ydata, m, c, sigma_y):
	chi_sq = np.sum(np.power(((ydata[:] - (m*xdata[:] + c)) / sigma_y), 2.))
	return(chi_sq)

### Chi-squared value at sig% significance level for dof degrees of freedom.
def chi_squared_significance(sig, dof):
	chi_sq_sig = stats.chi2.ppf(0.95, dof)
	return(chi_sq_sig)

### Compute values of c at which chi-squared well reaches a value n times its minimum. ### This is for case when only c is varied and m is a specified constant value ###
def chi_squared_acceptance_levels(xdata, ydata, m, c, sigma_y, chi_sq_min, n):
	M = np.size(xdata)
	plus_c = (sigma_y / M) * (-1. * np.sum((ydata[:] - (m*xdata[:] + c)) / sigma_y) + np.sqrt( (np.power( np.sum((ydata[:] - (m*xdata[:] + c)) / sigma_y ), 2.)) + M*(n-1.)*chi_sq_min))
	minus_c = (sigma_y / M) * (-1. * np.sum((ydata[:] - (m*xdata[:] + c)) / sigma_y) - np.sqrt( (np.power( np.sum((ydata[:] - (m*xdata[:] + c)) / sigma_y ), 2.)) + M*(n-1.)*chi_sq_min))
	return(plus_c, minus_c)

### Calculate coefficients a, b and d for equation chi_squared = ac^2 + bc + d, where c is fitted gradient. These expressions are valid only for when only the intercept c is varied and the gradient m is a specified constant value.
def chi_squared_well_direct_expression(xdata, ydata, m, c, sigma_y):
	a = np.size(ydata) * np.power(sigma_y, -2.)
	b = 2.*np.sum((m*xdata[:] - ydata[:]) * np.power(sigma_y, -2.))
	d = np.sum((np.power(ydata[:], 2.) - 2.*m*ydata[:]*xdata[:] + np.power((m*xdata[:]), 2.)) * np.power(sigma_y, -2.))
	return(a, b, d)
