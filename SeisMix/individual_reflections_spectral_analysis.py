import numpy as np
from scipy import signal
from scipy import optimize
from scipy import stats
from scipy import interpolate
from mtspec import mtspec

import gsw

import matplotlib.pyplot as plt

from python_mixing_functions import compute_gm76_horizontal_spectra
from python_mixing_functions import compute_gk91_horizontal_spectra
from python_mixing_functions import gm_model_fitting


######------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###### Functions for spectral estimation ######
### Take spectra of midpoints using MTFT.
def calculate_spectra(midpoints, sound_speed_f, dx, jres_mtft_f, K_mtft_f):
    ### Detrend twts
    depths = signal.detrend(0.5 * sound_speed_f * midpoints[:, 1], type="linear")
    power_spectrum, kx, jackknife, _, _ = mtspec(
        data=depths,
        delta=float(dx),
        time_bandwidth=jres_mtft_f,
        number_of_tapers=K_mtft_f,
        nfft=None,
        statistics=True,
    )
    #     r = pymutt.mtft(depths, kind=ki, npi=jres_mtft, nwin=K_mtft, dt=dx)
    #     power_spectrum = r['power']
    #     kx = np.arange(r['nspec']) * r['df']
    slope_spectrum = power_spectrum * (np.power((2.0 * np.pi * kx), 2.0))
    return (kx, power_spectrum, slope_spectrum)


def compute_spectral_confidence_limits(sig_level, degfree):
    ### Expression from Percival and Walden, Chapter 6, Equation 258. Formulated for log base 10.
    p = (1.0 - sig_level) / 2.0
    Qvp = stats.chi2.ppf(p, degfree)
    Qvp1 = stats.chi2.ppf((1.0 - p), degfree)
    spectrum_confidence_limit_high = 0.4343 * np.log10(degfree / Qvp)
    spectrum_confidence_limit_low = 0.4343 * np.log10(degfree / Qvp1)
    confidence_interval_width = 0.4343 * np.log10(Qvp1 / Qvp)
    return (spectrum_confidence_limit_low, spectrum_confidence_limit_high)


### Estimate error in logamp using expressions of Percival, 1994 ######
def estimate_sigma_logamp(K_mtft, sig_level):
    ### Degrees of freedom are estimated as 2K_mtft.
    degfree = 2 * K_mtft
    ### Use two-sided significance level.
    p = (1.0 - sig_level) / 2.0
    Qvp = stats.chi2.ppf(p, degfree)
    Qvp1 = stats.chi2.ppf((1.0 - p), degfree)
    upper_bound = np.log10(degfree / Qvp)
    lower_bound = np.log10(degfree / Qvp1)
    # 	print upper_bound
    # 	print lower_bound
    # 	x = np.linspace(stats.chi2.ppf(0.01, degfree), stats.chi2.ppf(0.99, degfree), 100)
    # 	plt.plot(x, stats.chi2.pdf(x, degfree))
    # 	plt.plot([Qvp, Qvp], [np.min(stats.chi2.pdf(x, degfree)), np.max(stats.chi2.pdf(x, degfree))])
    # 	plt.plot([Qvp1, Qvp1], [np.min(stats.chi2.pdf(x, degfree)), np.max(stats.chi2.pdf(x, degfree))])
    # 	plt.show()
    # 	print "here"
    ### For now, take larger of two error bounds on logamp. This is not correct. TODO think about this. ###### Need to include effect of multiplication by (2pikx)^2.
    sigma_logamp = np.max(abs(upper_bound), abs(lower_bound))
    return sigma_logamp


### Include effects of local sound speed in estimated error of displacement wavenumber spectrum.
def estimate_sigma_logphi(
    spectrum_confidence_limit_low,
    spectrum_confidence_limit_high,
    sound_speed_sigma_sound_speed,
):
    sigma_logphi = np.power(
        (
            (
                np.power(sigma_logamp, 2.0)
                + np.power(
                    (2.0 * np.log10(np.exp(1.0)) * (sigma_sound_speed / sound_speed)),
                    2.0,
                )
            )
        ),
        0.5,
    )
    return sigma_logphi


######------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###### Functions for four-parameter model fitting ######
def inst_model(logkx, par):
    iw = par[3] * (logkx - par[0]) + par[1]
    turb = 1.0 / 3.0 * (logkx - par[0]) + par[1]
    noise = 2.0 * (logkx - (par[0] + par[2])) + (par[1] + 1.0 / 3.0 * par[2])
    return np.fmax(iw, np.fmax(turb, noise))


def misfit(par):
    return np.nansum(np.power((logphi - inst_model(logkx, par)), 2))


### Three-parameter model. Specify internal wave to have -1 gradient.
def inst_model_three_parameter(logkx, par):
    iw = -1.0 * (logkx - par[0]) + par[1]
    turb = 1.0 / 3.0 * (logkx - par[0]) + par[1]
    noise = 2.0 * (logkx - (par[0] + par[2])) + (par[1] + 1.0 / 3.0 * par[2])
    return np.fmax(iw, np.fmax(turb, noise))


def misfit_three_component(par):
    return np.nansum(np.power((logphi - inst_model_three_parameter(logkx, par)), 2))


######------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###### Functions for analysis of turbulent subrange ######
def fit_intercept_with_fixed_gradient(xdata, ydata, sigma_y, m, x0):
    ### Find best-fit intercept with least squares. Equivalent to finding mean of ydata - ymodel.
    def func(c, x, y):
        return y - m * x - c

    line = optimize.leastsq(func, x0, args=(xdata, ydata), full_output=1)
    ### Find error in best-fit intercept. Equivalent to finding standard error in the mean.
    best_c = line[0]
    # 	print best_c
    cov_matrix = line[1]
    v11 = cov_matrix[0]
    sigma_c = sigma_y * np.sqrt(v11)
    # 	print "root v11 * sigma_logphi = ", sigma_c

    ### Find parameters using direct expressions:
    # 	print 'Direct'
    # 	c = np.mean(ydata[:] - m*xdata[:])
    # 	print c
    # 	sigma = sigma_y / np.sqrt(np.size(xdata))
    # 	print " direct expression sigma c = ", sigma

    return (best_c, sigma_c)


### Estimate logepsilon and logkt using the best fitting intercept, c, of the straight line with gradient 1/3.
def epsilon_kt_from_turb_intercept(c, gamma, CT, N_rads):
    logepsilon = (
        1.5 * c
        - 1.5 * np.log10(gamma)
        - 1.5 * np.log10(CT)
        + 3.0 * np.log10(N_rads)
        - 3.5 * np.log10(2.0)
        - 2.0 * np.log10(np.pi)
    )
    logkt = (
        1.5 * c
        - 0.5 * np.log10(gamma)
        - 1.5 * np.log10(CT)
        + np.log10(N_rads)
        - 3.5 * np.log10(2.0)
        - 2.0 * np.log10(np.pi)
    )
    return (logepsilon, logkt)


### Estimate error in logepsilon and logkt using error propagation formula.
def epsilon_kt_error_propagation(c2_turb, sigma_c2_turb, N_rads, sigma_N_rads):
    sigma_logepsilon = np.power(
        (
            np.power((1.5 * sigma_c2_turb), 2.0)
            + np.power((3.0 * np.log10(np.exp(1.0)) * (sigma_N_rads / N_rads)), 2.0)
        ),
        0.5,
    )
    sigma_logkt = np.power(
        (
            np.power((1.5 * sigma_c2_turb), 2.0)
            + np.power((np.log10(np.exp(1.0)) * (sigma_N_rads / N_rads)), 2.0)
        ),
        0.5,
    )
    return (sigma_logepsilon, sigma_logkt)


### Compare to Matts calculation.
def matts_calc(par, gamma, CT, N_rads):
    top = pow(10.0, par[1]) * N_rads * N_rads
    bottom = 4.0 * np.pi * gamma * CT * (pow(10, par[0]) * 2 * np.pi) ** 0.333333333
    logepsilon = np.log10((top / bottom) ** (1.5))
    logkt = np.log10((gamma * (top / bottom) ** (1.5)) / (N_rads * N_rads))
    return (logepsilon, logkt)


######------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###### Functions for analysis of internal wave subrange ######
### Straight line fitting (with both fitted gradient and fitted intercept) for isolated internal wave subrange.
def fit_straight_line(xdata, ydata, sigma_y, x0):
    def func(params, R, ydata):
        return ydata - np.dot(R, params)

    ### Organise x data into response matrix format.
    R = np.column_stack((np.ones(np.size(xdata)), xdata))
    line = optimize.leastsq(func, x0, args=(R, ydata), full_output=1)
    c = line[0][0]
    sigma_c = sigma_y * np.sqrt(line[1][0, 0])
    m = line[0][1]
    sigma_m = sigma_y * np.sqrt(line[1][1, 1])
    return (m, sigma_m, c, sigma_c)


def estimate_spectral_variance(kx, spectrum, kx_min_cpm, kx_max_cpm):
    kx_min_arg = np.argmax(kx >= kx_min_cpm)
    kx_max_arg = np.argmin(kx <= kx_max_cpm)
    kx_min = kx[kx_min_arg]
    kx_max = kx[kx_max_arg]
    kx_integrand = kx[kx_min_arg:kx_max_arg]
    spectrum_integrand = spectrum[kx_min_arg:kx_max_arg]
    spectrum_variance = np.trapz(spectrum_integrand, kx_integrand)
    return (kx_integrand, spectrum_integrand, spectrum_variance)


def compute_potential_energy(variance, N_rads, rho_kg_m3):
    gpe = 0.5 * rho_kg_m3 * np.power(N_rads, 2.0) * variance
    return gpe


### Estimate F0, epsilon0 and K0 using Polzin fine-scale parametrisation.
def estimate_K0(E, N0_rads, b, jstar, N_GM_rads, f_GM, Rf, gm_option):
    ### Estimate reference turbulent production (F0) and dissipation (epsilon0) and diapycnal diffusivity (K0) using the high-wavenumber limit of the GM spectrum
    ### Expression for F is Equation 27 of Polzin et al. (2014). The factor of (b*N0_rads)^4 is to account for the different definitions of E (Polzin's dimensional versus the dimensionless definition used here).
    ### The factor of np.power((2.*np.pi), 2.) is to convert from cyclical wavenumbers (used here) to circular wavenumbers (used in Polzin et al.)
    kzstar = jstar / (2.0 * b) * (N_GM_rads / N0_rads)
    if gm_option == "gk91":
        F0 = (
            (6.0 * f_GM)
            / (10.0 * np.pi)
            * np.power((E / N0_rads), 2.0)
            * np.power(kzstar, 2.0)
            * np.arccosh(N_GM_rads / f_GM)
            * np.power((b * N0_rads), 4.0)
            * np.power((2.0 * np.pi), 2.0)
        )
        epsilon0 = (1.0 - Rf) * F0
        K0 = Rf * F0 / np.power(N_GM_rads, 2.0)
        ### TODO Check values
    elif gm_option == "gm76":
        F0 = (
            (6.0 * f_GM)
            / (10.0 * np.pi)
            * np.power((E / N0_rads), 2.0)
            * np.power((2.0 / np.pi), 2.0)
            * np.power(kzstar, 2.0)
            * np.arccosh(N_GM_rads / f_GM)
            * np.power((b * N0_rads), 4.0)
            * np.power((2.0 * np.pi), 2.0)
        )
        epsilon0 = (1.0 - Rf) * F0
        ### TODO Check values
        K0 = Rf * F0 / np.power(N_GM_rads, 2.0)
    else:
        print("gm_option not specified - exiting()")
        exit()
    return (epsilon0, K0)


def estimate_K(
    Rw,
    RwGM,
    N_rads,
    f,
    N_GM_rads,
    f_GM,
    epsilon0,
    K0,
    spectrum_variance,
    gm_spectrum_variance,
):
    shear_strain_term = (
        Rw * (Rw + 1.0) / (6.0 * np.power(2.0, 0.5) * np.power((Rw - 1.0), 0.5))
    )
    ### TODO Check appropriate definitions of N.
    f_term = (f * np.arccosh(N_rads / f)) / (f_GM * np.arccosh(N_GM_rads / f_GM))
    epsilon = (
        epsilon0
        * np.power((N_rads / N_GM_rads), 2.0)
        * np.power(spectrum_variance, 2.0)
        / np.power(gm_spectrum_variance, 2.0)
        * shear_strain_term
        * f_term
    )
    K = (
        K0
        * np.power(spectrum_variance, 2.0)
        / np.power(gm_spectrum_variance, 2.0)
        * shear_strain_term
        * f_term
    )
    spectrum_variance_ratio = spectrum_variance / gm_spectrum_variance
    # 	print K, "K"
    # 	print gamma * epsilon / np.power(N_rads, 2.), "gamma * epsilon / np.power(N_rads, 2.)"
    return (spectrum_variance_ratio, epsilon, K)


def estimate_standard_deviation_spectral_ratio(
    spectrum_integrand, reference_spectrum_integrand
):
    integrand_ratio = spectrum_integrand / reference_spectrum_integrand
    integrand_ratio_std = np.std(integrand_ratio)
    return integrand_ratio_std


def estimate_logK_sigma(
    integrand_ratio_std,
    spectrum_integrand,
    reference_spectrum_integrand,
    spectrum_variance,
    reference_spectrum_variance,
):
    variance_ratio = spectrum_variance / reference_spectrum_variance
    logK_sigma = 2.0 * 0.4343 * integrand_ratio_std / variance_ratio
    return logK_sigma


### Estimation of epsilon and kt from internal wave subrange using Gregg-Henyey parametrisation.
def epsilon_kt_from_gregg_henyey(iw_logkx, iw_logphi, N_rads, gamma):
    iw_kx = np.power(10.0, iw_logkx)
    iw_phi = np.power(10.0, iw_logphi)

    ### Set up parameters for Gregg-Henyey parametrisation.
    N_GM_cph = 3.0
    N_GM_rads = 2.0 * np.pi * N_GM_cph / 3600.0
    Rw_GM = 3.0
    # 	f_GM = 7.3e-5
    epsilon_0 = 7.2 * np.power(10.0, -10.0)

    #     print N_GM_rads, "N_GM_rads old"
    #     print N_rads, "N_rads old"
    #     print f_GM, "f_GM old"
    #     print f, "f_old"

    J_Rw = (
        ((1.0 + 1.0 / Rw) / (1.0 + 1.0 / Rw_GM))
        * np.power(((1.0 - Rw_GM) / (1.0 - Rw)), 0.5)
        * np.power((Rw / Rw_GM), 2.0)
    )
    #     print J_Rw, "J_Rw"
    L_f = (f * np.arccosh(N_rads / f)) / (f_GM * np.arccosh(N_GM_rads / f_GM))
    #     print L_f, "L_f"

    ### Calculation using Garrett-Munk spectrum derived from Klymak toolbox.
    interp = interpolate.interp1d(logkx_reference, logphi_reference_slope)
    logphi_interp = interp(iw_logkx)
    phi_reference_interp = np.power(10.0, logphi_interp)
    epsilon = (
        epsilon_0
        * np.power((N_rads / N_GM_rads), 2.0)
        * np.power(np.mean(iw_phi / phi_reference_interp), 2.0)
        * J_Rw
        * L_f
    )
    logepsilon = np.log10(epsilon)

    # 	epsilon2 = prefactor * np.power((N_rads/N_GM_rads), 2.) * np.power((np.mean((np.power(10., (np.log10(iw_phi) - np.log10(phi)))))), 2.)
    # 	print epsilon, epsilon2

    #     print np.power(np.mean(iw_phi/phi_reference_interp), 2.), "integrand old"

    kt = (
        epsilon_0
        * gamma
        * np.power((1.0 / N_GM_rads), 2.0)
        * np.power(np.mean(iw_phi / phi_reference_interp), 2.0)
        * J_Rw
        * L_f
    )
    logkt = np.log10(kt)

    # 	print np.size(iw_logkx)
    # 	print np.size(logphi_interp)
    #
    # 	plt.plot(iw_logkx, iw_logphi)
    # 	plt.plot(iw_logkx, logphi_interp)
    # 	plt.plot(logkx, logphi)
    # 	plt.show()

    ### Calculation for power law fit.
    # 	a0 = 0.001
    # 	N_GM_cph = 3.
    # 	N_GM_rads = 2.*np.pi*N_GM_cph/3600.
    #
    # 	GM_amp = a0 * np.power(iw_kx, GM_p)
    #
    # 	prefactor = 7.*np.power(10., -10.)
    #
    # 	epsilon = prefactor * np.power((N_rads/N_GM_rads), 2.) * np.power(np.mean(iw_amp/GM_amp), 2.)
    # 	logepsilon = np.log10(epsilon)
    #
    # 	kt = prefactor * gamma * np.power((1./N_GM_rads), 2.) * np.power(np.mean(iw_amp/GM_amp), 2.)
    # 	logkt = np.log10(kt)

    return (logepsilon, logkt)


######------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
###### Functions for chi-squared analysis. ######
### Chi-squared function for straight line fitting with constant error sigma_y in ydata
def chi_squared(xdata, ydata, m, c, sigma_y):
    chi_sq = np.sum(np.power(((ydata[:] - (m * xdata[:] + c)) / sigma_y), 2.0))
    return chi_sq


### Chi-squared value at sig% significance level for dof degrees of freedom.
def chi_squared_significance(sig, dof):
    chi_sq_sig = stats.chi2.ppf(0.95, dof)
    return chi_sq_sig


### Compute values of c at which chi-squared well reaches a value n times its minimum. ### This is for case when only c is varied and m is a specified constant value ###
def chi_squared_acceptance_levels(xdata, ydata, m, c, sigma_y, chi_sq_min, n):
    M = np.size(xdata)
    plus_c = (sigma_y / M) * (
        -1.0 * np.sum((ydata[:] - (m * xdata[:] + c)) / sigma_y)
        + np.sqrt(
            (np.power(np.sum((ydata[:] - (m * xdata[:] + c)) / sigma_y), 2.0))
            + M * (n - 1.0) * chi_sq_min
        )
    )
    minus_c = (sigma_y / M) * (
        -1.0 * np.sum((ydata[:] - (m * xdata[:] + c)) / sigma_y)
        - np.sqrt(
            (np.power(np.sum((ydata[:] - (m * xdata[:] + c)) / sigma_y), 2.0))
            + M * (n - 1.0) * chi_sq_min
        )
    )
    return (plus_c, minus_c)


### Calculate coefficients a, b and d for equation chi_squared = ac^2 + bc + d, where c is fitted gradient. These expressions are valid only for when only the intercept c is varied and the gradient m is a specified constant value.
def chi_squared_well_direct_expression(xdata, ydata, m, c, sigma_y):
    a = np.size(ydata) * np.power(sigma_y, -2.0)
    b = 2.0 * np.sum((m * xdata[:] - ydata[:]) * np.power(sigma_y, -2.0))
    d = np.sum(
        (
            np.power(ydata[:], 2.0)
            - 2.0 * m * ydata[:] * xdata[:]
            + np.power((m * xdata[:]), 2.0)
        )
        * np.power(sigma_y, -2.0)
    )
    return (a, b, d)


def iw_slope_uncertainty(logkx, logphi):
    best_model_iw = best_par[3] * (logkx - best_par[0]) + best_par[1]
    # 	line = np.polyfit(logkx, logphi, 1, full=False)
    ### Calculate root-mean-squared error in straight line.
    residual_squared = np.power((logphi[:] - best_model_iw[:]), 2.0)
    sigma_estimate = np.power(
        (np.sum(np.power(residual_squared, 2.0)) / (np.size(residual_squared) - 2.0)),
        0.5,
    )
    return sigma_estimate
