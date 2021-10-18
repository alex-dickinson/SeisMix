import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt


def estimate_spectral_variance(kx, spectrum, kx_min_cpm, kx_max_cpm):

    kx_min_arg = np.argmax(kx >= kx_min_cpm)
    kx_max_arg = np.argmin(kx <= kx_max_cpm)

    kx_min = kx[kx_min_arg]
    kx_max = kx[kx_max_arg]
    kx_integrand = kx[kx_min_arg:kx_max_arg]

    spectrum_integrand = spectrum[kx_min_arg:kx_max_arg]

    spectrum_variance = np.trapz(spectrum_integrand, kx_integrand)

    return (kx_integrand, spectrum_integrand, spectrum_variance)


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


####### Estimation of epsilon and kt from internal wave subrange using Gregg-Henyey parametrisation.
####def epsilon_kt_from_gregg_henyey(iw_logkx, iw_logphi, N_rads, gamma):
####	iw_kx = np.power(10., iw_logkx)
####	iw_phi = np.power(10., iw_logphi)
####
####	### Set up parameters for Gregg-Henyey parametrisation.
####	N_GM_cph = 3.
####	N_GM_rads = 2.*np.pi*N_GM_cph/3600.
####	Rw_GM = 3.
#####	f_GM = 7.3e-5
####	epsilon_0 = 7.2*np.power(10., -10.)
####
####	print N_GM_rads, "N_GM_rads old"
####	print N_rads, "N_rads old"
####	print f_GM, "f_GM old"
####	print f, "f_old"
####
####	J_Rw = ((1. + 1./Rw)/(1. + 1./Rw_GM)) * np.power(((1. - Rw_GM)/(1. - Rw)), 0.5) * np.power((Rw / Rw_GM), 2.)
####	print J_Rw, "J_Rw"
####	L_f = (f*np.arccosh(N_rads/f)) / (f_GM*np.arccosh(N_GM_rads/f_GM))
####	print L_f, "L_f"
####
####	### Calculation using Garrett-Munk spectrum derived from Klymak toolbox.
####	interp = interpolate.interp1d(logkx_reference, logphi_reference)
####	logphi_interp = interp(iw_logkx)
####	phi_reference_interp = np.power(10., logphi_interp)
####	epsilon = epsilon_0 * np.power((N_rads/N_GM_rads), 2.) * np.power(np.mean(iw_phi/phi_reference_interp), 2.) * J_Rw * L_f
####	logepsilon = np.log10(epsilon)
####
#####	epsilon2 = prefactor * np.power((N_rads/N_GM_rads), 2.) * np.power((np.mean((np.power(10., (np.log10(iw_phi) - np.log10(phi)))))), 2.)
#####	print epsilon, epsilon2
####
####	print np.power(np.mean(iw_phi/phi_reference_interp), 2.), "integrand old"
####
####	kt = epsilon_0 * gamma * np.power((1./N_GM_rads), 2.) * np.power(np.mean(iw_phi/phi_reference_interp), 2.) * J_Rw * L_f
####	logkt = np.log10(kt)
####
#####	print np.size(iw_logkx)
#####	print np.size(logphi_interp)
#####
#####	plt.plot(iw_logkx, iw_logphi)
#####	plt.plot(iw_logkx, logphi_interp)
#####	plt.plot(logkx, logphi)
#####	plt.show()
####
####
####
####
####	### Calculation for power law fit.
#####	a0 = 0.001
#####	N_GM_cph = 3.
#####	N_GM_rads = 2.*np.pi*N_GM_cph/3600.
#####
#####	GM_amp = a0 * np.power(iw_kx, GM_p)
#####
#####	prefactor = 7.*np.power(10., -10.)
#####
#####	epsilon = prefactor * np.power((N_rads/N_GM_rads), 2.) * np.power(np.mean(iw_amp/GM_amp), 2.)
#####	logepsilon = np.log10(epsilon)
#####
#####	kt = prefactor * gamma * np.power((1./N_GM_rads), 2.) * np.power(np.mean(iw_amp/GM_amp), 2.)
#####	logkt = np.log10(kt)
####
####	return(logepsilon, logkt)


### TODO Clear up this script and check through


def estimate_epsilon_from_internal_wave_subrange(
    all_logkx,
    all_logphi,
    all_sigma_phi,
    all_sigma_logphi,
    iw_logkx,
    iw_logphi,
    iw_sigma_logphi,
    iw_kx_lower_limit_integration,
    iw_kx_upper_limit_integration,
    logkx_reference,
    logphi_reference,
    E,
    N0_rads,
    b,
    jstar,
    N_GM_rads,
    Rf,
    Rw,
    RwGM,
    N_rads,
    f,
    f_GM,
    gm_option,
):

    all_kx = np.power(10.0, all_logkx)
    all_phi = np.power(10.0, all_logphi)
    all_sigma_phi_from_logphi = np.power(10.0, all_sigma_logphi)
    iw_kx = np.power(10.0, iw_logkx)
    iw_kx_min = iw_kx[0]
    iw_kx_max = iw_kx[-1]
    iw_phi = np.power(10.0, iw_logphi)

    all_phi_low = all_phi - all_sigma_phi
    all_phi_high = all_phi + all_sigma_phi
    all_logphi_low = all_logphi - all_sigma_logphi
    all_logphi_high = all_logphi + all_sigma_logphi

    all_phi_from_logphi_low = np.power(10.0, all_logphi_low)
    all_phi_from_logphi_high = np.power(10.0, all_logphi_high)

    # 	print all_phi
    # 	print all_sigma_phi
    # 	print all_sigma_phi_from_logphi
    # 	print all_phi_low
    # 	print all_phi_high
    #
    all_phi_variance = np.power(all_sigma_phi, 2.0)
    # 	print all_phi_variance

    # 	print all_logphi_low
    # 	print all_logphi_high
    # 	print all_phi_low
    # 	print all_phi_high
    # 	print all_phi_from_logphi_low
    # 	print all_phi_from_logphi_high

    # 	plt.plot(all_logkx, all_logphi)
    # 	plt.plot(iw_logkx, iw_logphi)
    # 	plt.plot(all_logkx, all_logphi_low)
    # 	plt.plot(all_logkx, all_logphi_high)
    # 	plt.show()

    (
        iw_kx_integrand,
        iw_spectrum_integrand,
        iw_spectrum_variance,
    ) = estimate_spectral_variance(
        all_kx, all_phi, iw_kx_lower_limit_integration, iw_kx_max
    )
    # 	iw_kx_integrand_low, iw_spectrum_integrand_low, iw_spectrum_variance_low = estimate_spectral_variance(all_kx, all_phi_low, iw_kx_lower_limit_integration, iw_kx_max)
    # 	iw_kx_integrand_high, iw_spectrum_integrand_high, iw_spectrum_variance_high = estimate_spectral_variance(all_kx, all_phi_high, iw_kx_lower_limit_integration, iw_kx_max)
    (
        iw_kx_integrand_low,
        iw_spectrum_integrand_low,
        iw_spectrum_variance_low,
    ) = estimate_spectral_variance(
        all_kx, all_phi_from_logphi_low, iw_kx_lower_limit_integration, iw_kx_max
    )
    (
        iw_kx_integrand_high,
        iw_spectrum_integrand_high,
        iw_spectrum_variance_high,
    ) = estimate_spectral_variance(
        all_kx, all_phi_from_logphi_high, iw_kx_lower_limit_integration, iw_kx_max
    )

    ### Integrate variance to quantify error
    # 	iw_variance_variance = estimate_error_in_variance(all_kx, all_phi_variance, iw_kx_lower_limit_integration, iw_kx_max)

    # 	print iw_kx_integrand

    kx_min_arg = np.argmax(all_kx >= iw_kx_lower_limit_integration)
    kx_max_arg = np.argmin(all_kx <= iw_kx_max)
    iw_variance_variance = np.sum(all_phi_variance[kx_min_arg:kx_max_arg]) * np.power(
        (iw_kx_integrand[1] - iw_kx_integrand[0]), 2.0
    )
    iw_variance_error = np.power(iw_variance_variance, 0.5)

    # 	iw_kx_integrand_variance, iw_spectrum_integrand_variance, iw_spectrum_variance_variance = estimate_spectral_variance(all_kx, all_phi_variance, iw_kx_lower_limit_integration, iw_kx_max)
    # 	iw_kx_integrand_spacing = iw_kx_integrand_variance[1::] - iw_kx_integrand_variance[0:-1]
    # 	iw_spectrum_variance_actual = iw_kx_integrand_spacing * iw_spectrum_variance_variance[0:-1]
    # 	iw_variance_error = np.power(iw_spectrum_variance_actual ,0.5)

    # 	print iw_spectrum_variance, "iw_spectrum_variance"
    # 	print iw_spectrum_variance_low, "iw_spectrum_variance_low"
    # 	print iw_spectrum_variance_high, "iw_spectrum_variance_high"
    ##	print iw_spectrum_variance_variance, "iw_spectrum_variance_variance"
    # 	print iw_variance_error, "iw_variance_error"
    #
    # 	print iw_spectrum_variance - iw_spectrum_variance_low
    # 	print iw_spectrum_variance - iw_spectrum_variance_high
    #
    # 	print np.log10(iw_spectrum_variance), "np.log10(iw_spectrum_variance)"
    # 	print np.log10(iw_spectrum_variance_low), "np.log10(iw_spectrum_variance_low)"
    # 	print np.log10(iw_spectrum_variance_high), "np.log10(iw_spectrum_variance_high)"
    # 	exit()

    ### Calculate reference GM variance
    interp = interpolate.interp1d(logkx_reference, logphi_reference)
    logphi_interp = interp(all_logkx)
    reference_spectrum_interp = np.power(10.0, logphi_interp)

    (
        gm_kx_integrand,
        gm_spectrum_integrand,
        gm_spectrum_variance,
    ) = estimate_spectral_variance(
        all_kx, reference_spectrum_interp, iw_kx_lower_limit_integration, iw_kx_max
    )

    ### Estimate dissipation rate and diffusivity for reference GM spectrum
    epsilon0, K0 = estimate_K0(E, N0_rads, b, jstar, N_GM_rads, f_GM, Rf, gm_option)
    logepsilon0 = np.log10(epsilon0)
    logK0 = np.log10(K0)

    ### Estimate dissipation rate and diffusivity for observed spectrum
    spectral_variance_ratio, epsilon, K = estimate_K(
        Rw,
        RwGM,
        N_rads,
        f,
        N_GM_rads,
        f_GM,
        epsilon0,
        K0,
        iw_spectrum_variance,
        gm_spectrum_variance,
    )
    logK = np.log10(K)
    logepsilon = np.log10(epsilon)

    ### Estimate error on logK
    spectral_ratio_std = estimate_standard_deviation_spectral_ratio(
        iw_spectrum_integrand, gm_spectrum_integrand
    )
    ### TODO Estimate this properly
    sigma_iw_logK = estimate_logK_sigma(
        spectral_ratio_std,
        iw_spectrum_integrand,
        gm_spectrum_integrand,
        iw_spectrum_variance,
        gm_spectrum_variance,
    )

    ### Estimate error on logK using standard deviation of average Welch spectrum
    spectral_variance_ratio_low, epsilon_low, K_low = estimate_K(
        Rw,
        RwGM,
        N_rads,
        f,
        N_GM_rads,
        f_GM,
        epsilon0,
        K0,
        iw_spectrum_variance_low,
        gm_spectrum_variance,
    )
    spectral_variance_ratio_high, epsilon_high, K_high = estimate_K(
        Rw,
        RwGM,
        N_rads,
        f,
        N_GM_rads,
        f_GM,
        epsilon0,
        K0,
        iw_spectrum_variance_high,
        gm_spectrum_variance,
    )

    logK_low = np.log10(K_low)
    logK_high = np.log10(K_high)

    logK_sigma_low = logK - logK_low
    logK_sigma_high = logK_high - logK

    ### Estimate error on logK using variance of integral
    ### TODO Add in errors in N
    # 	print np.log10(np.exp(1.))
    logK_variance = 2.0 * 0.4343 / iw_spectrum_variance * iw_variance_error
    logK_error = np.power(logK_variance, 0.5)

    print(logK, logK_error, "logK, logK_error")
    # 	print logK_low, "logK_low"
    # 	print logK_high, "logK_high"
    # 	print logK_sigma_low, "logK_sigma_low"
    # 	print logK_sigma_high, "logK_sigma_high"
    # 	print sigma_iw_logK, "sigma_iw_log"
    #
    # 	print logK_variance, "logK_variance"
    # 	print logK_error, "logK_error"

    ### TODO Fix errors properly. Error on logepsilon needs to include error in N as well.
    sigma_logepsilon_low = sigma_iw_logK
    sigma_logepsilon_high = sigma_iw_logK

    sigma_logK_low = sigma_iw_logK
    sigma_logK_high = sigma_iw_logK

    return (
        logepsilon0,
        logK0,
        logepsilon,
        logK,
        logK_error,
        sigma_logepsilon_low,
        sigma_logepsilon_high,
        sigma_logK_low,
        sigma_logK_high,
    )
