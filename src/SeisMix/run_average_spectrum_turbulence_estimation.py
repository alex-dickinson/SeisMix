import numpy as np
from scipy import optimize

import gsw

import matplotlib.pyplot as plt


from python_mixing_functions import four_component_model_subrange_fitting
from python_mixing_functions import straight_line_fitting
from python_mixing_functions import estimate_turbulent_epsilon
from python_mixing_functions import estimate_internal_wave_epsilon
from python_mixing_functions import compute_gm76_horizontal_spectra
from python_mixing_functions import compute_gk91_horizontal_spectra
from python_mixing_functions import gm_model_fitting


###### Master wrapper script which calls other python modules to identify spectral subranges, fit straight lines to internal wave subrange, estimate values of epsilon from both turbulent and internal wave subranges, estimate values of K from values of epsilon ######
def estimate_mixing(
    avespec_array,
    N_cph,
    sigma_N_cph,
    latitude_degrees,
    gamma,
    CT,
    gm_option,
    N0_cph,
    N_GM_cph,
    E,
    b,
    jstar,
    RwGM,
    Rw,
    latitude0,
    Rf,
    Nkz,
    logkz_gm_min,
    logkz_gm_max,
    Nkx,
    logkx_gm_min,
    logkx_gm_max,
    NkH,
    gm_rolloff_option,
    min_turb_length,
    min_iw_length,
    averaging_option,
):

    ### Convert cph to rads
    N_rads = 2.0 * np.pi * N_cph / 3600.0
    sigma_N_rads = 2.0 * np.pi * sigma_N_cph / 3600.0
    N0_rads = 2.0 * np.pi * N0_cph / 3600.0
    N_GM_rads = 2.0 * np.pi * N_GM_cph / 3600.0

    ### Set up latitude-dependent variables
    f = gsw.f(latitude_degrees)
    latitude0 = float(latitude0)
    f0 = gsw.f(latitude0)
    f_GM = f0

    ### Load data
    all_kx = avespec_array[:, 0]
    all_phi = avespec_array[:, 1]
    all_sigma_phi = avespec_array[:, 2]
    all_logkx = avespec_array[:, 3]

    ### Specify whether to use log10 of spectra averaged in linear space or average spectrum averaged in log-space.
    if averaging_option == "log_mean_linear":
        all_logphi = avespec_array[:, 4]
        ### TODO Define logphi_std properly
        all_logphi_std = avespec_array[:, 5]
    elif averaging_option == "mean_log":
        all_logphi = avespec_array[:, 6]
        ### TODO Define logphi_std properly
        all_logphi_std = avespec_array[:, 7]

    ### Remove first 2 values to eliminate spectra roll-off - TODO make variable
    ### TODO Decide whether to take log10 of mean or mean of log10
    logkx = all_logkx[2::]
    logphi = all_logphi[2::]
    logphi_std = all_logphi_std[2::]
    kx = all_kx[2::]
    ### TODO phi and sigma_phi are currently average in linear space. If using logspace average, do need to use 10^(all_logphi)? Used in IW calculation.
    phi = all_phi[2::]
    sigma_phi = all_sigma_phi[2::]

    ### Identify spectral subranges by fitting four-component model to empirical spectrum
    par = [-2, -3, 0.5, -1.0]  # Initial guess of model parameters
    fitted_pars = four_component_model_subrange_fitting.subrange_fitting(
        logkx, logphi, par
    )

    print(fitted_pars)

    ###### Isolate internal wave, turbulent and noise subranges.
    ### Four parameter model.
    # print fitted_pars
    iw_turb_x = fitted_pars[0]
    iw_turb_arg = np.argmin(logkx < iw_turb_x)
    turb_noise_x = iw_turb_x + fitted_pars[2]
    iw_logkx = logkx[0:iw_turb_arg]
    iw_logphi = logphi[0:iw_turb_arg]
    iw_logphi_std = logphi_std[0:iw_turb_arg]

    if turb_noise_x < np.max(logkx):
        turb_noise_arg = np.argmin(logkx < turb_noise_x)
        turb_logkx = logkx[iw_turb_arg:turb_noise_arg]
        turb_logphi = logphi[iw_turb_arg:turb_noise_arg]
        turb_logphi_std = logphi_std[iw_turb_arg:turb_noise_arg]
        noise_logkx = logkx[turb_noise_arg::]
        noise_logphi = logphi[turb_noise_arg::]
        noise_logphi_std = logphi_std[turb_noise_arg::]
    else:
        turb_logkx = logkx[iw_turb_arg::]
        turb_logphi = logphi[iw_turb_arg::]
        turb_logphi_std = logphi_std[iw_turb_arg::]
        noise_logkx = []
        noise_logphi = []
        noise_logphi_std = []

    plt.plot(logkx, logphi)
    plt.plot(turb_logkx, turb_logphi)
    #     plt.plot(turb_logkx, c2_turb + 1./3.*turb_logkx)
    plt.show()

    ### Analyse turbulent subrange if it contains min_turb_length points or more.
    if np.size(turb_logkx) >= min_turb_length:

        ### Fit straight line through turbulent subrange using least squares method.
        ### x0 is guess of initial values for gradient and intercept.
        ### TODO x0 can be defined in module
        x0 = np.array([0.0, 1.0])

        (
            m1_turb_fit,
            sigma_m1_turb_fit,
            c1_turb_fit,
            sigma_c1_turb_fit,
        ) = straight_line_fitting.fit_straight_line(
            turb_logkx, turb_logphi, turb_logphi_std, x0
        )

        #         print(m1_turb_fit, "m1_turb_fit")
        #         print(sigma_m1_turb_fit, "sigma_m1_turb_fit")

        plt.plot(logkx, logphi)
        plt.plot(turb_logkx, turb_logphi)

        ### TODO Estimate logK from turbulent subrange if fitted straight line has gradient of 1/3 within error
        # 	if m1_turb_fit - sigma_m1_turb_fit <= 0.33 and m1_turb_fit + sigma_m1_turb_fit >= 0.33 and sigma_m1_turb_fit <= 0.33:
        if m1_turb_fit - 100.0 <= 1.0 / 3.0 and m1_turb_fit + 100.0 >= 1.0 / 3.0:
            #             print "turb"

            ### Estimate epsilon from turbulent subrange by fitting straight line of gradient 1/3 ### TODO Check straight line fitting
            (
                c2_turb,
                sigma_c2_turb,
                logepsilon_turb,
                sigma_logepsilon_turb,
                logK_turb,
                sigma_logK_turb,
            ) = estimate_turbulent_epsilon.estimate_epsilon_from_turbulent_subrange(
                turb_logkx,
                turb_logphi,
                turb_logphi_std,
                N_rads,
                sigma_N_rads,
                gamma,
                CT,
            )

            # 		plt.plot(logkx, logphi)
            # 		plt.plot(turb_logkx, turb_logphi)
            plt.plot(turb_logkx, c2_turb + 1.0 / 3.0 * turb_logkx)
            # 		plt.show()
            if m1_turb_fit - 0.3 <= 1.0 / 3.0 and m1_turb_fit + 0.3 >= 1.0 / 3.0:
                good_turb_marker = "good"
            else:
                good_turb_marker = "bad"

        # 		print logepsilon_turb, sigma_logepsilon_turb, "logepsilon_turb", "sigma_logepsilon_turb"
        #             print logK_turb, sigma_logK_turb, "logK_turb", "sigma_logK_turb"
        else:
            c2_turb = np.nan
            sigma_c2_turb = np.nan
            logepsilon_turb = np.nan
            sigma_logepsilon_turb = np.nan
            logK_turb = np.nan
            sigma_logK_turb = np.nan

    # 	plt.show()
    else:
        m1_turb_fit = np.nan
        sigma_m1_turb_fit = np.nan
        c1_turb_fit = np.nan
        sigma_c1_turb_fit = np.nan
        c2_turb = np.nan
        sigma_c2_turb = np.nan
        logepsilon_turb = np.nan
        sigma_logepsilon_turb = np.nan
        logK_turb = np.nan
        sigma_logK_turb = np.nan
        good_turb_marker = np.nan

    ### Analyse internal wave subrange
    if np.size(iw_logkx) >= min_iw_length:
        # 	print "iw"
        # 	print np.size(iw_logkx)

        ### Fit spectral models to internal wave spectral subrange
        ### Fit straight line through internal wave subrange using least squares method.
        ### x0 is guess of initial values for gradient and intercept.
        ### ### TODO x0 can be defined in module
        x0 = np.array([0.0, 1.0])
        (
            m1_iw_fit,
            sigma_m1_iw_fit,
            c1_iw_fit,
            sigma_c1_iw_fit,
        ) = straight_line_fitting.fit_straight_line(
            iw_logkx, iw_logphi, iw_logphi_std, x0
        )

        # 	print m1_iw_fit, "m1_iw_fit"
        # 	print sigma_m1_iw_fit, "sigma_m1_iw_fit"
        # 	print c1_iw_fit, "c1_iw_fit"
        # 	print sigma_c1_iw_fit, "sigma_c1_iw_fit"

        ### Fit GM spectral form to internal wave subrange. gm_fit_initial is initial guess of values for A_gm_fit, s_gm_fit, t_gm_fit, kstar_gm_fit based on GM75 model parameters.
        gm_fit_initial = np.array([100.0, 1.0, 2.5, 0.0001])
        (
            A_gm_fit,
            s_gm_fit,
            t_gm_fit,
            kstar_gm_fit,
        ) = gm_model_fitting.fit_gm75_kx_parametrisation_leastsq(
            iw_logkx, iw_logphi, iw_logphi_std, gm_fit_initial
        )

        ### Estimate epsilon from internal wave subrange using fine-scale shear-strain parametrisation with fixed shear-strain ratio and GM76 as reference spectrum.
        ### TODO If using constant buoyancy frequency for all pycnocline average spectra then this GM spectrum can be computed only once, which will be computationally more efficient.
        ### Compute reference GM spectrum at local buoyancy frequency.
        if gm_option == "gm76":
            (
                logkx_reference,
                power_spectrum_reference,
                logphi_reference,
            ) = compute_gm76_horizontal_spectra.run_kx_from_kxkz(
                N_rads,
                N0_rads,
                f,
                Nkz,
                logkz_gm_min,
                logkz_gm_max,
                Nkx,
                logkx_gm_min,
                logkx_gm_max,
                NkH,
                E,
                b,
                jstar,
                gm_rolloff_option,
            )
        elif gm_option == "gk91":
            (
                logkx_reference,
                power_spectrum_reference,
                logphi_reference,
            ) = compute_gk91_horizontal_spectra.run_kx_from_kxkz(
                N_rads,
                N0_rads,
                f,
                Nkz,
                logkz_gm_min,
                logkz_gm_max,
                Nkx,
                logkx_gm_min,
                logkx_gm_max,
                NkH,
                E,
                b,
                jstar,
                gm_rolloff_option,
            )
        else:
            print("gm_option not specified - exiting")
            exit()

        ### Make fitted spectrum for output
        logphi_fitted_gm = gm_model_fitting.make_fitted_spectrum(
            logkx_reference, A_gm_fit, s_gm_fit, t_gm_fit, kstar_gm_fit
        )

        ### TODO Currently integrating over entire internal wave subrange for fine-scale shear-strain parametrisation. May want to change.
        iw_kx_lower_limit_integration = np.min(iw_logkx)
        iw_kx_upper_limit_integration = np.max(iw_logkx)
        ### TODO Clear up this script and check through
        ### TODO Estimate error on logK_iw correctly
        (
            logepsilon0,
            logK0,
            logepsilon_iw,
            logK_iw,
            logK_error_iw,
            sigma_logepsilon_iw_low,
            sigma_logepsilon_iw_high,
            sigma_logK_iw_low,
            sigma_logK_iw_high,
        ) = estimate_internal_wave_epsilon.estimate_epsilon_from_internal_wave_subrange(
            logkx,
            logphi,
            sigma_phi,
            logphi_std,
            iw_logkx,
            iw_logphi,
            iw_logphi_std,
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
        )

    #         print(logepsilon_iw, sigma_logepsilon_iw_low, "logepsilon_iw", "sigma_logepsilon_iw")
    #         print(logK_iw, sigma_logK_iw_low, "logK_iw", "sigma_logK_iw")
    #         print(logK_iw, logK_error_iw, "logK_iw", "sigma_logK_iw")

    else:
        m1_iw_fit = np.nan
        sigma_m1_iw_fit = np.nan
        c1_iw_fit = np.nan
        sigma_c1_iw_fit = np.nan
        A_gm_fit = np.nan
        t_gm_fit = np.nan
        s_gm_fit = np.nan
        kstar_gm_fit = np.nan
        iw_kx_lower_limit_integration = np.nan
        iw_kx_upper_limit_integration = np.nan
        logepsilon0 = np.nan
        logK0 = np.nan
        logepsilon_iw = np.nan
        logK_iw = np.nan
        logK_error_iw = np.nan
        sigma_logepsilon_iw_low = np.nan
        sigma_logepsilon_iw_high = np.nan
        sigma_logK_iw_low = np.nan
        sigma_logK_iw_high = np.nan
        logkx_reference = np.zeros(2)
        logkx_reference[:] = np.nan
        logphi_reference = np.zeros(2)
        logphi_reference[:] = np.nan
        logphi_fitted_gm = np.zeros(2)
        logphi_fitted_gm[:] = np.nan

    min_logkx = np.min(logkx)


#     Append all this guff to dictionaries

#     ### Save all data
#     ### TODO Remove outfile
#     if os.path.isfile(outfile_name + '.dat') == True:
#         os.remove(outfile_name + '.dat')
#     outfile = open(str(outfile_name) + '.dat', 'a')
#     outfile.write("> -Z ")
#     outfile.write("%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %s\n" % (N_rads, sigma_N_rads, min_logkx, fitted_pars[0], fitted_pars[1], fitted_pars[2], fitted_pars[3], min_turb_length, m1_turb_fit, sigma_m1_turb_fit, c1_turb_fit, sigma_c1_turb_fit, c2_turb, sigma_c2_turb, logepsilon_turb, sigma_logepsilon_turb, logK_turb, sigma_logK_turb, min_iw_length, m1_iw_fit, sigma_m1_iw_fit, c1_iw_fit, sigma_c1_iw_fit, A_gm_fit, s_gm_fit, t_gm_fit, kstar_gm_fit, iw_kx_lower_limit_integration, iw_kx_upper_limit_integration, logepsilon0, logK0, logepsilon_iw, logK_iw, logK_error_iw, sigma_logepsilon_iw_low, sigma_logepsilon_iw_high, sigma_logK_iw_low, sigma_logK_iw_high, good_turb_marker))
#     for i in range(np.size(all_logkx[:])):
#         outfile.write("%f %f %f\n" % (all_logkx[i], all_logphi[i], all_logphi_std[i]))
#     outfile.close()


#     ### TODO Write outfile containing all arguments and their positions
#     comments_outfile = open(str(comments_name) + '.dat', 'w')
#     comments_outfile.write("### All spectra headers\n")
#     comments_outfile.write(">, - Z, N_rads, sigma_N_rads, min_logkx, fitted_pars[0], fitted_pars[1], fitted_pars[2], fitted_pars[3], min_turb_length, m1_turb_fit, sigma_m1_turb_fit, c1_turb_fit, sigma_c1_turb_fit, c2_turb, sigma_c2_turb, logepsilon_turb, sigma_logepsilon_turb, logK_turb, sigma_logK_turb, min_iw_length, m1_iw_fit, sigma_m1_iw_fit, c1_iw_fit, sigma_c1_iw_fit, A_gm_fit, s_gm_fit, t_gm_fit, kstar_gm_fit, iw_kx_lower_limit_integration, iw_kx_upper_limit_integration, logepsilon0, logK0, logepsilon_iw, logK_iw, logK_error_iw, sigma_logepsilon_iw_low, sigma_logepsilon_iw_high, sigma_logK_iw_low, sigma_logK_iw_high, good_turb_marker")
#     comments_outfile.close()


#     ### Write outfile containing reference GM spectrum
#     ### TODO Remove outfile
#     if os.path.isfile(gm_reference_spectra_file + '.dat') == True:
#         os.remove(gm_reference_spectra_file + '.dat')
#     gm_reference_spectra_outfile = open(str(gm_reference_spectra_file) + '.dat', 'a')
#     gm_reference_spectra_outfile.write("> Reference GM spectrum\n")
#     gm_reference_spectra_outfile.write(">")
#     for i in range(np.size(logkx_reference[:])):
#         gm_reference_spectra_outfile.write("%f %f\n" % (logkx_reference[i], logphi_reference[i]))
#     gm_reference_spectra_outfile.close()


#     ### Write outfile containing fitted GM spectrum
#     ### TODO Remove outfile
#     if os.path.isfile(gm_fitted_spectra_file + '.dat') == True:
#         os.remove(gm_fitted_spectra_file + '.dat')
#     gm_fitted_spectra_outfile = open(str(gm_fitted_spectra_file) + '.dat', 'a')
#     gm_fitted_spectra_outfile.write("> Fitted GM spectrum\n")
#     gm_fitted_spectra_outfile.write("> ")
#     gm_fitted_spectra_outfile.write("%f %f %f %f\n" % (A_gm_fit, t_gm_fit, s_gm_fit, kstar_gm_fit))
#     for i in range(np.size(logkx_reference[:])):
#         gm_fitted_spectra_outfile.write("%f %f\n" % (logkx_reference[i], logphi_fitted_gm[i]))
#     gm_fitted_spectra_outfile.close()
