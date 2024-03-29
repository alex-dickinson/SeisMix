import numpy as np
from scipy.signal import hilbert
import matplotlib.pyplot as plt


def contour_reflections(
    data_select, contour_level, cmp_grid_select, twt_s_grid_select, extent
):
    """
    Convert seismic image to cosine of the instantaneous phase and draw contours at a constant value.

    Parameters
    ----------
    data_select : numpy.ndarray
        2D array containing seismic amplitudes within selected portion of image
    contour_level : float
        Value of cosine of the instantaneous phase to contour. Must lie between 0 and 1
    cmp_grid_select : numpy.ndarray
        2D array listing common midpoint number at every point in selected portion of seismic image
    twt_s_grid_select : numpy.ndarray
        2D array listing two-way travel time (in seconds) at every point in selected portion of seismic image
    extent : list
        List containing points defining boundaries of subregion of image within which contours are drawn

    Returns
    ----------
    cosine_ipa : numpy.ndarray
        2D array containing cosine of the instantaneous phase angle within selected portion of image
    all_contours : list
        List containing matplotlib.path.Path objects that define contours within subregion of image specified by `extent`
    """

    # Convert to cosine of instantaneous phase angle
    cosine_ipa = np.cos(np.unwrap(np.angle(hilbert(data_select))))
    # Remove outer traces of cosine of instantaneous phase so that contours are closed
    cosine_ipa[:, 0] = 0
    cosine_ipa[:, -1] = 0

    # Contour cosine of instantaneous phase angle at specified level
    cs = plt.contour(
        cmp_grid_select,
        twt_s_grid_select,
        cosine_ipa,
        levels=[-1 * contour_level, contour_level],
        extent=extent,
        colors="yellow",
        linestyles="solid",
        linewidths=1,
    )
    plt.close()

    # Extract list of contour arrays from LineCollection object
    positive_contours = cs.collections[0].get_paths()
    negative_contours = cs.collections[1].get_paths()
    all_contours = positive_contours + negative_contours

    return (cosine_ipa, all_contours)


def contour_to_midpoints(
    contour, min_length, cmp_gap, contours_twt_gap, midpoints_twt_gap
):
    if np.size(contour[:, 0]) < min_length:
        return ()

    # Round cmp values to nearest 0.001 of a cmp to avoid problems caused by precision
    contour[:, 0] = np.round(contour[:, 0], 3)

    ### Remove all rows with non-integer CMP.
    contour_new = np.zeros(contour.shape)
    for i in range(np.size(contour[:, 0])):
        if (contour[i, 0]).is_integer() == True:
            contour_new[i, :] = contour[i, :]
    contour = contour_new[~np.all(contour_new == 0, axis=1)]

    ### Sort contour based on CMP and then TWT
    contour = np.array(sorted(contour, key=lambda x: (x[0], x[1])))

    ### Find maximum number of TWTs which contour crosses at one CMP.
    max_number = 1
    k = 1
    for i in range(contour.shape[0] - 1):
        if contour[i, 0] == contour[i + 1, 0]:
            k = k + 1
        elif contour[i, 0] != contour[i + 1, 0]:
            if k >= max_number:
                max_number = k
            k = 1
    if max_number <= 1:
        return ()

    ### Arrange data into 2d array with each row corresponding to one CMP. $1 records the CMP, other columns list the TWTs at which the contour crosses this CMP.
    contour_array = np.zeros([contour.shape[0], max_number + 2])
    for i in range(contour.shape[0] - max_number):
        if contour[i, 0] == contour[i + 1, 0] and contour[i, 0] != contour[i - 1, 0]:
            row = i
            contour_array[row, 0] = contour[i, 0]
            contour_array[row, 1] = contour[i, 1]

            for j in range(1, max_number):
                if contour[i + j, 0] == contour[i, 0]:
                    contour_array[row, 1 + j] = contour[i + j, 1]

            ### Remove all rows with odd number of TWTs.
            number_twts = sum(x > 0 for x in contour_array[row, :]) - 1
            contour_array[row, -1] = number_twts
            if number_twts % 2 != 0:
                contour_array[row, :] = 0

    ### Remove all zero rows.
    contour_array = contour_array[~np.all(contour_array == 0, axis=1)]

    ### Split array whenever change in number of TWTs per CMP
    last_break = 0
    label = 0
    all_contours = {}
    for i in range(np.size(contour_array[:, 0]) - 1):
        if contour_array[i, -1] != contour_array[i + 1, -1]:
            previous = contour_array[last_break : i + 1, :]
            all_contours[str(label)] = previous
            label = label + 1
            last_break = i + 1
        if i == np.size(contour_array[:, 0]) - 2:
            if last_break != i + 1:
                previous = contour_array[last_break : i + 1, :]
                all_contours[str(label)] = previous

    ### Remove all values for which contour values at one CMP are separated by > contours_twt_gap
    all_midpoints = {}
    for key in all_contours:
        section = all_contours[str(key)]
        number = int(section[0, -1] / 2.0)
        for i in range(number):
            code = str(key) + "_" + str(i)
            divided_section = np.zeros([np.size(section[:, 0]), 3])
            for j in range(np.size(section[:, 0])):
                if abs(section[j, 1] - section[j, 2]) < contours_twt_gap:
                    divided_section[j, 0] = section[j, 0]
                    divided_section[j, 1] = section[j, 1 + 2 * i]
                    divided_section[j, 2] = section[j, 2 + 2 * i]
            divided_section = divided_section[~np.all(divided_section == 0, axis=1)]
            if np.size(divided_section[:, 0]) > 0:
                midpoints = np.zeros([np.size(divided_section[:, 0]), 2])
                midpoints[:, 0] = divided_section[:, 0]
                midpoints[:, 1] = (divided_section[:, 1] + divided_section[:, 2]) / 2.0
                all_midpoints[str(code)] = midpoints

    ### Link sections back together if appropriate
    if len(all_midpoints.keys()) == 0:
        return ()
    last_key = list(all_midpoints.keys())[-1]
    done = False

    def combine_midpoints(done):
        for i in all_midpoints:
            for j in all_midpoints:
                if i != j:
                    if (
                        abs(all_midpoints[str(i)][0, 0] - all_midpoints[str(j)][-1, 0])
                        <= cmp_gap
                        and abs(
                            all_midpoints[str(i)][0, 1] - all_midpoints[str(j)][-1, 1]
                        )
                        <= midpoints_twt_gap
                    ):
                        combined = np.append(
                            all_midpoints[str(j)], all_midpoints[str(i)], axis=0
                        )
                        all_midpoints[str(i)] = combined
                        all_midpoints.pop(str(j), None)
                        return done
                    elif (
                        abs(all_midpoints[str(j)][0, 0] - all_midpoints[str(i)][-1, 0])
                        <= cmp_gap
                        and abs(
                            all_midpoints[str(j)][0, 1] - all_midpoints[str(i)][-1, 1]
                        )
                        <= midpoints_twt_gap
                    ):
                        combined = np.append(
                            all_midpoints[str(i)], all_midpoints[str(j)], axis=0
                        )
                        all_midpoints[str(j)] = combined
                        all_midpoints.pop(str(i), None)
                        return done
                elif i == last_key and j == last_key:
                    done = True
                    return done

    while done == False:
        done = combine_midpoints(done)

    ### Break midpoints apart whenever adjacent CMP separation > cmp_gap or adjacent TWT separation > midpoints_twt_gap
    broken_midpoints = {}
    for padlock in all_midpoints:
        first_line = 0
        tag = 0
        midpoints = all_midpoints[str(padlock)]
        length_midpoints = np.size(midpoints[:, 0])
        for i in range(length_midpoints):
            if (
                i < (length_midpoints - 1)
                and abs(midpoints[i, 0] - midpoints[i + 1, 0]) > cmp_gap
            ):
                spray = str(padlock) + "_" + str(tag)
                broken_midpoints[str(spray)] = midpoints[first_line : i + 1, :]
                first_line = i + 1
                tag = tag + 1
            elif (
                i < (length_midpoints - 1)
                and abs(midpoints[i, 1] - midpoints[i + 1, 1]) > midpoints_twt_gap
            ):
                spray = str(padlock) + "_" + str(tag)
                broken_midpoints[str(spray)] = midpoints[first_line : i + 1, :]
                first_line = i + 1
                tag = tag + 1
            elif i == length_midpoints - 1:
                spray = str(padlock) + "_" + str(tag)
                broken_midpoints[str(spray)] = midpoints[first_line : i + 1, :]

    ### Only keep midpoints longer than min_length
    final_midpoints = {}
    for padlock in broken_midpoints:
        if np.size(broken_midpoints[str(padlock)][:, 0]) >= min_length:
            final_midpoints[str(padlock)] = broken_midpoints[str(padlock)][:, :]

    return final_midpoints


def compute_midpoints(
    all_contours, min_length, cmp_gap, contours_twt_gap, midpoints_twt_gap
):
    """
    Compute midpoints of contours computed by function `contour_reflections`. See figures in document ? for visualisation of strategy for computing midpoints.

    Parameters
    ----------
    all_contours : list
        List containing matplotlib.path.Path objects that define contours
    min_length : int
        Minimum length of computed midpoints to keep (in units of common midpoints)
    cmp_gap : int
        Maximum allowable horizontal gap (in units of common midpoints) between adjacent points along a particular computed midpoint. If tracked reflections are to be used for spectral analysis, cmp_gap should be set to 1
    contours_twt_gap : float
        Maximum allowable gap (in units of seconds) between the upper and lower surface of a single contour at a particular common midpoint. Midpoints will not be calculated at points at which the upper and lower surfaces of the contour are separated by a gap greater than contours_twt_gap
    midpoints_twt_gap : float
        Maximum allowable vertical gap (in units of seconds) between adjacent points along a particular computed midpoint

    Returns
    ----------
    all_contours_dict : dict
        Dictionary containing all input contours
    all_midpoints_dict : dict
        Dictionary containing all computed midpoints
    """

    all_contours_dict = {}
    all_midpoints_dict = {}

    contour_counter = 0

    for contour_path in all_contours:
        contour_counter = contour_counter + 1
        contour = contour_path.vertices  # Extract contour coordinates
        all_contours_dict[
            str(contour_counter)
        ] = contour  # Write contours to dictionary with contours numerically labelled
        if np.size(contour[:, 0]) >= min_length:
            contour_midpoints = contour_to_midpoints(
                contour, min_length, cmp_gap, contours_twt_gap, midpoints_twt_gap
            )
            for midpoints_label in contour_midpoints:
                all_midpoints_dict_label = (
                    str(contour_counter) + "_" + str(midpoints_label)
                )
                all_midpoints_dict[all_midpoints_dict_label] = contour_midpoints[
                    midpoints_label
                ]

    return (all_contours_dict, all_midpoints_dict)


def calculate_reflection_statistics(all_midpoints_dict, cmp_spacing_m):
    """
    Calculate statistics describing dataset of tracked reflections.

    Parameters
    ----------
    all_midpoints_dict : dict
        Dictionary containing midpoints of all tracked reflections
    cmp_spacing_m : float
        Distance between adjacent common midpoints in metres

    Returns
    ----------
    number_reflections : int
        number of tracked reflections
    reflection_length_array : numpy.ndarray
        array of length number_reflections, with each element equal to the length of an individual reflection (in units of common midpoints)
    overall_reflection_length_m : float
        cumulative length of all tracked reflections in metres
    mean_reflection_length_m : float
        mean length of all tracked reflections in metres
    std_reflection_length_m : float
        standard deviation of lengths of all tracked reflections in metres
    median_reflection_length_m : float
        median length of all tracked reflections in metres
    """

    number_reflections = len(all_midpoints_dict)
    reflection_length_array = np.zeros(number_reflections)

    index = 0
    for i in all_midpoints_dict:
        reflection_length_array[index] = np.size(all_midpoints_dict[i][:, 0])
        index += 1

    overall_reflection_length_m = cmp_spacing_m * np.sum(reflection_length_array)
    mean_reflection_length_m = np.mean(cmp_spacing_m * reflection_length_array)
    std_reflection_length_m = np.std(cmp_spacing_m * reflection_length_array)
    median_reflection_length_m = np.median(cmp_spacing_m * reflection_length_array)

    return (
        number_reflections,
        reflection_length_array,
        overall_reflection_length_m,
        mean_reflection_length_m,
        std_reflection_length_m,
        median_reflection_length_m,
    )
