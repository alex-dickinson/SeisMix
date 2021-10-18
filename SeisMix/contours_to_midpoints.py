import numpy as np

# import pickle
# import json

from scipy import signal
from scipy.optimize import fmin

import matplotlib.pyplot as plt
# import pylab

from matplotlib.font_manager import FontProperties


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
