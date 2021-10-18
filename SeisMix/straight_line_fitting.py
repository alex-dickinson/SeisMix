import numpy as np
from scipy import optimize


### Straight line fitting using least squares method. Returns both fitted gradient and fitted intercept, together with errors estimated from covariance matrix.
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
