import numpy as np

# ----------------------------------------------------------------

# Round values to number of significant figures determined by uncertainty
def round_to_sf(x, sigma_x):
    if x == 0 or sigma_x == 0:
        x_round = x
    else:
        sigma_x_power = -int(np.floor(np.log10(abs(sigma_x))))
        sigma_x_first_digit = round(sigma_x, sigma_x_power) * np.power(
            10, float(sigma_x_power)
        )

        if sigma_x_first_digit < 3:
            n = 2
        else:
            n = 1

        sigma_x_round = round(sigma_x, sigma_x_power + (n - 1))
        x_round = round(x, sigma_x_power + (n - 1))

    return (x_round, sigma_x_round)


# ----------------------------------------------------------------
