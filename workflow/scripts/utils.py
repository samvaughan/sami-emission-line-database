import numpy as np
import scipy.stats as stats
from astropy.io import fits


def Kewley_2006_SF(x_ratio_name, ratio):
    """Classification scheme from Kewley et al. 2006 equations 1-3. SF galaxies lie below each of these lines

    Args:
        x_ratio_name (str): The choice of x-value for the Kewley 2001 equation. Must be one of 'NII', 'SII' or 'OI'
        ratio (array-like): The x values in the equation

    Returns:
        array-like: The values of log(OIII/Hbeta) on the Kewley 2006 lines
    """
    assert x_ratio_name in [
        "NII",
        "SII",
        "OI",
    ], "The name of the x-axis must be one of NII, SII or OI"
    if x_ratio_name == "SII":
        return 0.72 / (ratio - 0.32) + 1.30
    elif x_ratio_name == "NII":
        return 0.61 / (ratio - 0.05) + 1.30
    elif x_ratio_name == "OI":
        return 0.73 / (ratio + 0.59) + 1.33
    else:
        raise NameError(f"x_ratio_name of {x_ratio_name} not understood")


def Kewley_2006_Composite(x_ratio_name, ratio):
    """Classification scheme from Kewley et al. 2006 equations 4 and 5. Composite galaixes lie between these two lines.

    Args:
        ratio (array-like): The x values in the equation, in this case [NII/Halpha]

    Returns:
        array-like: The values of log(OIII/Hbeta) on the Kewley 2006 lines
    """

    assert (
        x_ratio_name == "NII"
    ), "The name of the x-axis must be NII for this classification"

    left_limit = 0.61 / (ratio - 0.05) + 1.3
    right_limit = 0.61 / (ratio - 0.47) + 1.19

    return left_limit, right_limit


def Kewley_2006_Seyfert_LINER(x_ratio_name, ratio):
    assert x_ratio_name in [
        "SII",
        "OI",
    ], "The name of the x-axis must be one of SII or OI for this classification"
    if x_ratio_name == "SII":
        return 1.89 * ratio + 0.76
    elif x_ratio_name == "OI":
        return 1.18 * ratio + 1.30


def classify_spaxel_BPT(log10_OIII_Hbeta_value, log10_SII_Halpha_value):
    """
    Classify a **single** value of emission line ratios from a spectrum.
    Follow the classification scheme of Kewley et al. 2006 (https://ui.adsabs.harvard.edu/abs/2006MNRAS.372..961K/abstract) but _only_ using their SII/Halpha based categories
    This follows the work of Belfiore et al. 2016 (https://ui.adsabs.harvard.edu/abs/2016MNRAS.461.3111B/abstract) and makes things **much** simpler.

    The classification scheme is as follows:

    * Take the log_10(OIII/Hbeta) value of a spaxel
    * Check if it is to the LEFT of the Kewley 2006 SF line in the OIII/HBeta vs SII/Halpha diagram.
        --> If yes, the spaxel is STAR FORMING
    * If it is to the right of the Kewley 2006 line, see whether it is ABOVE the Kewley et al. 2006 Seyfert/Liner line
        --> If it is ABOVE, the spaxel is a SEYFERT
        --> If it is BELOW, the spaxel is a LINER

    We return 0 for SF spaxels, 1 for LINERS and 2 for SEYFERTs. Missing data gets np.NaN

    Args:
        log10_OIII_Hbeta_value (float): The value of log10(OIII/Halpha) flux in a spaxel
        log10_OIII_Hbeta_value (float): The value of log10(OIII/Halpha) flux in a spaxel

    Returns:
        int: The classification of the spaxel- 0, 1, 2 or NaN
    """

    if not np.all(np.isfinite([log10_OIII_Hbeta_value, log10_SII_Halpha_value])):
        return np.nan
    elif log10_OIII_Hbeta_value < Kewley_2006_SF("SII", log10_SII_Halpha_value):
        # The Kewley_2006_SF diving line is only valid for log10_SII_Halpha_value ratios less than 0.32 (check the equation)
        if log10_SII_Halpha_value < 0.32:
            return 0
        else:
            pass
    elif log10_OIII_Hbeta_value < Kewley_2006_Seyfert_LINER(
        "SII", log10_SII_Halpha_value
    ):
        # Spaxel is LINER
        return 1
    else:
        # Spaxel is SEYFERT
        return 2


def classify_spaxel_BPT_vectorised(log10_OIII_Hbeta_values, log10_SII_Halpha_values):
    """
    Classify an array of emission line ratios from a spectrum.
    Follow the classification scheme of Kewley et al. 2006 (https://ui.adsabs.harvard.edu/abs/2006MNRAS.372..961K/abstract) but _only_ using their SII/Halpha based categories
    This follows the work of Belfiore et al. 2016 (https://ui.adsabs.harvard.edu/abs/2016MNRAS.461.3111B/abstract) and makes things **much** simpler.

    The classification scheme is as follows:

    * Take the log_10(OIII/Hbeta) value of a spaxel
    * Check if it is to the LEFT of the Kewley 2006 SF line in the OIII/HBeta vs SII/Halpha diagram.
        --> If yes, the spaxel is STAR FORMING
    * If it is to the right of the Kewley 2006 line, see whether it is ABOVE the Kewley et al. 2006 Seyfert/Liner line
        --> If it is ABOVE, the spaxel is a SEYFERT
        --> If it is BELOW, the spaxel is a LINER

    We return 0 for SF spaxels, 1 for LINERS and 2 for SEYFERTs. Missing data gets np.NaN

    Args:
        log10_OIII_Hbeta_values (array-like): The value of log10(OIII/Halpha) flux in a spaxel
        log10_OIII_Hbeta_values (array-like): The value of log10(OIII/Halpha) flux in a spaxel

    Returns:
        array-like: The classification of the spaxel- 0, 1, 2 or NaN
    """

    BPT_classifications = np.full(len(log10_SII_Halpha_values), -99.0)
    Kewley_SF_values = Kewley_2006_SF("SII", log10_SII_Halpha_values)
    Kewley_LINER_values = Kewley_2006_Seyfert_LINER("SII", log10_SII_Halpha_values)

    # Do the classifications

    # Bad/missing data
    BPT_classifications[
        ~(np.isfinite(log10_OIII_Hbeta_values) & (np.isfinite(log10_SII_Halpha_values)))
    ] = np.nan
    # Star-forming
    BPT_classifications[
        (log10_OIII_Hbeta_values < Kewley_SF_values) & (log10_SII_Halpha_values < 0.32)
    ] = 0
    # Liner
    BPT_classifications[
        (
            (log10_OIII_Hbeta_values > Kewley_SF_values)
            | (log10_SII_Halpha_values > 0.32)
        )
        & (log10_OIII_Hbeta_values < Kewley_LINER_values)
    ] = 1
    # Seyfert
    BPT_classifications[
        (
            (log10_OIII_Hbeta_values > Kewley_SF_values)
            | (log10_SII_Halpha_values > 0.32)
        )
        & (log10_OIII_Hbeta_values > Kewley_LINER_values)
    ] = 2

    return BPT_classifications


# Errors on a ratio
# z = x/y
# If x and y are normally distributed, with means mu_x, mu_y and sigma_x, sigma_y
# See https://en.wikipedia.org/wiki/Ratio_distribution#Normal_ratio_distributions
def a(z, sigma_x, sigma_y):
    a_z = np.sqrt(z**2 / sigma_x**2 + 1.0 / sigma_y**2)
    return a_z


def b(z, mu_x, sigma_x, mu_y, sigma_y):
    b_z = mu_x * z / sigma_x**2 + mu_y / sigma_y**2
    return b_z


def c(mu_x, sigma_x, mu_y, sigma_y):
    return mu_x**2 / sigma_x**2 + mu_y**2 / sigma_y**2


def d(z, mu_x, sigma_x, mu_y, sigma_y):
    a_z = a(z, sigma_x, sigma_y)
    b_z = b(z, mu_x, sigma_x, mu_y, sigma_y)
    c_z = c(mu_x, sigma_x, mu_y, sigma_y)

    return np.exp((b_z**2 - c_z * a_z**2) / (2 * a_z**2))


def p(z, mu_x, sigma_x, mu_y, sigma_y):
    a_z = a(z, sigma_x, sigma_y)
    b_z = b(z, mu_x, sigma_x, mu_y, sigma_y)
    c_z = c(mu_x, sigma_x, mu_y, sigma_y)
    d_z = d(z, mu_x, sigma_x, mu_y, sigma_y)
    p_z_part_one = (
        ((b_z * d_z) / a_z**3)
        * 1.0
        / (np.sqrt(2 * np.pi) * sigma_x * sigma_y)
        * (stats.norm.cdf(b_z / a_z) - stats.norm.cdf(-b_z / a_z))
    )
    p_z_part_two = 1.0 / (a_z**2 * np.pi * sigma_x * sigma_y) * np.exp(-c_z / 2)

    return p_z_part_one + p_z_part_two


def load_gas_quantity(catid, quantity, data_dir, cube="default", type="1-comp"):
    try:
        hdu = fits.open(f"{data_dir}/{catid}_A_{quantity}_{cube}_{type}.fits")
        img = hdu[0].data
        error = hdu[1].data
    except FileNotFoundError:
        shape = (50, 50)
        if quantity in ["Halpha", "sfr", "sfr-dens"]:
            shape = (2, 50, 50)
        img = np.full(shape, np.nan)
        error = np.full(shape, np.nan)

    return img, error


def load_stellar_quantity(
    catid, quantity, data_dir, cube="default", type="four-moment"
):
    hdu = fits.open(f"{data_dir}/{catid}_A_{quantity}_{cube}_{type}.fits")
    img = hdu[0].data
    error = hdu[1].data
    flux = hdu[2].data
    flux_error = hdu[3].data
    SN = hdu[4].data

    return img, error, flux, flux_error, SN
