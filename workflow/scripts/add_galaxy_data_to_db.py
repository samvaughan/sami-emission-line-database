import pandas as pd
from astropy.io import fits
import numpy as np
from astropy.wcs import WCS
from pathlib import Path
import astropy.wcs.utils as wcs_utils
import sqlite3
import utils

smk = snakemake  # noqa
conn = sqlite3.connect(smk.input.database)

# The data directory for each galaxy
data_dir = smk.params.data_dir
catid = smk.wildcards.CATID

hdu = fits.open(f"{data_dir}/{catid}_A_Halpha_default_1-comp.fits")

# Get the WCS
wcs = WCS(hdu[0].header)

# Index for each spaxel
x, y = np.indices((50, 50))
coords = wcs_utils.pixel_to_skycoord(x, y, wcs)

RA = coords.ra.ravel()
DEC = coords.dec.ravel()

halpha, halpha_error = utils.load_gas_quantity(catid, "Halpha", data_dir)
hbeta, hbeta_error = utils.load_gas_quantity(catid, "Hbeta", data_dir)
n_II, n_II_error = utils.load_gas_quantity(catid, "NII6583", data_dir)
o_I, o_I_error = utils.load_gas_quantity(catid, "OI6300", data_dir)
o_II, o_II_error = utils.load_gas_quantity(catid, "OII3728", data_dir)
o_III, o_III_error = utils.load_gas_quantity(catid, "OIII5007", data_dir)
s_II_a, s_II_a_error = utils.load_gas_quantity(catid, "SII6716", data_dir)
s_II_b, s_II_b_error = utils.load_gas_quantity(catid, "SII6731", data_dir)
extinct_corr, extinct_corr_error = utils.load_gas_quantity(
    catid, "extinct-corr", data_dir
)
gas_vdisp, gas_vdisp_error = utils.load_gas_quantity(catid, "gas-vdisp", data_dir)
gas_vel, gas_vel_error = utils.load_gas_quantity(catid, "gas-velocity", data_dir)
sfr, sfr_error = utils.load_gas_quantity(catid, "sfr", data_dir)
sfr_density, sfr_density_error = utils.load_gas_quantity(catid, "sfr-dens", data_dir)
(
    stellar_vdisp,
    stellar_vdisp_error,
    stellar_flux,
    stellar_flux_error,
    stellar_SN,
) = utils.load_stellar_quantity(catid, "stellar-velocity-dispersion", data_dir)
(
    stellar_vel,
    stellar_vel_error,
    _,
    _,
    _,
) = utils.load_stellar_quantity(catid, "stellar-velocity", data_dir)

(
    h3,
    h3_error,
    _,
    _,
    _,
) = utils.load_stellar_quantity(catid, "stellar-velocity-h3", data_dir)

(
    h4,
    h4_error,
    _,
    _,
    _,
) = utils.load_stellar_quantity(catid, "stellar-velocity-h4", data_dir)

table = pd.DataFrame(
    data=dict(
        CATID=np.full(RA.shape, catid),
        RA=RA,
        DEC=DEC,
        cube_x=x.ravel(),
        cube_y=y.ravel(),
        stellar_flux=stellar_flux.ravel(),
        e_stellar_flux=stellar_flux_error.ravel(),
        stellar_SN=stellar_SN.ravel(),
        v_stars=stellar_vel.ravel(),
        e_v_stars=stellar_vel_error.ravel(),
        sigma_stars=stellar_vdisp.ravel(),
        e_sigma_stars=stellar_vdisp_error.ravel(),
        h3=h3.ravel(),
        e_h3=h3_error.ravel(),
        h4=h4.ravel(),
        e_h4=h4_error.ravel(),
        v_gas=gas_vel.ravel(),
        e_v_gas=gas_vel_error.ravel(),
        sigma_gas=gas_vdisp.ravel(),
        e_sigma_gas=gas_vdisp_error.ravel(),
        sfr=sfr[0].ravel(),
        e_sfr=sfr_error[0].ravel(),
        sfr_density=sfr_density[0].ravel(),
        e_sfr_density=sfr_density_error[0].ravel(),
        halpha_flux=halpha[0].ravel(),
        e_halpha_flux=halpha_error[0].ravel(),
        hbeta_flux=hbeta.ravel(),
        e_hbeta_flux=hbeta_error.ravel(),
        NII6583_flux=n_II.ravel(),
        e_NII6583_flux=n_II_error.ravel(),
        OI6300_flux=o_I.ravel(),
        e_OI6300_flux=o_I_error.ravel(),
        OII3728_flux=o_II.ravel(),
        e_OII3728_flux=o_II.ravel(),
        OIII5007_flux=o_III.ravel(),
        e_OIII5007_flux=o_III_error.ravel(),
        SII6716_flux=s_II_a.ravel(),
        e_SII6716_flux=s_II_a_error.ravel(),
        SII6731_flux=s_II_b.ravel(),
        e_SII6731_flux=s_II_b_error.ravel(),
        extinct_corr=extinct_corr.ravel(),
        e_extinct_corr=extinct_corr_error.ravel(),
    )
)

table["log10_NII_Halpha"] = np.log10(table.NII6583_flux / table.halpha_flux)
table["log10_OI_Halpha"] = np.log10(table.OI6300_flux / table.halpha_flux)
table["log10_OIII_Hbeta"] = np.log10(table.OIII5007_flux / table.hbeta_flux)
table["log10_SII_Halpha"] = np.log10(table.SII6731_flux / table.halpha_flux)
table["BPT_class"] = utils.classify_spaxel_BPT_vectorised(
    table["log10_OIII_Hbeta"], table["log10_SII_Halpha"]
)

table.to_sql("sami", conn, index=False, if_exists="append")
Path(smk.output.flag_file).touch()
