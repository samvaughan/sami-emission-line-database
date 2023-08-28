import sqlite3
from astropy.table import Table
from astropy.io import fits

hdu = fits.open(
    "/Users/samvaughan/Science/SAMI/data/Tables/jvds_stelkin_cat_v012_mge_seecorr_kh20_v20220604_private.fits"
)
df = Table(hdu[1].data).to_pandas()

con = sqlite3.connect("")