# SAMI emission lines database for the 7th COIN residence programme

This repo contains code which collates various SAMI data products into a single SQLite database, which will be used for the COIN residence program. 

This repo _doesn't_ contain the individual SAMI dataproducts or the final database, since they're each a few Gb in size. The data products can be downloaded from [the bulk download tool at Data Central](https://datacentral.org.au/services/download/) and the database itself can be found on the COIN server (location tbd).If you've downloaded the SAMI data, you can run all the scripts and build the database using [snakemake](https://snakemake.readthedocs.io/en/stable/) but I don't imagine that this will be necessary!

## Scripts

The scripts to build the database are in `workflow/scripts`. They do the following:

* `make_database.py` sets up a new SQLite database at `results/SAMI_spaxel_data.db` with the appropriate columns.
* `add_galaxy_data_to_db.py` reads in the various gas and star measurements for a single SAMI galaxy, adds some extra columns (i.e. calculates some important emission line ratios and adds the BPT classification for each spaxel) and then adds this data to the master database.
* `utils.py` contains utility functions which are used by `add_galaxy_data_to_db.py`. In particular, it contains the functions which do the BPT classifications based on different emission line ratios. It also contains functions which read in the raw data, as well as a few things I've been playing around with (and aren't used yet). 

## Full database schema

The database has one table called `sami`, which contains measured values for every spectrum for every galaxy in the SAMI survey. Each row in the table corresponds to one spectrum. 

The table `sami` has the following columns:

* __id__ (INTEGER PRIMARY KEY): The row number in the table.
* __CATID__ (INTEGER): The ID number of the SAMI galaxy. There are 3068 unique CATID values.
* __RA__ (FLOAT): The Right Ascension of the spectrum (J2000).
* __DEC__ (FLOAT): The Declination of the spectrum (J2000).
* __cube_x__ (INTEGER): The x index of the spectrum within the datacube corresponding to this galaxy. Will be a number between 0 and 50, with the centre of the galaxy at (25, 25).
* __cube_y__ (INTEGER): The y index of the spectrum within the datacube corresponding to this galaxy. Will be a number between 0 and 50, with the centre of the galaxy at (25, 25).
* __stellar_flux__ (FLOAT): The stellar continuum flux of the spectrum, found within the $V$-band filter (I _think_)
* __e_stellar_flux__ (FLOAT): The uncertainty on the stellar continuum flux.
* __stellar_SN__ (FLOAT): The signal-to-noise value of the stellar continuum.
* __v_stars__ (FLOAT): The recessional velocity of the stellar continuum component along the line of sight. Measured in km/s with respect to the centre of the galaxy.
* __e_v_stars__ (FLOAT): Uncertainty on the recessional velocity. Measured in km/s.
* __sigma_stars__ (FLOAT): The velocity dispersion of the stellar continuum component (i.e. how broad the absorption lines are). Measured in km/s.
* __e_sigma_stars__ (FLOAT): Uncertainty on the velocity dispersion. Measured in km/s.
* __h3__ (FLOAT): Third moment of the line-of-sight velocity distribution. Unitless number. 
* __e_h3__ (FLOAT): Uncertainty on the `h3` measurement.
* __h4__ (FLOAT): Fourth moment of the line-of-sight velocity distribution. Unitless number.
* __e_h4__ (FLOAT): Uncertainty on the `h4` measurement.
* __v_gas__ (FLOAT): Recessional velocity of the single Gaussian emission line fit to the H $\alpha$ emission line (I _think_- it might be the average velocity of all the gas components). Measured in km/s with respect to the galaxy's centre.
* __e_v_gas__ (FLOAT): Uncertainty in the gas recessional velocity. Measured in km/s.
* __sigma_gas__ (FLOAT): Velocity dispersion of the H $\alpha$ emission line component (although see caveat above). Measured in km/s.
* __e_sigma_gas__ (FLOAT): Uncertainty in the velocity dispersion. Measured in km/s.
* __sfr__ (FLOAT): Star-formation rate of the spaxel calculated from the H $\alpha$ flux. See [Medling et al. 2018](https://ui.adsabs.harvard.edu/abs/2018MNRAS.475.5194M/abstract). Measured in solar masses / year.
* __e_sfr__ (FLOAT): Uncertainty in SFR, measured in solar masses / year.
* __sfr_density__ (FLOAT): Star-formation rate surface density of the spaxel, equal to SFR divided by the surface area over which it is measured. Measured in solar masses / year / kpc. 
* __e_sfr_density__ (FLOAT): Uncertainty in the star-formation rate surface density. 
* __halpha_flux__ (FLOAT): Flux of the H $\alpha$ emission line, found from a _single_ Gaussian component fit. Measured in $10^{-16} \mathrm{erg\,s}^{-1} \mathrm{cm}^{-2} \mathrm{A}^{-1}$.
* __e_halpha_flux__ (FLOAT): Uncertainty of the H $\alpha$ flux. Measured in $10^{-16} \mathrm{erg\,s}^{-1} \mathrm{cm}^{-2} \mathrm{A}^{-1}$.
* __hbeta_flux__ (FLOAT): Flux of the H $\beta$ emission line, found from a _single_ Gaussian component fit. Measured in $10^{-16} \mathrm{erg\,s}^{-1} \mathrm{cm}^{-2} \mathrm{A}^{-1}$.
* __e_hbeta_flux__ (FLOAT): Uncertainty of the H $\beta$ flux. Measured in $10^{-16} \mathrm{erg\,s}^{-1} \mathrm{cm}^{-2} \mathrm{A}^{-1}$.
* __NII6583_flux__ (FLOAT): Combined flux of the two Nitrogen emission lines at 6549.86 $\mathrm{A}$ and 6585.27 $\mathrm{A}$. The flux is derived from a _single_ Gaussian component fit to each line (with their ratios fixed according to that from atomic physics) and then summed. Measured in $10^{-16} \mathrm{erg\,s}^{-1} \mathrm{cm}^{-2} \mathrm{A}^{-1}$.
* __e_NII6583_flux__ (FLOAT): Uncertainty on the NII flux. Measured in $10^{-16} \mathrm{erg\,s}^{-1} \mathrm{cm}^{-2} \mathrm{A}^{-1}$.
* __OI6300_flux__ (FLOAT): Flux of the \[OI] emission line at 6300 $\mathrm{A}$, found from a _single_ Gaussian component fit. Measured in $10^{-16} \mathrm{erg\,s}^{-1} \mathrm{cm}^{-2} \mathrm{A}^{-1}$.
* __e_OI6300_flux__ (FLOAT): Uncertainty on the \[OI] emission line flux. Measured in $10^{-16} \mathrm{erg\,s}^{-1} \mathrm{cm}^{-2} \mathrm{A}^{-1}$.
* __OII3728_flux__ (FLOAT): Combined flux of the \[OII] doublet emission lines at 3729.875 $\mathrm{A}$ and 3729.875 $\mathrm{A}$. The flux is derived from a _single_ Gaussian component fit to each line and then summed. Measured in $10^{-16} \mathrm{erg\,s}^{-1} \mathrm{cm}^{-2} \mathrm{A}^{-1}$.
* __e_OII3728_flux__ (FLOAT): Uncertainty on the \[OII] emission line flux. Measured in $10^{-16} \mathrm{erg\,s}^{-1} \mathrm{cm}^{-2} \mathrm{A}^{-1}$.
* __OIII5007_flux__ (FLOAT): Combined flux of the \[OIII] doublet emission lines at 4960.295 $\mathrm{A}$ and 5008.240 $\mathrm{A}$. The flux is derived from a _single_ Gaussian component fit to each line and then summed. Measured in $10^{-16} \mathrm{erg\,s}^{-1} \mathrm{cm}^{-2} \mathrm{A}^{-1}$.
* __e_OIII5007_flux__ (FLOAT): Uncertainty on the \[OIII] emission line flux. Measured in $10^{-16} \mathrm{erg\,s}^{-1} \mathrm{cm}^{-2} \mathrm{A}^{-1}$.
* __SII6716_flux__ (FLOAT): Flux of the \[SII] emission line at 6716$\mathrm{A}$. The flux is derived from a _single_ Gaussian component fit. Measured in $10^{-16} \mathrm{erg\,s}^{-1} \mathrm{cm}^{-2} \mathrm{A}^{-1}$.
* __e_SII6716_flux__ (FLOAT). Uncertainty on the \[SII] emission line flux. Measured in $10^{-16} \mathrm{erg\,s}^{-1} \mathrm{cm}^{-2} \mathrm{A}^{-1}$.
* __SII6731_flux__ (FLOAT): Flux of the \[SII] emission line at 6731 $\mathrm{A}$. The flux is derived from a _single_ Gaussian component fit. Measured in $10^{-16} \mathrm{erg\,s}^{-1} \mathrm{cm}^{-2} \mathrm{A}^{-1}$.
* __e_SII6731_flux__ (FLOAT): Uncertainty on the \[SII] emission line flux. Measured in $10^{-16} \mathrm{erg\,s}^{-1} \mathrm{cm}^{-2} \mathrm{A}^{-1}$.
* __extinct_corr__ (FLOAT): The measured dust extinction in the spectrum, derived using the Balmer decrement (the ratio of H $\beta$ to H $\alpha$ flux). Units are magnitudes. 
* __e_extinct_corr__ (FLOAT): Uncertainty in the exitinction correction. 
* __log10_NII_Halpha__ (FLOAT): Log base 10 of the ratio between the \[NII] and H $\alpha$ fluxes.
* __log10_OI_Halpha__ (FLOAT): Log base 10 of the ratio between the \[OI] and H $\alpha$ fluxes.
* __log10_OIII_Hbeta__ (FLOAT): Log base 10 of the ratio between the \[OIII] and H $\beta$ fluxes.
* __log10_SII_Halpha__ (FLOAT): Log base 10 of the ratio between the \[SII] and H $\alpha$ fluxes.
* __BPT_class__ (INT): Classification of the spectrum into a BPT class of 0 (star-forming), 1 (a "Low-Ionisation Emission-line Region", due to sources like hot old stars, diffuse ionised gas, etc) or 2 (ionisation from an Active Galactic Nucleus). Can also be NaN, for spectra without a measurement.



## Using the database to get the data

Once you have access to `SAMI_spaxel_data.db`, you can interact with it like any other SQL database. If you use `python` and the `pandas` package, one way to interact with it is as follows:

```python
import sqlite3
import pandas as pd

con = sqlite3.connect("path/to/SAMI_spaxel_data.db")
# Note that this selects all rows and columns and will take ~1 minute
# You can use any valid SQL statement here, and the results are returned into a pandas dataframe.
df = pd.read_sql("select * from sami", con)
```

Or in `R`:

```R
library(tidyverse)
library(dbplyr)
library(RSQLite)

database  <- "/path/to/SAMI_spaxel_data.db"
con  <-  DBI::dbConnect(RSQLite::SQLite(), database)

# This 'lazily' loads, so should be very quick
sami  <- tbl(con, 'sami')
```

You can then interact with the `sami` table as you would a normal tibble. The `R` package `dbplyr` cleverly translates your usual `tidyverse` commands into `SQL`, then only queries the database right at the end. See the introduction [here](https://dbplyr.tidyverse.org/articles/dbplyr.html).

### Plotting a 2D map 

To plot a 2D map, we can select all spaxels from a single galaxy (by its CATID) and then reshape those values into an image. All SAMI cubes have 50 pixels in the $x$ and $y$ directions, which makes our task simple:

```python
# Pick a random galaxy
CATID = 623679

galaxy = pd.read_sql(f"select * from sami where CATID == {CATID}", con)
# or, if we've already got the master dataframe from above:
# galaxy = df.loc[df.CATID == CATID]

# Plot the H-alpha flux for this galaxy, alongside its BPT classification for each pixel
import matplotlib as mpl
import matplotlib.pyplot as plt
fig, axs = plt.subplots(ncols=2, figsize=(14, 5))

# Plot the flux map
axs[0].imshow(galaxy['halpha_flux'].values.reshape(50, 50), cmap=cmap)

# Make stuff for the BPT colorbar
cmap = mpl.cm.plasma
bounds = [0, 1, 2, 3]
norm = mpl.colors.BoundaryNorm(bounds, cmap.N, extend='neither')

# Plot the classes
im = axs[1].imshow(galaxy['BPT_class'].values.reshape(50, 50), cmap=cmap, norm=norm)

# Add some extra titles and a colorbar
fig.colorbar(im, ax=axs[1], label='BPT Class')
fig.suptitle(f"{CATID}")
axs[0].set_title("H-alpha Flux")
axs[1].set_title("BPT Classification")
```

which gives this:

![example H-alpha flux map for galaxy 623679](resources/images/example_image.png)

## Getting rid of NaNs 

You'll notice that a large chunk of rows for each galaxy are just NaNs. This is due to the fact that the reduced data for each galaxy square arrays but the footprint of the SAMI instrument is circular (hexagonal, actually, but the final galaxies have round cutouts after the dithering and stacking of 7 observations). The data is therefore padded with NaNs to make it fit into a square allocation in memory.

To drop these NaN values, you can use the following:

```python
# Select the columns you want to drop NaN values from
df_no_NaNs = df.dropna(subset=['halpha_flux', 'hbeta_flux', 'BPT_classification', 'etc...']
```

However you'll see that you _can't_ easily plot a 2D image from this dataframe, because the spectra for each individual galaxies no longer make a nice square `numpy` array:

```python
# This won't work!
plt.imshow(df_no_NaNs[df_no_NaNs['CATID'] == CATID, 'halpha_flux'].values.reshape(50, 50))
```
```bash
ValueError: cannot reshape array of size 1234 into shape (50,50)
```

What you can do is plot the RA/DEC or $x$ and $y$ positions of each spectrum, and colour the points based on some other attribute:

```python
galaxy = df_no_NaNs[df_no_NaNs['CATID'] == CATID]
# This will work
plt.scatter(galaxy['RA'], galaxy['DEC'], c=galaxy['halpha_flux'])
```
