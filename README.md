# Code to add SAMI dataproducts into a simple SQLite database

This repo contains code which collates various SAMI data products into a single SQLite database, which will be used for the COIN residence program. 

This repo _doesn't_ contain the individual SAMI dataproducts or the final database, since they're each a few Gb in size. The data products can be downloaded from [the bulk download tool at Data Central](https://datacentral.org.au/services/download/) and the database itself can be found on the COIN server (location tbd).If you've downloaded the SAMI data, you can run all the scripts and build the database using [snakemake](https://snakemake.readthedocs.io/en/stable/) but I don't imagine that this will be necessary!

## Scripts

The scripts to build the database are in `workflow/scripts`. They do the following:

    * `make_database.py` sets up a new SQLite database with the appropriate columns

    * `add_galaxy_data_to_db.py` reads in the various gas and star measurements for a single SAMI galaxy, adds some extra columns (i.e. calculates some important emission line ratios and adds the BPT classification for each spaxel) and then adds this data to the master database.
    
    * `utils.py` contains utility functions which are used by `add_galaxy_data_to_db.py`. In particular, it contains the functions which do the BPT classifications based on different emission line ratios. It also contains functions which read in the raw data, as well as a few things I've been playing around with (and aren't used yet). 

## Using the database to get the data




