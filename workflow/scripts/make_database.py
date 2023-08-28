import sqlite3

smk = snakemake  # noqa

conn = sqlite3.connect(smk.output.database)
c = conn.cursor()

c.execute(
    """DROP TABLE IF EXISTS sami
        """
)
c.execute(
    """
        CREATE TABLE IF NOT EXISTS sami
        ([id] INTEGER PRIMARY KEY,
        [CATID] INTEGER,
        [RA] FLOAT,
        [DEC] FLOAT,
        [cube_x] INTEGER,
        [cube_y] INTEGER,
        [stellar_flux] FLOAT,
        [e_stellar_flux] FLOAT,
        [stellar_SN] FLOAT,
        [v_stars] FLOAT,
        [e_v_stars] FLOAT,
        [sigma_stars] FLOAT,
        [e_sigma_stars] FLOAT,
        [h3] FLOAT,
        [e_h3] FLOAT,
        [h4] FLOAT,
        [e_h4] FLOAT,
        [v_gas] FLOAT,
        [e_v_gas] FLOAT,
        [sigma_gas] FLOAT,
        [e_sigma_gas] FLOAT,
        [sfr] FLOAT,
        [e_sfr] FLOAT,
        [sfr_density] FLOAT,
        [e_sfr_density] FLOAT,
        [halpha_flux] FLOAT,
        [e_halpha_flux] FLOAT,
        [hbeta_flux] FLOAT,
        [e_hbeta_flux] FLOAT,
        [NII6583_flux] FLOAT,
        [e_NII6583_flux] FLOAT,
        [OI6300_flux] FLOAT,
        [e_OI6300_flux] FLOAT,
        [OII3728_flux] FLOAT,
        [e_OII3728_flux] FLOAT,
        [OIII5007_flux] FLOAT,
        [e_OIII5007_flux] FLOAT,
        [SII6716_flux] FLOAT,
        [e_SII6716_flux] FLOAT,
        [SII6731_flux] FLOAT,
        [e_SII6731_flux] FLOAT,
        [extinct_corr] FLOAT,
        [e_extinct_corr] FLOAT,
        [log10_NII_Halpha] FLOAT,
        [log10_OI_Halpha] FLOAT,
        [log10_OIII_Hbeta] FLOAT,
        [log10_SII_Halpha] FLOAT,
        [BPT_class] INT)
        """
)

conn.commit()
