import numpy as np
from os import path
from astropy import table
import pandas as pd

_basepath = path.join(path.expanduser('~'), 'projects', 'karl')


def prep_allVisit(apogee_file=None):
    if apogee_file is None:
        apogee_file = path.join(path.expanduser('~'), 'obsdata', 'dr14_beta',
                                'allVisit-l31c.1.fits')
    apogee_df = table.Table.read(apogee_file).to_pandas()
    apogee_df['APOGEE_ID'] = pd.Series([i.strip() for
                                        i in apogee_df['APOGEE_ID']],
                                       index=apogee_df.index)
    apogee_df.set_index(['APOGEE_ID', 'JD'], inplace=True)
    return apogee_df


def find_apogee_match(twomassid, JD, apogee_cat):
    """
    Given 2MASS ID and JD of INSYNC observation, find corresponding
    visit data in apogee catalog.

    Parameters
    ----------------

    twomassid : string
        2MASS ID string
    JD : float
        Julian Date of observation
    apogee_cat : dataframe
        APOGEE allVisit data in pandas dataframe
    Returns
    ----------------
    Row of APOGEE allVisit file corresponding to input data.

    """
    # Get JD of rows corresponding to star

    obsdates = apogee_cat.loc[twomassid].index.values
    indx, delta_date = min(enumerate(np.abs(obsdates - JD)),
                           key=lambda X: X[1])
    if delta_date < 0.009:
        apogee_JD = obsdates[indx]
    else:
        print('INSYNC data: {} does not match'.format(JD))
        print('APOGEE dates: {}'.format(obsdates))
        return None

    return apogee_cat.loc[twomassid, apogee_JD]
