from tools import apvisit_interface
import numpy as np
from os import path
from astropy import table
from matplotlib import pyplot as plt
import pandas as pd

_basepath = path.join(path.expanduser('~'), 'projects', 'karl')

insync_file = path.join(_basepath,
                        'insync_apogee_all_star_masked_pleiades_data.csv')
insync_df = table.Table.read(insync_file).to_pandas()


dr13_file = path.join(path.expanduser('~'), 'obsdata',
                      'dr13/apogee/spectro/redux/r6/stars/l30e/l30e.2',
                      'allVisit-l30e.2.fits')
# Default is dr14 file
dr14_allVisit = apvisit_interface.prep_allVisit()
dr13_allVisit = apvisit_interface.prep_allVisit(dr13_file)

insync_rv = insync_df['vrad_epoch'].values
insync_twomassid = insync_df['twomassid_epoch'].values
insync_jd = insync_df['raw_jd'].values

dr14_synthvhelio = np.empty_like(insync_rv)
dr14_vhelio = np.empty_like(insync_rv)
dr14_obsvhelio = np.empty_like(insync_rv)

dr14_flags = []

dr14rows = []
dr13rows = []

for ii, (twomassid, jd) in enumerate(zip(insync_twomassid, insync_jd)):
    dr14cat_row = apvisit_interface.find_apogee_match(twomassid, jd,
                                                      dr14_allVisit)
    dr13cat_row = apvisit_interface.find_apogee_match(twomassid, jd,
                                                      dr13_allVisit)
    dr14rows.append(dr14cat_row)
    dr13rows.append(dr13cat_row)

dr14_data = pd.concat(dr14rows)
dr13_data = pd.concat(dr13rows)
dr14_data = dr14_data.replace(999999.0, np.nan)
dr13_data = dr13_data.replace(999999.0, np.nan)


def plot_comparison(is_rv, ap_rv, rv_error=None, cb_label=None):
    fig, ax = plt.subplots()
    xx = np.linspace(-20,10,50)
    if rv_error is None:
        ax.scatter(is_rv, ap_rv, s=25, alpha=0.7)
    else:
        f = ax.scatter(is_rv, ap_rv, c=rv_error, vmin=0, vmax=1, s=25,
                       alpha=0.7)
        cb = plt.colorbar(f)
        cb.set_label(cb_label)
    ax.plot(xx,xx, lw=1)
    ax.set_ylim(0,8)
    ax.set_xlim(0,8)
    return fig, ax


def NDRV_allVisit(dr_data):
    bystar = dr_data.groupby(level='APOGEE_ID')
    maxRV_idx = bystar.VHELIO.idxmax()
    minRV_idx = bystar.VHELIO.idxmin()
    deltaRV = (dr_data.loc[maxRV_idx]['VHELIO'].values -
               dr_data.loc[minRV_idx]['VHELIO'].values)
    print(deltaRV)
    normDRV = np.sqrt((0.1 + dr_data.loc[maxRV_idx]['VRELERR'].values)**2 +
                      (0.1 + dr_data.loc[minRV_idx]['VRELERR'].values)**2)
    print(normDRV)
    return deltaRV / (normDRV)


def NDRV_insync(insync_data):
    bystar = insync_data.groupby('twomassid_epoch')
    maxRV_idx = bystar.vrad_epoch.idxmax()
    minRV_idx = bystar.vrad_epoch.idxmin()
    deltaRV = (insync_data.loc[maxRV_idx]['vrad_epoch'].values -
               insync_data.loc[minRV_idx]['vrad_epoch'].values)
    print(deltaRV)
    normDRV = np.sqrt(insync_data.loc[maxRV_idx]['dvrad_epoch'].values**2 +
                      insync_data.loc[minRV_idx]['dvrad_epoch'].values**2)
    print(normDRV)
    return deltaRV / normDRV

if __name__ == '__main__':
    plot_dir = 'plots'

    cblabel = "DR14 VRELERR [km/s]"
    fig, ax = plot_comparison(insync_rv, dr14_data['VHELIO'],
                              dr14_data['VRELERR'], cb_label=cblabel)
    ax.set_xlabel('INSYNC RV [km/s]')
    ax.set_ylabel('DR14 VHELIO [km/s]')
    ax.set_title('DR14 INSYNC RV Comparison')
    basename = path.join(plot_dir, 'RV_insync_vs_dr14')
    plt.savefig(basename + '.png', format='png')
    plt.savefig(basename + '.eps', format='eps')

    cblabel = "DR13 VRELERR [km/s]"
    fig, ax = plot_comparison(insync_rv, dr13_data['VHELIO'],
                              dr13_data['VRELERR'], cb_label=cblabel)
    ax.set_xlabel('INSYNC RV [km/s]')
    ax.set_ylabel('DR13 VHELIO [km/s]')
    ax.set_title('DR13 INSYNC RV Comparison')
    basename = path.join(plot_dir, 'RV_insync_vs_dr13')
    plt.savefig(basename + '.png', format='png')
    plt.savefig(basename + '.eps', format='eps')

    cblabel = "DR14 VRELERR [km/s]"
    fig, ax = plot_comparison(dr14_data['VHELIO'], dr13_data['VHELIO'],
                              dr14_data['VRELERR'], cb_label=cblabel)
    ax.set_xlabel('DR14 VHELIO [km/s]')
    ax.set_ylabel('DR13 VHELIO [km/s]')
    ax.set_title('DR13 vs. DR14 Pleiades RVs')
    basename = path.join(plot_dir, 'RV_dr14_vs_dr13')
    plt.savefig(basename + '.png', format='png')
    plt.savefig(basename + '.eps', format='eps')

    cblabel = "DR14 VRELERR [km/s]"
    fig, ax = plot_comparison(dr14_data['VHELIO'] - dr13_data['VHELIO'],
                              dr14_data['RV_LOGG'] - dr13_data['RV_LOGG'],
                              dr14_data['VRELERR'], cb_label=cblabel)
    ax.set_xlabel('DR14 - DR13 VHELIO [km/s]')
    ax.set_ylabel('DR14 - DR13 RV_LOGG')
    ax.set_title(r'DR13, DR14 $\Delta$RV vs. $\Delta$RV_LOGG')
    ax.set_xlim(-6,6)
    ax.set_ylim(-2.1,2.1)
    basename = path.join(plot_dir, 'RV_dr14_vs_dr13_RVLOGG')
    plt.savefig(basename + '.png', format='png')
    plt.savefig(basename + '.eps', format='eps')

    fig, ax = plot_comparison(dr14_data['VHELIO'] - dr13_data['VHELIO'],
                              dr14_data['CHISQ'] - dr13_data['CHISQ'])
    ax.set_xlabel('DR14 - DR13 VHELIO [km/s]')
    ax.set_ylabel('DR14 - DR13 CHISQ')
    ax.set_title(r'DR13, DR14 $\Delta$RV vs. $\Delta$CHISQ')
    ax.set_xlim(-6,6)
    ax.set_ylim(-6,6)
    basename = path.join(plot_dir, 'RV_dr14_vs_dr13_CHISQ')
    plt.savefig(basename + '.png', format='png')
    plt.savefig(basename + '.eps', format='eps')

    # import seaborn as sns
    # sns.kdeplot(dr14_data['VHELIO'], clip=(-25,15), label='DR14')
    # sns.kdeplot(dr13_data['VHELIO'], clip=(-25,15), label='DR13')
    # plt.xlabel('VHELIO [km/s]')
    # basename = path.join(plot_dir, 'RV_PDFs_dr14_vs_dr13')
    # plt.savefig(basename + '.png', format='png')
    # plt.savefig(basename + '.eps', format='eps')
