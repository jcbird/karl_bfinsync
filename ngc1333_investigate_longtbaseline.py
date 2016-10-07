import astropy.io.fits as pyfits
import corner
import numpy as np
from astropy.table import Table


data = Table.read('./per_epoch_data_file.fits').to_pandas()

print('hey')

#ngc1333 = Table.read(data[data['cluster']=='NGC1333']).to_pandas()


bycluster = data.groupby('cluster')

ngc1333 = bycluster.get_group('NGC1333')
ngc1333_bystar = ngc1333.groupby('twomassid')


med_teff = ngc1333_bystar['teff'].median()
med_logg = ngc1333_bystar['logg'].median()
med_vsini = ngc1333_bystar['vsini'].median()

ngc1333['dev_teff'] = ngc1333.apply(lambda X: np.abs(X['teff'] - med_teff[X['twomassid'].values]) )

ngc1333['dev_logg'] = ngc1333.apply(lambda X: np.abs(X['logg'] - med_logg[X['twomassid']]))
ngc1333['dev_vsini'] = ngc1333.apply(lambda X: np.abs(X['vsini'] - med_vsini[X['twomassid']]))

fig, axs = plt.subplots(6,7, figsize=(15,10))
aaxs = np.ravel(axs)
for ii, (key, grp) in enumerate(n1333_bystar):
    if ii<42:
        aaxs[ii].scatter(grp['JD']-2450000. , grp['corr_vrad'], c=grp['dev_teff'], vmax=150., vmin=50., cmap='coolwarm', s=30)
        aaxs[ii].plot(grp['JD']-2450000. , grp['corr_vrad'], c='k')
        aaxs[ii].set_title(key, fontsize=10)



drv = ngc1333_bystar['corr_vrad'].apply(lambda X: X.max() - X.min())

drv = ngc1333_bystar['corr_vrad'].apply(lambda X: X.max() - X.min())

timebaseline = ngc1333_bystar['JD'].apply(lambda X: X.max() - X.min())

max_Teff_resid = ngc1333_bystar['teff'].apply(lambda X: np.max(np.abs(X - X.median())))
max_logg_resid = ngc1333_bystar['logg'].apply(lambda X: np.max(np.abs(X - X.median())))
max_vsini_resid = ngc1333_bystar['vsini'].apply(lambda X: np.max(np.abs(X - X.median())))


mad_Teff_resid = ngc1333_bystar['teff'].apply(lambda X: np.median(np.abs(X - X.median())))
mad_logg_resid = ngc1333_bystar['logg'].apply(lambda X: np.median(np.abs(X - X.median())))
mad_vsini_resid = ngc1333_bystar['vsini'].apply(lambda X: np.median(np.abs(X - X.median())))

rvmax_Nmad_Teff = ngc1333_bystar.apply(lambda X: np.abs(X['teff'][X['corr_vrad'].argmax()] - np.median(X['teff']))/(np.median(np.abs(X['teff'] - X['teff'].median()))))
rvmin_Nmad_Teff = ngc1333_bystar.apply(lambda X: np.abs(X['teff'][X['corr_vrad'].argmin()] - np.median(X['teff']))/(np.median(np.abs(X['teff'] - X['teff'].median()))))


rvmax_Nmad_logg = ngc1333_bystar.apply(lambda X: np.abs(X['logg'][X['corr_vrad'].argmax()] - np.median(X['logg']))/(np.median(np.abs(X['logg'] - X['logg'].median()))))
rvmin_Nmad_logg = ngc1333_bystar.apply(lambda X: np.abs(X['logg'][X['corr_vrad'].argmin()] - np.median(X['logg']))/(np.median(np.abs(X['logg'] - X['logg'].median()))))

rvmax_Nmad_vsini = ngc1333_bystar.apply(lambda X: np.abs(X['vsini'][X['corr_vrad'].argmax()] - np.median(X['vsini']))/(np.median(np.abs(X['vsini'] - X['vsini'].median()))))

rvmin_Nmad_vsini = ngc1333_bystar.apply(lambda X: np.abs(X['vsini'][X['corr_vrad'].argmin()] - np.median(X['vsini']))/(np.median(np.abs(X['vsini'] - X['vsini'].median()))))

dd_A = np.column_stack( (drv, rvmax_Nmad_logg, rvmin_Nmad_logg, rvmax_Nmad_Teff, rvmin_Nmad_Teff, rvmax_Nmad_vsini, rvmin_Nmad_vsini, timebaseline))

corner.corner(dd_A, labels=['drv', 'max_Nmad_logg', 'min_Nmad_logg', 'max_Nmad_teff', 'min_Nmad_teff', 'max_Nmad_vsini', 'min_Nmad_vsini', 'time'], show_titles=True, quantiles=[.16, .5, .84])

fracmad_Teff_resid = ngc1333_bystar['teff'].apply(lambda X: np.median(np.abs(X - X.median()))/X.median())
fracmad_logg_resid = ngc1333_bystar['logg'].apply(lambda X: np.median(np.abs(X - X.median()))/X.median())
fracmad_vsini_resid = ngc1333_bystar['vsini'].apply(lambda X: np.median(np.abs(X - X.median()))/X.median())
