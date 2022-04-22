#!/usr/bin/env python
# To check average Qvalue
# By Yusuke Satoh
# On

import sys
import os
import gc
import numpy as np
import itertools
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.gridspec as gridspec
from multiprocessing import Pool
from netCDF4 import Dataset
from utiltools import day2month, fillinMasked, extractMasked, get_nearest_value

DUMMY = False
TEST = False
suffix = 'pdf'

#ghms = ['cwatm', 'h08', 'lpjml','matsiro', 'watergap2']
#gcms = ['hadgem2-es','ipsl-cm5a-lr','gfdl-esm2m','miroc5']
ghms = ['cwatm']
gcms = ['hadgem2-es']
# reference period of drought threshold
syear, eyear = 1861, 2005
years = range(syear, eyear+1)
nyear = len(years)
qvale_directory_main = '/data/rg001/sgec0017/data/isimip2b.drought'
dis_directory_main = '/data/rg001/sgec0017/data/isimip2b/out/nc/water_global.processed'
fig_directory_main = '/data/rg001/sgec0017/data/figure_box/drought_tpcd_isimip2b/estimate_average_Qvalue_for_reviewer3'

lndmskPath = '/data/rg001/sgec0017/data/mapmask/ISIMIP2b_landseamask_generic.nc4'
lndmask = np.ma.masked_not_equal(Dataset(lndmskPath).variables['LSM'][:][0], 1.0).mask
bm = Basemap(projection='cyl', llcrnrlat=-56.5, urcrnrlat=84.5, llcrnrlon=-180., urcrnrlon=180., resolution='l')
scn, soc, co2 = 'historical', 'histsoc', 'co2'

if TEST:
    ghms = ['cwatm']
    gcms = ['hadgem2-es']
    years = years[:3]
    nyear = len(years)

# ----------------------------------------------------------------------------------------------------------------------
def load_qvalue(ghm, gcm):
    srcpath = os.path.join(qvale_directory_main, ghm, gcm,
                           f'Qvalue_historical_histsoc_co2.{syear}_{eyear}',
                           'Qvalues', 'Q80win15.nc4'
                           )
    src = Dataset(srcpath)['Qvalue'][:]  # (365, ny, nx)
    print(f'load: {srcpath} {src.shape}')
    return extractMasked(src, lndmask)  # (365, nland)


def load_dis(scn, soc, co2, ghm, gcm, year):
    ncFileName = f'{ghm}_{gcm}_ewembi_{scn}_{soc}_{co2}_dis_global_daily_{year}.nc'
    srcPath = os.path.join(dis_directory_main, ghm, gcm, scn, ncFileName)
    if not os.path.isfile(srcPath): srcPath = srcPath[:-2]+'nc4'
    if not os.path.isfile(srcPath):
        raise FileNotFoundError(f'Error!! {srcPath} is not exist... Check!!')
    else:
        ncFile = Dataset(srcPath)
        aSrc = ncFile.variables['dis'][:]
        shp = aSrc.shape
        if shp[0] == 366:
            aSrc = np.concatenate((aSrc[:59], aSrc[60:]), axis=0)
        print(f'loaded: {srcPath} {shp}')  #, type(aSrc)
        return extractMasked(aSrc, lndmask)  # Extract only land  (nday, nland)


# ----------------------------------------------------------------------------------------------------------------------
def draw_maps(ghm, gcm, src):

    shp = src.shape
    src = np.ma.masked_equal(src, 1e+20)
    src = np.ma.masked_array(src, mask=np.resize(lndmask, src.shape))
    if len(shp) == 3:
         src = src[:, 11:293, :]
    elif len(shp) == 4:
         src = src[:, :, 11:293,:]

    vmin, vmax = 75, 100
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
    cmap = mpl.cm.rainbow


    for stat in ['average', 'min', 'max', 'hist']:

        if len(shp) == 3:
            if stat == 'average':
                _src = src.mean(axis=0)
            elif stat == 'min':
                _src = src.min(axis=0)
            elif stat == 'max':
                _src = src.max(axis=0)
            elif stat == 'hist':
                _src = src.compressed()
        elif len(shp) == 4:
            if stat == 'average':
                _src = src.mean(axis=(0,1))
            elif stat == 'min':
                _src = src.min(axis=1).mean(axis=0)
            elif stat == 'max':
                _src = src.max(axis=1).mean(axis=0)
            elif stat == 'hist':
                _src = src.compressed()
        else:
            raise ValueError
        if stat == 'hist': figsize, left, right, bottom, top = (5, 2.5), 0.14, 0.98, 0.09, 0.98
        else:              figsize, left, right, bottom, top = (6, 2.5), 0.01, 0.99, 0.05, 0.99
        print(f'_src.shape: {_src.shape}')

        fig = plt.figure(figsize=figsize, dpi=300)
        plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top)  #, wspace=0.1, hspace=0.15)
        #fig.suptitle(f'Q80\n{ghm} {gcm}')
        #gs = gridspec.GridSpec(nrows=1, ncols=1, figure=fig)
        #ax = fig.add_subplot(gs[0,0])
        ax = fig.add_subplot(111)
	
        if stat == 'hist':  # histgram
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            edges = np.arange(75, 102) - 0.5
            ax.hist(_src, bins=edges, color='skyblue', edgecolor='black', linewidth=0.5, cumulative=False, density=True)
            ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))
            ax.set_ylabel('Probability')
            #ax.set_xlabel('Fraction [%]')
            ax.set_xlabel('Re-estimated Qx')
            ax.set_xlim([74.5, 100.5])
        elif stat in ['average', 'min', 'max']:
            ax.axis('off')
            #ax.set_title(f'{stat}')
            bm.imshow(_src[::-1], vmin=vmin, vmax=vmax, cmap=cmap)
            #bm.drawcoastlines(linewidth=0.05, color='#808080')
            ax_pos = ax.get_position()
            ax11 = fig.add_axes([ax_pos.x0 + 0.01, ax_pos.y0 + 0.05, 0.25, 0.02])
            cb1 = mpl.colorbar.ColorbarBase(ax11, cmap=cmap, norm=norm, orientation='horizontal')
            cb1.outline.set_visible(False)
            cb1.ax.set_title(f'{stat} [%]')
        else:
            raise ValueError

        figname = f'{ghm}.{gcm}.relative_qvalue.{stat}.{suffix}'
        fig_directory = os.path.join(fig_directory_main)
        if not os.path.isdir(fig_directory): os.makedirs(fig_directory)
        figpath = os.path.join(fig_directory, figname)
        plt.savefig(figpath, dpi=300)
        print(f'savefig: {figpath}')
        plt.close(fig)
        gc.collect()


# ----------------------------------------------------------------------------------------------------------------------
def main():
    """
    The comment by the reveiwer#3:
    I would simply recommend assessing a posteriori the average (over the cells and days of the year) of the probability of exceedance
    over the reference period (1861-2005) that actually derives from the chosen calculations with a Q80 over a 31-day moving window. 
    The question is basically : is it more like 85 or 95?
    """

    nland = np.where(lndmask==0, 1, 0).sum()
    print(f'nland: {nland}')
    ensemble = []
    for ghm, gcm in itertools.product(ghms, gcms):

        if DUMMY:
            output = np.random.randint(75, 101, (365, nland))
        else:
            # load qvalue used
            qvalue = np.ma.masked_greater(load_qvalue(ghm, gcm), 1e+18)  # (365, nland)
            # load daily discharge
            discharge = np.ma.masked_greater(np.array([load_dis(scn, soc, co2, ghm, gcm, year) for year in years]), 1e+18)  # (nYear, 365, nland)

            # calc daily return period of qvalue (considering only year-to-year variability over 31-day moving window river discharge time series)
            window_size = 31
            weights = np.ones(window_size) / window_size
            discharge = discharge.reshape(-1, nland)                                                                                    # (nyear*365, nland)
            discharge = np.concatenate((discharge[-15:], discharge, discharge[:15]), axis=0)                                            # (15+nyear*365+15, nland) 
            discharge = np.array([np.convolve(timeseries_at_a_grid, weights, mode='valid') for timeseries_at_a_grid in discharge.T]).T  # (nyear*365, nland)
            discharge = discharge.reshape(-1, 365, nland)                                                                               # (nYear, 365, nland)
            
            output = np.full((365, nland), 1e+20)
            for iland, iday in itertools.product(range(nland), range(365)):
                dis = discharge[:, iday, iland]  # (nYear)
                dis.sort()                       # from smaller to larger
                dis = dis[::-1]                  # from larget to smaller
                nearest_dis = get_nearest_value(dis, qvalue[iday, iland])
                if nearest_dis != np.nan:
                    output[iday, iland] = dis.tolist().index(nearest_dis) / float(nyear)  * 100
                    #output[iday, iland] = (1 - dis.tolist().index(nearest_dis) / float(nyear)) * 100  # just in case there are same values for several years 
                if np.mod(iland, 1000) == 0 and iday == 0: print(f' iland = {iland}')
            del qvalue, discharge
            gc.collect()

        draw_maps(ghm, gcm, fillinMasked(output, lndmask))
        sys.exit()
        
        ensemble.append(output)

    ensemble = np.array(ensemble)  # (nmember, 365, nland)
    draw_maps('ensemble', 'info', fillinMasked(ensemble, lndmask))


if __name__ == '__main__':
    main()
