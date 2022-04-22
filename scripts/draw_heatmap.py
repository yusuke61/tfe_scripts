import os
import sys
import itertools
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker, cm
from mpl_toolkits.axes_grid1 import make_axes_locatable
from netCDF4 import Dataset
from utiltools import get_region_info
plt.rcParams['font.family'] = 'Arial'

experiment = 'basic'
#experiment = 'picontrol'
#experiment = 'co2exp'
#experiment = 'rcp26soc'


sample_type = 'Mean'
#sample_type = 'Average'

#region_type = 'AR6_regions'
region_type = 'HydroBASINS_lev1'

#process_type = 'original'
process_type = 'bootstrap'

trend_type = 'quadratic'

#trend_syear = 2005
#trend_syear = 1995
trend_syear = 1985
trend_syear = str(trend_syear)

#K = 1000
K = 100000

blocksize = 5

tChunk = 5
case = 'case_1'

# --- basic
rcps = ['rcp26', 'rcp85']
socs = ['2005soc']
co2s = ['co2']
gcms = ['hadgem2-es', 'ipsl-cm5a-lr', 'gfdl-esm2m', 'miroc5']  # full members
ghms = ['cwatm', 'h08', 'lpjml', 'matsiro', 'watergap2']  # full members
# --- experimental
if experiment == 'picontrol':
    rcps = ['picontrol', 'rcp85', 'rcp26']
elif experiment == 'co2exp':
    ghms = ['lpjml', 'matsiro']
    co2s = ['co2', '2005co2']
elif experiment == 'rcp26soc':
    ghms = ['cwatm', 'h08', 'lpjml', 'matsiro']
    rcps = ['rcp26']
    socs = ['2005soc', 'rcp26soc']

#seasons = ['DRY', 'WET', 'annual']
#seasons = ['DRY', 'annual']
seasons = ['DRY']


ticklabels_year = range(2020, 2100, 10)

suffixes = ['png', 'pdf']
dpi_figure = 80
dpi_savefig = 300

most_parent_directory = '/data/rg001/sgec0017/data/figure_box/drought_tpcd_isimip2b'
src_directory_parent = os.path.join(most_parent_directory, 'bootstrap')
figure_directory_parent = os.path.join(most_parent_directory, 'draw_tfe_heatmap')

# date of bootstrap
if region_type == 'AR6_regions':
    dict_input_date = {  
        ('Mean', 'basic',    'annual'): 20210707,
        ('Mean', 'basic',    'DRY'   ): 20210707,
        ('Mean', 'basic',    'WET'   ): 20210707,
        ('Mean', 'co2exp',   'DRY'   ): 20210708,
        ('Mean', 'rcp26soc', 'DRY'   ): 20210709,
        }
elif region_type == 'HydroBASINS_lev1':
    dict_input_date = {  # date of bootstrap
        ('Mean', 'basic',    'DRY'): 20220103,
        ('Mean', 'co2exp',   'DRY'): 20220103,
        ('Mean', 'rcp26soc', 'DRY'): 20220103,
        }
else:
    raise ValueError

# list of regions
region_map, dict_region = get_region_info(region_type, TEST=False)
regions = [f'{region_id:02}.{region_name}' for region_id, region_name in dict_region.items()]
_regions = [f'{region_id:02} {region_name}' for region_id, region_name in dict_region.items()]
del region_map


def draw_heatmap(aSrc, rcp, soc, co2, season):
    print('\nheatmap...')  # member heatmap of TPCD

    if experiment == 'co2exp': 
        width = 5
    elif experiment == 'rcp26soc': 
        width = 6.5
    else:  # basic
        width = 7
    if region_type == 'AR6_regions':
        height = 12
        left = 0.12
        gcm_ypos = 47
        yaxis_tick_pad = 65
    elif region_type == 'HydroBASINS_lev1':
        height = 14
        bottom = 0.14
        left = 0.26
        if experiment in ['basic', 'rcp26soc']:
            gcm_ypos = 67
        elif experiment == 'co2exp':
            gcm_ypos = 65
        yaxis_tick_pad = 134
    else:
        raise ValueError

    fig = plt.figure(num=10, figsize=(width,height), dpi=dpi_figure)
    plt.subplots_adjust(left=left, bottom=bottom, right=0.975, top=0.96)
    ax = fig.add_subplot(111)
    ax.set_title(f'{rcp} {soc} {co2}')

    im = ax.imshow(aSrc, vmin=2018, vmax=2098, cmap=cm.hot)
    
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.14)
    cbar = plt.colorbar(im, cax=cax, ticks=ticklabels_year, orientation='vertical')
    cbar.ax.tick_params(labelsize=8)

    xticklabels = ['{:<22}'.format('median')] + ['{:<7}'.format(ghm) for ghm in ghms] * len(gcms)
    xticks = np.arange(len(xticklabels))
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels, rotation=-90)
    for igcm, gcm in enumerate(gcms):
        if experiment == 'co2exp':
            ax.text(0.5+len(ghms)/2+len(ghms)*igcm, gcm_ypos+2, gcm, ha='center', rotation=-90)
        else:
            ax.text(0.5+len(ghms)/2+len(ghms)*igcm, gcm_ypos, gcm, ha='center', rotation=-90)
    if experiment == 'co2exp':
        minor_ticks = [0, 0.5, 2.5, 4.5, 6.5, 8.5]
    elif experiment == 'rcp26soc':
        minor_ticks = [0, 0.5, 4.5, 8.5, 12.5, 16.5]
    else:
        minor_ticks = [0, 0.5, 5.5, 10.5, 15.5, 20.5]
    ax.xaxis.set_minor_locator(ticker.FixedLocator(minor_ticks))
    #ax.xaxis.set_major_formatter(ticker.NullFormatter())
    #ax.xaxis.set_minor_locator(ticker.FixedLocator(np.arange(len(xticklabels))))
    #ax.xaxis.set_minor_formatter(ticker.FixedFormatter(xticklabels))
    #labels = ax.get_xticklabels(minor=True)
    #plt.setp(labels, fontsize=6, rotation=-90)

    ax.set_yticks(range(len(regions)))
    ax.set_yticklabels(_regions)
    for label in ax.yaxis.get_ticklabels():
            label.set_horizontalalignment('left')
    ax.yaxis.set_tick_params(pad=yaxis_tick_pad)
    #ax.yaxis.set_major_formatter(ticker.NullFormatter())
    #ax.yaxis.set_minor_locator(ticker.FixedLocator(np.arange(len(regions))))
    #ax.yaxis.set_minor_formatter(ticker.FixedFormatter([region for region in regions]))
    #labels = ax.get_yticklabels(minor=True)
    #plt.setp(labels, fontsize=6)

    ax.tick_params(axis='y',    which='major', direction='out', length=0.2, width=0.5, colors='black')
    ax.tick_params(axis='x',    which='major', direction='out', length=0.3)
    ax.tick_params(axis='both', which='minor', direction='out', length=70)

    ax.axvline(x=0.5, color='black', lw=0.5)
    for igcm, gcm in enumerate(gcms):
        if igcm != 0: ax.axvline(x=0.5+len(ghms)*igcm, color='black', ls='--', lw=0.5)

    figname = f'heatmap.{rcp}_{soc}_{co2}.{process_type}.{season}.'
    for suffix in suffixes:
        _fig_directory = os.path.join(figure_directory_parent, region_type, experiment)
        if not os.path.isdir(_fig_directory): os.makedirs(_fig_directory)
        figpath = os.path.join(_fig_directory, figname+suffix)
        plt.savefig(figpath, dpi=dpi_savefig)
        print(f'savefig: {figpath}')
    plt.close(10)


def main(*argv):

    for season, rcp, soc, co2 in itertools.product(seasons, rcps, socs, co2s):

        src_directory = os.path.join(src_directory_parent, region_type,
                                     f'{experiment}.{sample_type}.{trend_type}{trend_syear[2:]}99.{K}.block{blocksize}.{season}', 
                                     f'tChunk_{tChunk:02}', 
                                     str(dict_input_date[(sample_type, experiment, season)]), f'{rcp}_{soc}_{co2}')
        heatmap = np.full((len(regions), len(ghms)*len(gcms)+1), 1e+20)

        for iregion, region in enumerate(regions):

            nc_file = f'resampled_tfe.{region}.nc4'
            nc_path = os.path.join(src_directory, nc_file)
            tfes = Dataset(nc_path)['resampled_tfe'][:]  # (nmember, nK)

            if 'original' in process_type: 
                heatmap[iregion, 0] = np.median(tfes[:,0])
                tfes = tfes[:,0]
            elif 'bootstrap' in process_type: 
                heatmap[iregion, 0] = np.median(tfes)
                tfes = np.median(tfes, axis=1)

            for (igcm, gcm), (ighm, ghm) in itertools.product(enumerate(gcms), enumerate(ghms)):
                index = len(gcms)*ighm + igcm
                icol = len(ghms)*igcm + ighm
                heatmap[iregion, icol+1] = tfes[index]

        heatmap = np.ma.masked_equal(heatmap, 1e+20)
        draw_heatmap(heatmap, rcp, soc, co2, season)


if __name__=='__main__':
    main(*sys.argv)
