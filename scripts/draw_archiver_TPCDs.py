#! /usr/local/bin/python
import os, sys
import itertools
import datetime
import pandas as pd
import numpy  as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import mpl_toolkits.axisartist as axisartist
from netCDF4 import Dataset
from utiltools import get_region_info
plt.rcParams['font.family'] = 'Arial'
today = datetime.date.today().strftime('%Y%m%d')
hostname = os.uname()[1]


#ensemble_type = 'original'
ensemble_type = 'bootstrap'

experiment = 'basic'

sample_type = 'Mean'

trend_type = 'quadratic'
K = 100000

#trend_syear = 2005
#trend_syear = 1995
trend_syear = 1985
trend_syear = str(trend_syear)

blocksize = 5

#seasons = ['DRY', 'WET', 'annual']
seasons = ['DRY']


# OPTIONS ==================================================================================
# select consecutive chunk size for the TPCD analysis -------------------------------------
tChunks = [1, 4, 5, 6, 10, 15, 20, 99]
#tChunks = [1, 4, 5, 6, 10, 15, 20]
#tChunks = [4, 5, 6, 10, 15, 20, 99]
#tChunks = [4, 5, 6, 10, 15, 20, 25]
#tChunks = [3, 5, 10, 15, 20, 25, 30]
#tChunks = [5, 10, 15, 20, 25, 30]

# select stats value(s) to pick up from a kernel density function -------------------------
#smplTYPs = ['Average', 'Extreme']
#smplTYPs = ['Average']

# select a spatial scale you would like to estimate TPCD for ------------------------------
#region_type = 'AR6_regions'
region_type = 'HydroBASINS_lev1'
#region_type = 'BasinCountryUnit'
#region_type = 'Basin'
#region_type = 'Nation'

if region_type == 'AR6_regions':
    dict_input_date = {  # date of bootstrap
        #('Average', 'annual'): 20210701,
        #('Average', 'DRY'): 20210629,
        #('Average', 'WET'): 20210630,
        ('Mean', 'annual'): 20210707,
        ('Mean', 'DRY'   ): 20210707,
        ('Mean', 'WET'   ): 20210707,
        #('Mean', 'annual'): 20210726,  # with DryGridMask
        #('Mean', 'DRY'   ): 20210726,  # with DryGridMask
        #('Mean', 'WET'   ): 20210726,  # with DryGridMask
        }
elif region_type == 'HydroBASINS_lev1':
    dict_input_date = {  # date of bootstrap
        ('Mean', 'annual'): 20220103,
        ('Mean', 'DRY'   ): 20220103,
        ('Mean', 'WET'   ): 20220103,
        }
else:
    raise ValueError(f'check region_type: {region_type}')

# select how to deal with yealy sampling --------------------------------------------------
#kdfSampling = 'yearly'
#kdfSampling = 'window'; window = 2               # 5yr-sampling in total  (2+1+2)
#kdfSampling = 'window'; window = 3               # 7yr-sampling in total  (3+1+3)

# select how to sample from ensemble member for each grid. default is 'all' ---------------
#mmeTYPs = ['median','mean','all']
#mmeTYPs = ['all']
#mmeTYPs = ['median']
#mmeTYPs = ['mean']

# select threshold type -------------------------------------------------------------------
#autQvalue = True          # automatically specify, depending on soc type: 'Qvale_hist_%s'%soc
#autQvalue = False        # set qvaltype, below
#qvalTYPE = 'Qvale_hist_nosoc'
#qvalTYPE = 'Qvale_hist_pressoc'

# select which season you would like to analyze about --------------------------------------
#seasons = ['ALL','DJF','MAM','JJA','SON']
#seasons = ['ALL']

# select scenarios ------------------------------------------------------------------------
rcps = ['rcp26', 'rcp85']

# select socs. default is both for pressoc-nosoc comparison -------------------------------
#socs = ['2005soc', 'nosoc']
#socs = ['2005soc']
#socs = ['nosoc']
soc = '2005soc'
co2 = 'co2'

# drought indices --------------------------------------------------------------------------
#IDXs = ['nDayTot', 'dfcTot', 'nEvent', 'nOnset', 'avrDSL', 'maxDSL']
#IDXs = ['nDayTot']

# drought analysis parameters ---------------------------------------------------------------
# --- Qvalue
#Qs = [90, 85, 80, 75, 70, 65]
#Qs = [80]
# --- window size for Qvalue 
#wins = [15, 10, 7, 5]
#wins = [15]
# --- minimum drought days
#LENs = [180,90,60,30,14,7] 
#LENs = [30,14,60,90,180,7] 
#LENs = [30] 
# --- Pooled duration
#TAUs = [4]

suffixes = ['png', 'pdf']

# ----------------------------------------------------------------------------------------
data_directory = '/data/rg001/sgec0017/data'
bootstrap_directory = os.path.join(data_directory, 'figure_box', 'drought_tpcd_isimip2b', 'bootstrap')
mapmask_directory = os.path.join(data_directory, 'mapmask')

# a directory for figures
fig_directory_main = os.path.join(data_directory, 'figure_box', 'drought_tpcd_isimip2b', 'draw_archiver_TFEs')


# -----------------------------------------------------------------------------------------------------------
def load_tfe_data(tChunk, region, rcp, season):
    nc_directory = os.path.join(bootstrap_directory, region_type,
                                '{}.{}.{}{}99.{}.block{}.{}'.format(experiment, sample_type, trend_type, trend_syear[2:], K, blocksize, season), 
                                'tChunk_{:02}'.format(tChunk), 
                                str(dict_input_date[(sample_type,season)]), '{}_{}_{}'.format(rcp, soc, co2))
    nc_file = 'resampled_tfe.{}.nc4'.format(region)
    nc_path = os.path.join(nc_directory, nc_file)
    return Dataset(nc_path)['resampled_tfe'][:]  # (nmember, nsample)


def draw_clustered_bar_chart(df, rcp, season, fig_directory):

    def label_name(chunk):
        if chunk == 99:
            label_name = 'PE'
        else:
            label_name = chunk
        return label_name

    chunks = df.columns
    regions = df.index
    _regions = [f'{region[:2]}  {region[3:]}' for region in regions]
    print('chunks :', chunks)
    #print('regions:', regions)
    for _region in _regions: print(_region)

    ys     = [len(regions)*10 - (i_region*10+i_chunk) for i_region, region in enumerate(regions) for i_chunk, chunk in enumerate(chunks)]
    yticks = [len(regions)*10 - (i_region*10+3.5)       for i_region, region in enumerate(regions)]
    ylim_max = len(regions) * 10 + 3
    xticks = range(2010,2110,10)
    colors = [plt.get_cmap('jet_r')(i/len(chunks)) for i in range(len(chunks))]

    if region_type == 'AR6_regions':
        left = 0.12
        yaxis_tick_pad = 65
    elif region_type == 'HydroBASINS_lev1':
        left = 0.275
        yaxis_tick_pad = 203
    else:
        raise ValueError

    # maker setting, in accordance with the editor comment...
    markers =[
        'P',  # plus(filled)
        'X',  # x(filled)
        'o',  # circle 
        'v',  # triangle_down 
        '^',  # triangle_up
        'H',  # hexagon2 
        'D',  # Diamond 
        's',  # square
        ]
    if len(markers) < len(tChunks): 
        raise ValueError('len(markers) < len(tChunks). Add extra marker to markers.')

    # draw a bar chart
    fig = plt.figure(figsize=(12,15))
    plt.subplots_adjust(left=left, bottom=0.03, right=0.96, top=0.95)

    ax = fig.add_subplot(111)
    ax.xaxis.grid(True)
    ax.xaxis.set_ticks_position('top')
    ax.xaxis.set_label_position('top')
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    horizontal_line_color = '#868686'
    for i_region, region in enumerate(regions):
        xs = np.ma.masked_greater([df.loc[region, chunk] for chunk in chunks], 2099)
        for i_chunk, chunk in enumerate(chunks):
            ax.axhline(y=ys[i_region*len(chunks)+i_chunk], color=horizontal_line_color, linewidth=0.3, alpha=0.8)
        ax.plot(xs, ys[i_region*len(chunks):(i_region+1)*len(chunks)], linewidth=0.5, color='#626262') 
        # editted due to the edditor comment...
        #ax.scatter(xs, ys[i_region*len(chunks):(i_region+1)*len(chunks)], c=colors, s=40, edgecolor='k', lw=0.1)
        for i_chunk, chunk in enumerate(chunks):
            ax.scatter(xs[i_chunk], ys[i_region*len(chunks)+i_chunk], c=colors[i_chunk], marker=markers[i_chunk], s=40, edgecolor='k', lw=0.1)

    ax.set_yticks(yticks)
    ax.set_yticklabels(_regions, fontsize=15)
    ax.set_ylim([0.5, ylim_max])
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticks, fontsize=15)
    for label in ax.yaxis.get_ticklabels():
            label.set_horizontalalignment('left')
    ax.yaxis.set_tick_params(pad=yaxis_tick_pad)
    #ax.set_xlim([2006,2099])
    ax.set_xlim([2006,2110])

    legend_elements = [Line2D([0],[0], marker=markers[i_chunk], color='k', linewidth=0.05, markerfacecolor=colors[i_chunk], label=label_name(chunk), markersize=8) for i_chunk, chunk in enumerate(chunks)]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=13.5, title='minimum\nduration', title_fontsize=13.5)

    for suffix in suffixes:
        figname = 'barchart.{}.{}.{}.{}'.format(ensemble_type, season, rcp, suffix)
        figpath = os.path.join(fig_directory, figname)
        plt.savefig(figpath)
        print('savefig: {}'.format(figpath))
    plt.close()


# -----------------------------------------------------------------------------------------------------------
def main(*args):
    print('\n\nStart Job !!! d(^.^)\n')
    strTime = datetime.datetime.now()

    #=== preparations
    region_map, dict_region = get_region_info(region_type)
    del region_map
    regions = [f'{region_id:02}.{region_name}' for region_id, region_name in dict_region.items()]

    for season in seasons:
        
        info = '{}.{}{}.{}{}99.{}.{}.{}'.format(experiment, ensemble_type, K, trend_type, trend_syear[2:], blocksize, sample_type, season)
        fig_directory = os.path.join(fig_directory_main, region_type, info)  #, today)
        if not os.path.isdir(fig_directory): os.makedirs(fig_directory)
        print('fig_directory: {}'.format(fig_directory))

        output_excel_name = 'TFEs.archived.{}.{}.xlsx'.format(ensemble_type, season)
        outpath = os.path.join(fig_directory, output_excel_name)
        writer = pd.ExcelWriter(outpath)

        for rcp in rcps:
            print('\nrcp: {}'.format(rcp))

            df = pd.DataFrame(index=regions, columns=tChunks)
            for tChunk, region in itertools.product(tChunks, regions):
                print('tChunk {:02}  @{}'.format(tChunk, region))

                tfes = load_tfe_data(tChunk, region, rcp, season)

                if 'original' in ensemble_type:
                    tfe_median = np.median(tfes[:,0])
                elif 'bootstrap' in ensemble_type:
                    tfe_median = np.median(tfes)
                df.loc[region, tChunk] = tfe_median

            df = df.mask(df==9999.)
            print(df)
            df.to_excel(writer, sheet_name=rcp)
            draw_clustered_bar_chart(df, rcp, season, fig_directory)
        writer.save()
        print('write out: {}'.format(outpath))

    print('Successfully Done!!  d(^o^)b')

if __name__=='__main__':
    main(*sys.argv)
