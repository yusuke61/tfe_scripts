#!/usr/bin/env python
# To draw bar plots regarding total UD years for each regions
# By Yusuke Satoh
# On 2020.May.12

import os
import sys
import itertools
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from tqdm import tqdm
plt.rcParams['font.family'] = 'Arial'

# --- setting

#region_type = 'AR6_regions'
region_type = 'HydroBASINS_lev1'

experiment = 'basic'

sample_type = 'Mean'
#sample_type = 'Average'

trend_type = 'quadratic'
#trend_syear = 2005
#trend_syear = 1995
trend_syear = 1985
K = 100000
blocksize = 5
tChunk = 5

#median_type = 'original_member'
median_type = 'bootstrap_full_member'

#uncertainty_range = 'bootstrap_25-75_range'
uncertainty_range = 'bootstrap_05-95_range'
#uncertainty_range = 'quantile'
#uncertainty_range = '1sigma'
#uncertainty_range = 'minmax'
#uncertainty_range = 'comblined'

#seasons = ['annual', 'DRY', 'WET']
#seasons = ['DRY', 'WET', 'annual']
seasons = ['DRY']


# --- fixed parameters
file_directory = f'/data/rg001/sgec0017/data/figure_box/drought_tpcd_isimip2b/bootstrap/{region_type}'
figure_directory_main = f'/data/rg001/sgec0017/data/figure_box/drought_tpcd_isimip2b/draw_accumulated_UDyears/{region_type}'
if not os.path.isdir(figure_directory_main): os.makedirs(figure_directory_main)

rcps = ['rcp26', 'rcp85']
dict_rcp = {
    #'rcp26': '#0080ff',  # right blue    
    #'rcp85': '#ff6600',  # orange       
    'rcp26': '#61aac6',  # right blue 2  
    'rcp85': '#d98a30',  # brown       
}

ensemble_stats = ['median', 'BS05', 'BS95', 'BS25', 'BS75']
#ensemble_stats = ['median', 'max', 'min', '1sigma','75percentile', '25percentile', 'BS05', 'BS95']

if region_type == 'AR6_regions':
    dict_input_date = {  # date of bootstrap
        #('Average','annual'): 20210701,
        #('Average','DRY'): 20210625,
        #('Average','WET'): 20210630,
        ('Mean','annual'): 20210707,
        ('Mean','DRY'): 20210707,
        ('Mean','WET'): 20210707,
        }
elif region_type == 'HydroBASINS_lev1':
    dict_input_date = {  # date of bootstrap
        ('Mean','annual'): 20220103,
        ('Mean','DRY'):    20220103,
        ('Mean','WET'):    20220103,
        }

fontsize_large = 5
fontsize_small = 4
line_color_1 = '#353535'
line_color_2 = '#353535'
markeredgecolor = '#bfbfbf'
markeredgewidth = 0.2

#suffixes = ['png', 'eps', 'pdf']
suffixes = ['png', 'pdf']


# ----------------------------------------------------------------------------------------------------------------------
def get_regions_with_robust_median(season):
    
    bootstrap_tfe_summary_excel_path = os.path.join(
        '/data/rg001/sgec0017/data/figure_box/drought_tpcd_isimip2b/draw_tfe_map_from_bootstrap', region_type,
        f'{experiment}.bootstrap{K}.{trend_type}{str(trend_syear)[2:]}99.{blocksize}.{sample_type}.{season}',
        f'tfe_summary.{experiment}.bootstrap{K}.{trend_type}{str(trend_syear)[2:]}99.{blocksize}.{sample_type}.{season}.median.{tChunk}.xlsx'
        ) 
    df = pd.read_excel(bootstrap_tfe_summary_excel_path, sheet_name='medianTFE(bootstrap)', index_col=0)
    regions = list(df[(df['medi-rcp26']==1)|(df['medi-rcp85']==1)].index)  # if both or eigher of rcps are robust
    print(f'regions ({len(regions)}): {regions}')

    return regions


# ----------------------------------------------------------------------------------------------------------------------
def load_tfes_nc(region, rcp, season, LOUD=False):
    nc_file = f'n_recordbreaking_years.{region}.nc4'
    nc_path = os.path.join(file_directory, 
                           f'{experiment}.{sample_type}.{trend_type}{str(trend_syear)[2:]}99.{K}.block{blocksize}.{season}', 
                           f'tChunk_{tChunk:02}',
                           str(dict_input_date[(sample_type,season)]), f'{rcp}_2005soc_co2',
                           nc_file)
    if LOUD: print(f'loading... {nc_path}')
    src = Dataset(nc_path)['n_recordbreaking_years'][:]
    if LOUD: print(f' >>> {region} {src.min()}-{src.max()}  ({src.shape})\n')
    return src  # n_member, n_K


# ----------------------------------------------------------------------------------------------------------------------
def main():

    for season in seasons:
        
        # preparation
        regions = get_regions_with_robust_median(season)  # ex., ['08.North_America_2', ...]
        print(f'regions with robust change:\n{regions}\n')

        dict_src = {ensemble_stat: pd.DataFrame(index=regions, columns=rcps) for ensemble_stat in ensemble_stats}
        for rcp in rcps:
            print(f'\n{rcp}')
            # get ensemble stats
            for region in regions:
                Ns = load_tfes_nc(region, rcp, season, LOUD=True)  # (n_member x n_K)
                for ensemble_stat in ensemble_stats:
                    if ensemble_stat == 'median': 
                        if median_type == 'original_member':
                            dict_src['median'].loc[region, rcp] = np.median(Ns[:,0])  # median of original 20 member
                        elif median_type == 'bootstrap_full_member':
                            dict_src['median'].loc[region, rcp] = np.median(Ns)  # median of full bootstrap member
                    elif ensemble_stat == 'max': dict_src['max'].loc[region, rcp] = np.max(Ns)
                    elif ensemble_stat == 'min': dict_src['min'].loc[region, rcp] = np.min(Ns)
                    elif ensemble_stat == '1sigma': dict_src['1sigma'].loc[region, rcp] = np.std(Ns)
                    elif ensemble_stat == '75percentile': dict_src['75percentile'].loc[region, rcp] = np.quantile(Ns, 0.75)
                    elif ensemble_stat == '25percentile': dict_src['25percentile'].loc[region, rcp] = np.quantile(Ns, 0.25)
                del Ns

        # test: if rcps are significantly different.
        SEED = 1
        np.random.seed(SEED)
        rng = np.random.default_rng()
        regions_significantly_different = []
        print('\n\n--- start bootstrap test ---')
        for region in regions:
            # load src
            rcp26 = load_tfes_nc(region, 'rcp26', season)  # (20, 100001)
            rcp85 = load_tfes_nc(region, 'rcp85', season)  # (20, 100001)
            rcp26 = rcp26[:,:-1].flatten()         # (20, 100000) --> (2000000)
            rcp85 = rcp85[:,:-1].flatten()         # (20, 100000) --> (2000000)
            # random shuffle srcs
            rng.shuffle(rcp26)
            rng.shuffle(rcp85)
            # split src into 1000 and calculate medians
            # reshape ---> (1000 group, 2000 samples)
            rcp26 = np.median(rcp26.reshape(2000,1000), axis=1)  # (2000)
            rcp85 = np.median(rcp85.reshape(2000,1000), axis=1)  # (2000)
            # test their differences
            diff = rcp85 - rcp26  # (2000)
            tail_left, tail_right = np.quantile(diff, 0.05), np.quantile(diff, 0.95)
            if 0 < tail_left or tail_right < 0:
                print(f'{region}: o  Significant. ;)   tail_left={tail_left}, tail_right={tail_right}')
                regions_significantly_different.append(region)
            else:
                print(f'{region}: x  NOT significant.  tail_left={tail_left}, tail_right={tail_right}')
            if 'BS05' in dict_src.keys():
                dict_src['BS05'].loc[region, 'rcp26'] = np.quantile(rcp26, 0.05)
                dict_src['BS05'].loc[region, 'rcp85'] = np.quantile(rcp85, 0.05)
            if 'BS25' in dict_src.keys():
                dict_src['BS25'].loc[region, 'rcp26'] = np.quantile(rcp26, 0.25)
                dict_src['BS25'].loc[region, 'rcp85'] = np.quantile(rcp85, 0.25)
            if 'BS75' in dict_src.keys():
                dict_src['BS75'].loc[region, 'rcp26'] = np.quantile(rcp26, 0.75)
                dict_src['BS75'].loc[region, 'rcp85'] = np.quantile(rcp85, 0.75)
            if 'BS95' in dict_src.keys():
                dict_src['BS95'].loc[region, 'rcp26'] = np.quantile(rcp26, 0.95)
                dict_src['BS95'].loc[region, 'rcp85'] = np.quantile(rcp85, 0.95)
        print(f'---\nN of significant: {len(regions_significantly_different)}\n')

        # figure
        bar_width = 0.4
        xticks = np.arange(1, len(regions)+1)
        xticks_rcp26 = xticks - bar_width/2
        xticks_rcp85 = xticks + bar_width/2
        xticklabels = [f'{region[:2]}* {region[3:]}' if region in regions_significantly_different else f'{region[:2]}  {region[3:]}' for region in regions]
        yticklabels = range(0, 70 + 1, 10)
        yticks = yticklabels
        if region_type == 'AR6_regions':
            figsize = (3.3, 1.8)
            bottom = 0.18
        elif region_type == 'HydroBASINS_lev1':
            #figsize = (3.3, 2.3); left = 0.095
            figsize = (4.8, 2.3); left = 0.07
            bottom = 0.34
        else:
            raise ValueError

        fig = plt.figure(figsize=figsize, dpi=300)
        plt.subplots_adjust(left=left, right=0.99, bottom=bottom, top=0.99)
        ax = fig.add_subplot(111)
        ax.spines['top'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_linewidth(0.25)
        ax.set_xlim([xticks_rcp26[0]-0.3, xticks_rcp85[-1]+0.3])
        ax.set_xlabel('Regions', fontsize=fontsize_large)
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticklabels, fontsize=fontsize_small, rotation=-80)
        ax.set_ylabel('The total number of years\nunder unprecedented drought conditions', fontsize=fontsize_large)
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticklabels, fontsize=fontsize_small)
        ax.tick_params(width=0.3, length=1.5, pad=0.5)
        ax.grid(axis='y', color='#dbdbdb', linewidth=0.5, zorder=0)

        for _xticks, rcp in zip([xticks_rcp26, xticks_rcp85], rcps):
            print(f'\n{rcp}')
            median = dict_src['median'][rcp].values.tolist()
            print(f'median: {median} ({len(median)})')
            if rcp == 'rcp26': RCP = 'RCP2.6'
            elif rcp == 'rcp85': RCP = 'RCP8.5'
            ax.bar(_xticks, median,
                   width=bar_width, bottom=0, align='center', color=dict_rcp[rcp], alpha=0.6, edgecolor=dict_rcp[rcp], linewidth=0, zorder=3, label=RCP)
            # uncertainty information
            if uncertainty_range == '1sigma':
                sigma = dict_src['1sigma'].loc[region, rcp]; print(f'sigma: {sigma}')
                ax.plot([_xticks, _xticks], [median-sigma, median+sigma], lw=0.5, color='#d0d0d0', zorder=4)
            elif uncertainty_range == 'quantile':
                percentile_25 = dict_src['25percentile'].loc[:, rcp]; print(f'25percentile: {percentile_25.values.tolist()}')
                percentile_75 = dict_src['75percentile'].loc[:, rcp]; print(f'75percentile: {percentile_75.values.tolist()}')
                ax.plot([_xticks, _xticks], [percentile_25, percentile_75], 
                        lw=0.03, color=line_color_1, marker='_', markersize=1, markeredgecolor=line_color_1, markeredgewidth=markeredgewidth, zorder=4)
            elif uncertainty_range == 'minmax':
                _min = dict_src['min'].loc[:, rcp]; print(f'_min: {_min.values.tolist()}')
                _max = dict_src['max'].loc[:, rcp]; print(f'_max: {_max.values.tolist()}')
                ax.plot([_xticks, _xticks], [_min, _max], lw=0.5, color=line_color, zorder=4)
            elif uncertainty_range == 'comblined':
                percentile_25 = dict_src['25percentile'].loc[:, rcp]; print(f'25percentile: {percentile_25.values.tolist()}')
                _max = dict_src['max'].loc[:, rcp]; print(f'_max: {_max.values.tolist()}')
                ax.plot([_xticks, _xticks], [percentile_25, _max], lw=0.5, color=line_color, marker='_', markersize=3, markeredgecolor=markeredgecolor, zorder=4)
            elif uncertainty_range == 'bootstrap_25-75_range':
                percentile_25 = dict_src['BS25'].loc[:, rcp]; print(f'25percentile: {percentile_25.values.tolist()}')
                percentile_75 = dict_src['BS75'].loc[:, rcp]; print(f'75percentile: {percentile_75.values.tolist()}')
                ax.plot([_xticks, _xticks], [percentile_25, percentile_75], 
                        lw=0.05, color=line_color_1, marker='_', markersize=1.3, markeredgecolor=line_color_1, markeredgewidth=markeredgewidth, zorder=4)
            elif uncertainty_range == 'bootstrap_05-95_range':
                percentile_05 = dict_src['BS05'].loc[:, rcp]; print(f' 5percentile: {percentile_05.values.tolist()}')
                percentile_95 = dict_src['BS95'].loc[:, rcp]; print(f'95percentile: {percentile_95.values.tolist()}')
                ax.plot([_xticks, _xticks], [percentile_05, percentile_95], 
                        lw=0.1, color=line_color_1, marker='_', markersize=2.5, markeredgecolor=line_color_1, markeredgewidth=markeredgewidth, zorder=4)

        ax.legend(fontsize=fontsize_small, labelspacing=0.5, frameon=False)
        for suffix in suffixes:
            figure_name = f'totalnumber_UDyears.{season}.{uncertainty_range}.{tChunk:02}.{suffix}'
            figure_dir = os.path.join(figure_directory_main, 
                                      f'{experiment}.{sample_type}.{trend_type}{str(trend_syear)[2:]}99.{K}.{blocksize}.{season}')
            if not os.path.isdir(figure_dir): os.makedirs(figure_dir)
            figure_path = os.path.join(figure_dir, figure_name)
            plt.savefig(figure_path, dpi=300)
            print(f'savefig: {figure_path}')
        plt.close()

    print('Successfully DONE!!  d(^o^)b')

if __name__ == '__main__':
    main()
