#!/usr/bin/env python
# To plot time series of regional average FDD with statistical test
# By Yusuke Satoh
# On 20210628

import os
import sys
import itertools
import datetime
import pathlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from utiltools import get_region_info
plt.rcParams['font.family'] = 'Arial'

TEST = False

experiments = ['basic']
#experiments = ['picontrol']
#experiments = ['basic', 'picontrol']

#smplType  = 'Average'  # KDE
smplType  = 'Mean'  # area-weighted mean

#region_type = 'AR6_regions'
#region_type = 'HydroBASINS_lev1'
region_type = sys.argv[1]

seasons = [sys.argv[2]]
#seasons = ['annual']
#seasons = ['annual', 'DRY', 'WET']

ylim_types = [sys.argv[3]]
# ylim_types = ['full', 'auto']

gcms = ['hadgem2-es', 'ipsl-cm5a-lr', 'gfdl-esm2m', 'miroc5']  # full members
ghms = ['cwatm', 'h08', 'lpjml', 'matsiro', 'watergap2']  # full members


syear, eyear = 1865, 2099
years = range(syear, eyear+1)
nyear = len(years)
eyear_historical = 2005
years_historical = range(syear, eyear_historical+1)
years_future = range(eyear_historical, eyear+1)  # from 2005
index_hist_e = years.index(eyear_historical)

main_directory = '/data/rg001/sgec0017/data/figure_box/drought_tpcd_isimip2b'
figure_directory = os.path.join(main_directory, 'draw_plot_regional_avrFDD_timeseries')
mapmask_directory = '/data/rg001/sgec0017/data/mapmask'

if region_type == 'AR6_regions':
    regional_name_path = os.path.join(mapmask_directory, 'IPCC-AR6-WGI-reference-regions-v4_shapefile', 'regions.xlsx')
    df_region = pd.read_excel(regional_name_path, header=None)
    dict_regions = {
        1: 'NWN',   2: 'NEN',  3: 'WNA',   4: 'CNA',   5: 'ENA',
        6: 'NCA',   7: 'SCA',  8: 'CAR',   9: 'NWS',  10: 'NSA',
        11: 'NES',  12: 'SAM', 13: 'SWS',  14: 'SES',  15: 'SSA',
        16: 'NEU',  17: 'WCE', 18: 'EEU',  19: 'MED',  20: 'SAH',
        21: 'WAF',  22: 'CAF', 23: 'NEAF', 24: 'SEAF', 25: 'WSAF',
        26: 'ESAF', 27: 'MDG', 28: 'RAR',  29: 'WSB',  30: 'ESB',
        31: 'RFE',  32: 'WCA', 33: 'ECA',  34: 'TIB',  35: 'EAS',
        36: 'ARP',  37: 'SAS', 38: 'SEA',  39: 'NAU',  40: 'CAU',
        41: 'EAU',  42: 'SAU', 43: 'NZ',
        }
    dict_input_date = {  # date of TFE analysis
        #('basic', 'annual', 'Average'): 20210501,
        #('basic', 'DRY', 'Average'): 20210624,
        #('basic', 'WET', 'Average'): 20210630,
        ('basic', 'annual', 'Mean'): 20210722,
        ('basic', 'DRY', 'Mean'): 20210722,
        ('basic', 'WET', 'Mean'): 20210722,
        #('picontrol', 'annual', 'Average'): 20210505,
        #('picontrol', 'DRY', 'Average'): '',
        #('picontrol', 'WET', 'Average'): '',
        }
elif region_type == 'HydroBASINS_lev1':
    region_map, dict_regions = get_region_info(region_type, TEST=False)
    del region_map
    dict_input_date = {  # date of TFE analysis
        ('basic', 'annual', 'Mean'): 20220102,
        ('basic', 'DRY', 'Mean'):    20220102,
        ('basic', 'WET', 'Mean'):    20220102,
        }
else:
    raise ValueError(f'check region_type -> {region_type}')

drought_setting = 'Q80win15Len30tau4'
rangeType = 'IQR'

MANNWHITNEY = False
mannwhitney_sample = 'median'
#mannwhitney_sample = 'all'
#significance_level = 0.05
significance_level = 0.01


dict_scenarios = {
    'basic': ['rcp26', 'rcp85'],
    'picontrol': ['picontrol', 'rcp26', 'rcp85'],
}


dict_nday = {
    'annual': 365,
    'DRY': 91,
    'WET': 91,
    }

dict_color = {
    #'rcp26': '#0080ff',  # right blue
    #'rcp85': '#ff6600',  # orange  
    'rcp26': '#61aac6',  # right blue 2     
    'rcp85': '#d98a30',  #'#ff5500',  # orange  '#dc0d0d',
    'picontrol': '#808080'  #'#b0b0b0'  #
    }

dict_ylim = { 
    'full': [0,95], 
    'auto': [None,None],
    }

dict_yticks = {
    'full': range(10,95, 20),
    'auto': None,
    }


suffixes = ['png', 'pdf']

soc, co2  = '2005soc', 'co2'
index     = 'nDayTot'
mmeType   = 'all'
threshold = 'max'
dpi_figure, dpi_savefig = 300, 300

if TEST:
    dict_regions = {19: 'MED'}
    dpi_figure, dpi_savefig = 80, 80



def read_src_from_excel(experiment, season, scenario, region, ghm, gcm):
    #print(experiment, season, scenario, region, ghm, gcm)
    tfe_setting = f'{experiment}.abs.05.{smplType[:3]}.1865-2005.max.05.historical_histsoc_co2_1861-2005.1865-2099.plt1865-'
    excel_directory = os.path.join(main_directory, 'TPCD.estimator.scns.v2_for_1stRevision', region_type, 'all_all', 
                                   f'{drought_setting}{season}', tfe_setting, str(dict_input_date[(experiment, season, smplType)]))
    excel_file = f'TPCDs.nDayTot.005_{threshold}.{drought_setting}{season}.all.{smplType}.{scenario}_{soc}_{co2}.xlsx'
    excel_path = os.path.join(excel_directory, excel_file )
    df_input = pd.read_excel(excel_path, sheet_name='data', index_col=0)
    df_src = df_input.loc[f'{ghm}_{gcm}_{region}', syear:]
    df_src[df_src == '--'] = 1e+20 
    src = df_src.values
    src = np.ma.masked_equal((src.astype(np.float32)), 1e+20)
    return src


def mann_whitney_u_test(src_h, src_f):
    """
    To check if changes are significant.
    :param src_h: historical time series (full period)
    :param src_f: future time series
    :return: True or False
    """
    n_sample = 30
    weight = np.ones(n_sample)/n_sample
    if len(src_h.shape) == 1:
        index_max = np.argmax(np.convolve(src_f, weight, mode='valid')-src_h.mean())
        p_value = stats.mannwhitneyu(src_h, src_f[index_max: index_max+n_sample], alternative='two-sided')[1]
    elif len(src_h.shape) == 2:
        index_max = np.argmax(np.convolve(np.median(src_f, axis=0), weight, mode='valid')-src_h.mean())
        p_value = stats.mannwhitneyu(src_h.flatten(), src_f[:,index_max:index_max+n_sample].flatten(), alternative='two-sided')[1]
    significance = True if p_value < significance_level else False
    return significance, years_future[index_max]


def create_figureobject(region, ylim_type):
    #print('gen fig...')
    if not ylim_type == 'auto':
        figsize=(1.75, 1); left=0.12; bottom=0.03; right=0.98; top=0.95  # for Fig2
    else:  # ylim_type == 'auto'
        figsize=(1.5, 1); left=0.12; bottom=0.03; right=0.98; top=0.95  # for Fig2
    fig = plt.figure(num=1, figsize=figsize, dpi=dpi_figure)
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top)
    ax = fig.add_subplot(111)
    return fig, ax


def update_indevidual_timeseries(ax, median, srcs, scenario, season):
    #print('add plot...')

    lw = 0.5
    linestyle = '-'
    if co2 == '2005soc': linestyle = '--'

    if scenario == 'picontrol':
        color_hist, color_future = dict_color[scenario], dict_color[scenario]
        alpha_hist, alpha_future = 0.15, 0.15
    elif 'rcp' in scenario:
        color_hist, color_future = '#2b2b2b', dict_color[scenario]
        #alpha_hist, alpha_future  = 0.05, 0.25
        alpha_hist, alpha_future  = 0.1, 0.25
    else: raise ValueError

    median  = median / dict_nday[season] * 100  # ensemble timeseries (nYear)
    srcs = srcs / dict_nday[season] * 100  # each member (nGHM*nGCM,nYear)

    if   rangeType == 'minmax':
        aMinTS = np.min(srcs, axis=0)
        aMaxTS = np.max(srcs, axis=0)
    elif rangeType == 'IQR':
        aMinTS = np.percentile(srcs, 25., axis=0)
        aMaxTS = np.percentile(srcs, 75., axis=0)
    elif rangeType == '2sigma':
        aSigma = np.std(srcs, axis=0)
        aMinTS = median - 2*aSigma
        aMaxTS = median + 2*aSigma

    ax.fill_between(years_historical,
                    aMinTS[:index_hist_e+1],
                    aMaxTS[:index_hist_e+1],
                    facecolor=color_hist,
                    alpha=alpha_hist)
    ax.fill_between(years_future,
                    aMinTS[index_hist_e:],
                    aMaxTS[index_hist_e:],
                    facecolor=color_future,
                    alpha=alpha_future)
    ax.plot(years_historical,
            median[:index_hist_e+1],
            linestyle=linestyle, lw=lw, color=color_hist)
    ax.plot(years_future,
            median[index_hist_e:],
            linestyle=linestyle, lw=1.3*lw, color=color_future,
            label=scenario.upper())


def custum_indevidual_timeseries_axis(ax, region, ylim_type):
    #print('custum axis...')

    years_on_xaxis = [syear, 2006, 2050, eyear]
    spine_color = '#525252'
    lw = 0.45
    if not ylim_type == 'auto':  # full
        fontsize = 8
    else:  # ylim_type == 'auto'
        fontsize = 7

    # x-axis
    ax.set_xlim([syear, eyear])
    ax.set_xticks(years_on_xaxis)
    ax.set_xticklabels(['']*len(years_on_xaxis))
    ax.tick_params(axis='x', direction='in', length=3.5, pad=1, labelsize=fontsize)
    # y-axis 
    if not ylim_type == 'auto':
        ax.set_ylim(dict_ylim[ylim_type] )
        ax.set_yticks(dict_yticks[ylim_type])
        ax.set_yticklabels(dict_yticks[ylim_type])
    else:  # ylim_type == 'auto'
        ax.yaxis.get_major_locator().set_params(integer=True)
    ax.tick_params(axis='y', direction='in', length=3, pad=1, labelsize=fontsize)
    # spines
    ax.patch.set_alpha(0.7)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_color(spine_color)
    ax.spines['left'].set_linewidth(lw)

    _y = 1; va = 'top'
    _x = 0.06; ha = 'left'
    #region_id, region_abbrev = region.split('.')
    #region_id = int(region_id)
    region_id = int(region[:2])
    region_abbrev = region[3:]
    if region_type == 'AR6_regions':
        region_full_name = df_region[df_region[3]==region_abbrev][2].values[0]
        if region_id == 12: region_full_name = 'South-American\n   -Monsoon'
        ax.text(_x, _y, '{}.{}'.format(region_id, region_full_name), ha=ha, va=va, fontsize=fontsize, transform=ax.transAxes)
    elif region_type == 'HydroBASINS_lev1':
        _region_name = dict_regions[region_id]
        if region_id == 32: _region_name = _region_name.replace('.', '\n     ')
        ax.text(_x, _y, '{} {}'.format(region_id, _region_name), ha=ha, va=va, fontsize=fontsize, transform=ax.transAxes)
    #ax.text(0, 1.05, '[%]', ha='right', va='bottom', fontsize=fontsize, transform=ax.transAxes)
    #if regionID == 3:
    #    ax.legend(prop={'size' : 9},
    #              bbox_to_anchor=(0, 1.08), loc='upper left',
    #              labelspacing=0.2,
    #              frameon=False)


def close_figureobject(region, season, ylim_type, fig_directory):

    figname = f'ts.abs.{rangeType}.{index}.{drought_setting}{season}.Emsenble_all.{smplType}.{region}.'
    plt.figure(num=1)
    for suffix in suffixes:
        _fig_directory = os.path.join(fig_directory, smplType, ylim_type, suffix)
        if not os.path.isdir(_fig_directory): os.makedirs(_fig_directory)
        figpath = os.path.join(_fig_directory, figname+suffix)
        plt.savefig(figpath, dpi=dpi_savefig)
        if suffix == suffixes[0]:
            print(f'savefig: {figpath}')
    plt.close()
    print('\n')



# ----------------------------------------------------------------------------------------------------------------------
def main():
    stime = datetime.datetime.now()

    regions = [f'{region_id:02}.{region_name}' for region_id, region_name in dict_regions.items()]

    for experiment, season, ylim_type in itertools.product(experiments, seasons, ylim_types):
        if MANNWHITNEY:
            df_mannwhitney = pd.DataFrame(index=regions,
                                          columns=[f'{item}.{scenario}'
                                                   for item in ['test', 'syear_of_the30yrs']
                                                   for scenario in dict_scenarios[experiment]])
        for region in regions:

            fig_directory = os.path.join(figure_directory, region_type, experiment, season)
            if not os.path.isdir(fig_directory): os.makedirs(fig_directory)
            tfe_setting = f'{experiment}.abs.05.{smplType:3}.1865-2005.max.05.historical_histsoc_co2_1861-2005.1865-2099.plt1865-'
            pathlib.Path(os.path.join(fig_directory, 'input_from_'+tfe_setting)).touch()  # crease a empty file just for a memo

            # generate figure and axis
            fig, ax = create_figureobject(region, ylim_type)

            for scenario in dict_scenarios[experiment]:
                print(f'{scenario}...')

                #srcs = np.array([read_src_from_excel(experiment, season, scenario, region, ghm, gcm) for ghm in ghms for gcm in gcms])  # (nmember, nyear)

                tfe_setting = f'{experiment}.abs.05.{smplType[:3]}.1865-2005.max.05.historical_histsoc_co2_1861-2005.1865-2099.plt1865-'
                excel_directory = os.path.join(main_directory, 'TPCD.estimator.scns.v2_for_1stRevision', region_type, 'all_all', 
                                               f'{drought_setting}{season}', tfe_setting, str(dict_input_date[(experiment, season, smplType)]))
                excel_file = f'TPCDs.nDayTot.005_{threshold}.{drought_setting}{season}.all.{smplType}.{scenario}_{soc}_{co2}.xlsx'
                excel_path = os.path.join(excel_directory, excel_file )
                df_input = pd.read_excel(excel_path, sheet_name='data', index_col=0)
                #df_input[df_input == '--'] = 1e+20 
                df_input = df_input.mask(df_input == '--', 1e+20) 

                srcs = np.array([df_input.loc[f'{ghm}_{gcm}_{region}', syear:].values.astype(np.float32) for ghm in ghms for gcm in gcms])  # (nmember, nyear)
                srcs = np.ma.masked_equal(srcs, 1e+20)
                median = np.median(srcs, axis=0)

                if MANNWHITNEY:
                    print('Mann Whitney U-test...')
                    if mannwhitney_sample == 'median':
                        significance, syear_of_the30yrs = mann_whitney_u_test(median[:index_hist_e+1], median[index_hist_e+1:])  # True or False
                    elif mannwhitney_sample == 'all':
                        significance, syear_of_the30yrs = mann_whitney_u_test(srcs[:,:index_hist_e+1], srcs[:,index_hist_e+1:])  # True or False
                    df_mannwhitney.loc[region, f'test.{scenario}'] = significance
                    df_mannwhitney.loc[region, f'syear_of_the30yrs.{scenario}'] = syear_of_the30yrs
                    print(f'{experiment} {season} {region} {scenario} ---> {significance}  {syear_of_the3}~{syear_of_the30yrs+29}')

                # add a plot to the axis
                update_indevidual_timeseries(ax, median, srcs, scenario, season)

            # custum
            custum_indevidual_timeseries_axis(ax, region, ylim_type)

            # save
            close_figureobject(region, season, ylim_type, fig_directory)

        if MANNWHITNEY:
            mannwhitney_excel_file = f'mann_whitney_u_test.{mannwhitney_sample}.{significance_level}.xlsx'
            mannwhitney_excel_path = os.path.join(fig_directory, mannwhitney_excel_file)
            df_mannwhitney.to_excel(mannwhitney_excel_path)
            print(f'save: {mannwhitney_excel_path}')

    print('the process took {} minutes in total...'.format((datetime.datetime.now() - stime).seconds / 60))
    print('Successfully DONE!! d(^o^)b')


if __name__ == '__main__':
    main()
