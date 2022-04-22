#!/usr/bin/env python
# To draw CDF plot from bootstrap
# By Yusuke Satoh
# On 2021/06/24
import os
import sys
import math
import datetime
import itertools
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from utiltools import interpolated_intercepts, get_region_info
plt.rcParams['font.family'] = 'Arial'
today = datetime.date.today().strftime('%Y%m%d')

TEST = False

#ensemble_type = 'original'
ensemble_type = 'bootstrap'

experiment = 'basic'

#sample_type = 'Average'
sample_type = 'Mean'

trend_type = 'quadratic'

blocksize = 5

K = 100000

#region_type = 'AR6_regions'
#region_type = 'HydroBASINS_lev1'
region_type = sys.argv[1]

trend_syear = int(sys.argv[2])
#trend_syear = 2005

figure_type = sys.argv[3]
#figure_type = 'forMain'
#figure_type = 'ExtendedFig5'

seasons = [sys.argv[4]]
#seasons = ['annual', 'DRY', 'WET']
#seasons = ['annual']

target_chunks = [sys.argv[5]]
#target_chunks = [5, 10]

ghms = ['cwatm', 'h08', 'lpjml', 'matsiro', 'watergap2']
gcms = ['hadgem2-es', 'ipsl-cm5a-lr', 'gfdl-esm2m', 'miroc5']

#rcps = ['rcp26', 'rcp85']
#rcps = ['rcp26']
rcps = ['rcp85']

if TEST:
    ghms = ['cwatm']
    gcms = ['hadgem2-es']
    #rcps = ['rcp85']


dict_color = {
    #'rcp26': '#0080ff',  # right blue
    #'rcp85': '#ff6600',  # orange
    'rcp26': '#61aac6',  # right blue 2     
    'rcp85': '#d98a30',  #'#ff5500',  # orange  '#dc0d0d',
    }

# x-range
if figure_type == 'forMain':        low_dgmt, high_dgmt = 1., 5.0
elif figure_type == 'ExtendedFig5': low_dgmt, high_dgmt = 1., 6.0
else: raise ValueError
dict_dgmts = {
    'rcp85': np.arange(low_dgmt, high_dgmt+0.1, 0.5).tolist(),
    'rcp26': np.arange(low_dgmt, 2.6+0.1, 0.5).tolist(),
    }
counter_thresholds = [1.5, 2.0]
missing_value = 9999

main_directory = '/data/rg001/sgec0017/data/figure_box/drought_tpcd_isimip2b'
fig_directory = os.path.join(main_directory, 'draw_dgmt_CDF_from_bootstrap', region_type, sample_type)
if not os.path.isdir(fig_directory): os.makedirs(fig_directory)
srcdirectory_tdgmt =  os.path.join(main_directory, 'extract_temperature')
mapmask_directory = '/data/rg001/sgec0017/data/mapmask'

if region_type == 'AR6_regions':
    regional_name_path = os.path.join(mapmask_directory, 'IPCC-AR6-WGI-reference-regions-v4_shapefile', 'regions.xlsx')
    df_region = pd.read_excel(regional_name_path, header=None)
    dict_input_date = {  # date of boostrap (basic)
        ('Mean', 'annual'): 20210707,
        ('Mean', 'DRY'   ): 20210707,
        ('Mean', 'WET'   ): 20210707,
        }
elif region_type == 'HydroBASINS_lev1':
    dict_input_date = {  # date of boostrap (basic)
        ('Mean', 'annual'): 20220103,
        ('Mean', 'DRY'   ): 20220103,
        ('Mean', 'WET'   ): 20220103,
        }
    regional_name_path = os.path.join(mapmask_directory, 'HydroSHEDS', 'HydroBASINS', 
                                      'withLakes', 'Global_merge', 
                                      'hybas_lake____lev02_v1c_merge.xlsx'
                                      )
    df_region_longname = pd.read_excel(regional_name_path, header=0, usecols=[18, 19], index_col=0)  # {region_ID: long_name}
else:
    raise ValueError(f'chck region_type  {region_type}')

SEED = 1
np.random.seed(SEED)
rng = np.random.default_rng()


def load_tfes_nc(region, rcp, season, target_chunk):
    ncfile = f'resampled_tfe.{region}.nc4'
    input_directory = os.path.join(
        main_directory, 'bootstrap', region_type,
        f'{experiment}.{sample_type}.{trend_type}{str(trend_syear)[2:]}99.{K}.block{blocksize}.{season}',
        f'tChunk_{target_chunk:02}', str(dict_input_date[(sample_type,season)]), f'{rcp}_2005soc_co2'
        )
    nc_path = os.path.join(input_directory, ncfile)
    print(f'loading... {nc_path}')
    return Dataset(nc_path)['resampled_tfe'][:,:-1]  # (n_member, n_samples)


def conv_tfe2dgmt(tfes, rcp):  # tfes: (nmember, K)

    # load TdGMT excel file
    tdgmt_excelname = 'isimip2b_temperature_timeseries__baseperiod_1850-1900.xlsx'
    tdgmt_excel_path = os.path.join(srcdirectory_tdgmt, tdgmt_excelname)
    df_dgmt_org = pd.read_excel(tdgmt_excel_path, sheet_name='dGMT', index_col=0)

    _dgmts = []
    for (ighm, ghm), (igcm, gcm), k in itertools.product(enumerate(ghms), enumerate(gcms), range(K)):
        column = f'{rcp}_{gcm}'
        tfe = int(tfes[len(gcms) * ighm + igcm, k])
        if tfe <= 2084:
            dgmt = df_dgmt_org.loc[str(tfe), column]
            _dgmts.append(dgmt)
        elif 2085 <= tfe <= 2099:
            dgmt = df_dgmt_org.loc['2084', column]
            _dgmts.append(dgmt)
        else:
            dgmt = missing_value
            _dgmts.append(dgmt)

    return np.array(_dgmts)


def estimate_frequency(input_tfes, rcp):
    total = input_tfes.shape
    cumulative_density = []
    dgmts = dict_dgmts[rcp]
    for dgmt in dgmts:
        cumulative_density.append((np.where(input_tfes<=dgmt,1,0).sum()/total)[0])
    return cumulative_density


def bootstrap_medians(src):
    rng.shuffle(src)  # shuffle
    medians = np.median(src.reshape(2000, 1000), axis=1)  # (2000)
    return medians


def generate_container(region_id, region_abbrev, region_full_name):

    # fixed parameters
    figsize = (5, 6)
    left, bottom, right, top = 0.15, 0.052, 0.96, 0.945
    fontsize_region_name = 19
    fontsize_ticks = 17
    xlim_range = [low_dgmt-0.1, high_dgmt+0.1]

    fig = plt.figure(figsize=figsize)
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top)
    ax = fig.add_subplot(111)
    ax.set_title(f'{region_id}. {region_full_name}', fontsize=fontsize_region_name)
    for position in ['right', 'left']:
        ax.spines[position].set_visible(False)
    # y-axis
    yticks = [0, 0.33, 0.5, 0.75, 0.95, 1]
    yticklabels = yticks
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels, fontsize=fontsize_ticks)
    ax.set_ylim([0, 1])
    # x-axis
    xticks = np.arange(low_dgmt, high_dgmt+0.1, 0.5)
    xticklabels = xticks
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels, fontsize=fontsize_ticks)
    ax.set_xlim(xlim_range)
    # horizontal lines at 0.05, 0.5, and 0.95
    #ax.axhline(y=0.05, color='#808080', linewidth=0.8, linestyle='-')  # gray
    ax.axhline(y=0.33, color='#808080', linewidth=0.8, linestyle='-')  # gray
    ax.axhline(y=0.50, color='#808080', linewidth=0.8, linestyle='-')  # gray
    ax.axhline(y=0.66, color='#808080', linewidth=0.8, linestyle='-')  # gray
    ax.axhline(y=0.95, color='#808080', linewidth=0.8, linestyle='-')  # gray
    # virtical lines at 2050
    ax.axvline(x=1.5, color='#808080', linewidth=0.8, linestyle='-')  # gray
    ax.axvline(x=2.0, color='#808080', linewidth=0.8, linestyle='-')  # gray
    # grid
    ax.grid(color='#e0e0e0', linewidth=0.6)
    # ---
    ax_pos = ax.get_position()
    ax2 = fig.add_axes([ax_pos.x0, ax_pos.y0-0.2, ax_pos.x1-ax_pos.x0, 0.2])
    ax2.axis('off')
    ax2.set_ylim([-0.2, 0])
    ax2.set_xlim(xlim_range)

    return fig, ax, ax2, ax_pos


def add_uncertainty_plot(ax, rcp, bs_cdf_max, bs_cdf_min):

    # plot CDF uncertainty
    xs = dict_dgmts[rcp]
    ax.fill_between(xs, bs_cdf_max, bs_cdf_min, color=dict_color[rcp], alpha=0.1)
    ax.plot(xs, bs_cdf_max, color=dict_color[rcp], linestyle='-', linewidth=0.5, markersize=0)
    ax.plot(xs, bs_cdf_min, color=dict_color[rcp], linestyle='-', linewidth=0.5, markersize=0)
    ## median error bars
    #if rcp == 'rcp85': y_median_errorbar = -0.06
    #elif rcp == 'rcp26': y_median_errorbar = -0.1
    #ax2.plot([robust_median_05, robust_median_95], [y_median_errorbar, y_median_errorbar],
    #         color=dict_color[rcp], linewidth=1.5, linestyle='--', marker='|')
    ## overall error bars
    #if rcp == 'rcp85': y_overall_errorbar = -0.08
    #elif rcp == 'rcp26': y_overall_errorbar = -0.12
    #ax2.plot([robust_overall_05, robust_overall_95], [y_overall_errorbar, y_overall_errorbar],
    #         color=dict_color[rcp], linewidth=1.5, linestyle='-', marker='|')

def add_main_plot(ax, rcp, cdf):
    xs = dict_dgmts[rcp]
    ax.plot(xs, cdf, color=dict_color[rcp], linewidth=2.5, markersize=0)

def add_table(fig, ax_pos, cell_text):
    # table
    ax3 = fig.add_axes([ax_pos.x0+0.503, ax_pos.y1-0.75, 0.3, 0.15])
    ax3.axis('off')
    ax3.set_title('Likelihood')
    columns = [u'by1.5\u00B0C', u'by2\u00B0C']
    #cell_text.reverse()
    table = ax3.table(cellText=cell_text,
                      rowLabels=rcps, colLabels=columns, 
                      cellColours=None, loc='lower right',
                      edges='open',
                      cellLoc='center')
    table.auto_set_font_size(False)
    table.set_fontsize(15)
    table.scale(1,2)
    ax3.axhline(y=1, color='#808080') 

def savefig(target_chunk, region, season):
    for suffix in ['png', 'pdf']:
        figname = f'dgmt_cdf_indevidual_split.chunk{target_chunk:03}.{season}.{region}.{suffix}'
        figdir = os.path.join(
                    fig_directory, 
                    f"{experiment}.{sample_type}.{trend_type}{str(trend_syear)[2:]}99.{K}.{blocksize}.{season}.{''.join([rcp[-2:] for rcp in rcps])}",
                    figure_type, suffix
                    )
        if not os.path.isdir(figdir): os.makedirs(figdir)
        figpath = os.path.join(figdir, figname)
        plt.savefig(figpath)
        print(f'savefig: {figpath}')
    plt.close()
    
    
def gen_writer_df(topic, target_chunk, season, regions):
    mmm = ['median', 'min', 'max']
    excel_file = f'{topic}_summary.{ensemble_type}.{target_chunk}.{season}.xlsx'
    figdir = os.path.join(
                fig_directory, 
                f"{experiment}.{sample_type}.{trend_type}{str(trend_syear)[2:]}99.{K}.{blocksize}.{season}.{''.join([rcp[-2:] for rcp in rcps])}.{figure_type}",
                )
    if not os.path.isdir(figdir): os.makedirs(figdir)
    excel_path = os.path.join(figdir, excel_file)
    writer = pd.ExcelWriter(excel_path)
    if topic == 'median_tfe':
        df = pd.DataFrame(index=[region for region in regions], columns=rcps)
        return writer, df
    elif topic == 'robust_level':
        columns2 = [f'by{counter_threshold}.{rcp}.{m}' for rcp in rcps
                                                       for counter_threshold in counter_thresholds
                                                       for m in mmm]
        colmuns3 = [f'{rcp}.{m}' for rcp in rcps for m in mmm]
        df2 = pd.DataFrame(index=regions, columns=columns2)
        df3 = pd.DataFrame(index=regions, columns=colmuns3)
        return writer, df2, df3

def save_writer(sheet_name, writer, df):
    df.to_excel(writer, sheet_name=sheet_name)
    writer.save()
    writer.close()


# ----------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------
def main():

    region_map, dict_region = get_region_info(region_type, TEST=TEST)
    regions = [f'{region_id:02}.{region_name}' for region_id, region_name in dict_region.items()]

    for target_chunk, season in itertools.product(target_chunks, seasons): 
        target_chunk = int(target_chunk)
        writer1, df1 = gen_writer_df('median_tfe', target_chunk, season, regions)
        writer2, df2, df3 = gen_writer_df('robust_level', target_chunk, season, regions)

        for region in regions:
            #region_id, region_abbrev = region.split('.')
            #region_id = int(region_id)
            region_id = int(region[:2])
            region_abbrev = region[3:]
            if region_type == 'AR6_regions':
                region_full_name = df_region[df_region[3]==region_abbrev][2].values[0]
            else:
                region_full_name = region_abbrev
            print(f'\n\n >>> @{region} << {region_full_name}')

            # --- figure ---
            fig, ax, ax2, ax_pos = generate_container(region_id, region_abbrev, region_full_name)

            cell_text = []
            dict_cdf = {}
            for rcp in rcps:

                # prepare empty dictionaries...
                dict_countor = {counter_threshold: 0 for counter_threshold in counter_thresholds}

                # load bootstrap regional results
                full_tfes = load_tfes_nc(region, rcp, season, target_chunk)  # TFE (member, K)
                # convert TFE into dGMT
                full_dgmts = conv_tfe2dgmt(full_tfes, rcp).flatten()  # dGMT: (member x K)

                # calculate stats
                bs_medians = bootstrap_medians(full_dgmts)  # (2000)
                #robust_overall_05, robust_overall_95 = np.percentile(full_dgmts, 5), np.percentile(full_dgmts, 95)
                #robust_median_05, robust_median_95 = np.percentile(bs_medians, 5), np.percentile(bs_medians, 95)
                #print('robust_overall_05={}, robust_overall_95={}'.format(robust_overall_05, robust_overall_95))
                #print('robust_median_05={}, robust_median_95={}'.format(robust_median_05, robust_median_95))

                # calc frac_tfe_member by 2050 and 2100
                for dgmt in full_dgmts:
                    for counter_threshold in counter_thresholds:
                        if dgmt <= counter_threshold:
                            dict_countor[counter_threshold] += 1
                _cell_text = []
                for counter_threshold in counter_thresholds:
                    frac_tfe_member = dict_countor[counter_threshold] / full_dgmts.shape[0]
                    _cell_text.append(np.round(frac_tfe_member, 2))
                    df2.loc[region, f'by{counter_threshold}.{rcp}.median'] = frac_tfe_member
                cell_text.append(_cell_text)

                # main CDF
                dict_cdf[rcp] = estimate_frequency(full_dgmts, rcp)  

                # 2000 sub-groups for an CDF uncertainty estimate
                rng.shuffle(full_dgmts)  # shuffle!!
                bs_cdfs = np.array([estimate_frequency(dgmt, rcp) for dgmt in full_dgmts.reshape(2000, 1000)])  # (2000, 19)
                bs_cdf_max = np.max(bs_cdfs, axis=0)
                bs_cdf_min = np.min(bs_cdfs, axis=0)
                for counter_threshold in counter_thresholds:
                    df2.loc[region, f'by{counter_threshold}.{rcp}.min'] = np.min(bs_cdfs[:, dict_dgmts[rcp].index(counter_threshold)])
                    df2.loc[region, f'by{counter_threshold}.{rcp}.max'] = np.max(bs_cdfs[:, dict_dgmts[rcp].index(counter_threshold)])
                # at which dgmt, CDF reach 95%
                target_dgmt = np.array([0.95]*len(dict_dgmts[rcp]))
                _dgmts, dummy = interpolated_intercepts(dict_dgmts[rcp], dict_cdf[rcp], target_dgmt)
                if not _dgmts.shape[0] == 0:
                    df3.loc[region, f'{rcp}.median'] = _dgmts[0]
                _dgmts, dummy = interpolated_intercepts(dict_dgmts[rcp], bs_cdf_min, target_dgmt)
                if not _dgmts.shape[0] == 0:
                    df3.loc[region, f'{rcp}.min'] = _dgmts[0]
                _dgmts, dummy = interpolated_intercepts(dict_dgmts[rcp], bs_cdf_max, target_dgmt)
                if not _dgmts.shape[0] == 0:
                    df3.loc[region, f'{rcp}.max'] = _dgmts[0]

                add_uncertainty_plot(ax, rcp, bs_cdf_max, bs_cdf_min)

            for rcp in rcps:
                add_main_plot(ax, rcp, dict_cdf[rcp])

            #add_table(fig, ax_pos, cell_text)

            savefig(target_chunk, region, season)

        df1.to_excel(writer1, sheet_name='medianTFE')
        writer1.save()
        writer1.close()

        df2.to_excel(writer2, sheet_name='CDatDGMT')
        df3.to_excel(writer2, sheet_name='DGMTat95%')
        writer2.save()
        writer2.close()

    print('Successfully DONE!! d(^o^)b')


if __name__ == '__main__':
    main()
