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

experiment = 'basic'

#ensemble_type = 'original'
ensemble_type = 'bootstrap'

#region_type = 'AR6_regions'
#region_type = 'HydroBASINS_lev1'
region_type = sys.argv[1]

sample_type = sys.argv[2]

trend_type = 'quadratic'

K = 100000

trend_syear = int(sys.argv[3])

blocksize = 5

figure_type = sys.argv[4]
#figure_type = 'main'
#figure_type = 'supplementary'

seasons = [sys.argv[5]]
#seasons = ['annual', 'DRY', 'WET']
#seasons = ['DRY']  #, 'WET']

target_chunks = [sys.argv[6]]
#target_chunks = [5, 10]

ghms = ['cwatm', 'h08', 'lpjml', 'matsiro', 'watergap2']
gcms = ['hadgem2-es',  'ipsl-cm5a-lr', 'gfdl-esm2m', 'miroc5']

rcps = ['rcp26', 'rcp85']
#rcps = ['rcp26']
#rcps = ['rcp85']

if TEST:
    ghms = ['cwatm']
    gcms = ['hadgem2-es']
    #rcps = ['rcp85']

threshold = 'max'

dict_color = {
    #'rcp26': '#0080ff',  # right blue
    #'rcp85': '#ff6600',  # orange
    'rcp26': '#61aac6',  # right blue 2     
    'rcp85': '#d98a30',  #'#ff5500',  # orange  '#dc0d0d',
    }

syear, eyear = 2005, 2099
years = range(syear, eyear+1)
counter_thresholds = [2050, 2100]
missing_value = 9999
years_5yrs = range(2010, 2100+1, 5)  # intervals for CDF estimation

main_directory = '/data/rg001/sgec0017/data/figure_box/drought_tpcd_isimip2b'
fig_directory_main = os.path.join(main_directory, 'draw_CDF_from_bootstrap')
mapmask_directory = '/data/rg001/sgec0017/data/mapmask'

# basic
if region_type == 'AR6_regions':
    dict_input_date = {                       # date of bootstrap
        ('Mean', 'annual', 1985): 20210707,
        ('Mean', 'DRY',    1985): 20210707,
        ('Mean', 'WET',    1985): 20210707,
        #('Mean', 'annual', 1985): 20210726,  # with DryGridMask
        #('Mean', 'DRY',    1985): 20210726,  # with DryGridMask
        #('Mean', 'WET',    1985): 20210726,  # with DryGridMask
        }
    regional_name_path = os.path.join(mapmask_directory, 'IPCC-AR6-WGI-reference-regions-v4_shapefile', 'regions.xlsx')
    df_region = pd.read_excel(regional_name_path, header=None)
elif region_type == 'HydroBASINS_lev1':
    dict_input_date = {                       # date of bootstrap
        ('Mean', 'annual', 1985): 20220103,
        ('Mean', 'DRY',    1985): 20220103,
        ('Mean', 'WET',    1985): 20220103,
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

# -----------------------------------------------------------------------------------------------------------------------
def load_tfes_nc(region, rcp, season, target_chunk):
    ncfile = f'resampled_tfe.{region}.nc4'
    input_directory = os.path.join(
        main_directory, 'bootstrap', region_type,
        f'{experiment}.{sample_type}.{trend_type}{str(trend_syear)[2:]}99.{K}.block{blocksize}.{season}',
        f'tChunk_{target_chunk:02}',
        str(dict_input_date[(sample_type,season,trend_syear)]), f'{rcp}_2005soc_co2'
        )
    nc_path = os.path.join(input_directory, ncfile)
    print(f'loading... {nc_path}')
    return Dataset(nc_path)['resampled_tfe'][:,:-1]  # (n_member, n_samples)


def estimate_frequency(input_tfes):
    total = input_tfes.shape
    cumulative_density = []
    for year_5yrs in years_5yrs:
        cumulative_density.append((np.where(input_tfes<=year_5yrs,1,0).sum()/total)[0])
    return cumulative_density


def bootstrap_medians(src):
    rng.shuffle(src)  # shuffle
    medians = np.median(src.reshape(2000, 1000), axis=1)  # (2000)
    return medians


# ----------------------------------------------------------------------------------------------------------------------
def main():

    region_map, dict_region = get_region_info(region_type, TEST=TEST)
    regions = [f'{region_id:02}.{region_name}' for region_id, region_name in dict_region.items()]

    for target_chunk, season in itertools.product(target_chunks, seasons): 

        target_chunk = int(target_chunk)

        info = f'{experiment}.{ensemble_type}{K}.{trend_type}{str(trend_syear)[2:]}99.{blocksize}.{sample_type}.{season}'
        figure_directory = os.path.join(fig_directory_main, region_type, info, figure_type)
        if not os.path.isdir(figure_directory): os.makedirs(figure_directory)

        excel_file = f'median_tfe_summary.{ensemble_type}.{season}.{target_chunk}.xlsx'
        excel_path = os.path.join(figure_directory, excel_file)
        writer = pd.ExcelWriter(excel_path)
        df = pd.DataFrame(index=regions, columns=rcps)

        robustlevel_file = f'robst_level_summary.{ensemble_type}.{season}.{target_chunk}.xlsx'
        excel_path_2 = os.path.join(figure_directory, robustlevel_file)
        writer_2 = pd.ExcelWriter(excel_path_2)
        # cumulative probability (member) reaching TFE by 2050 and 2099
        mmm = ['median', 'min', 'max']
        columns2 = [f'by{counter_threshold}.{rcp}.{m}' for rcp in rcps
                                                       for counter_threshold in counter_thresholds
                                                       for m in mmm]
        colmuns3 = [f'{rcp}.{m}' for rcp in rcps for m in mmm]
        df2 = pd.DataFrame(index=regions, columns=columns2)
        # time at which 95% member reach TFE
        df3 = pd.DataFrame(index=regions, columns=colmuns3)

        for region in regions:

            region_id = int(region[:2])
            region_abbrev = region[3:]
            if region_type == 'AR6_regions':
                region_longname = df_region[df_region[3]==region_abbrev][2].values[0]
            elif region_type == 'HydroBASINS_lev1':
                region_longname = df_region_longname.loc[region_id, 'long_name']
            else: raise ValueError
            print(f'\n\n >>> @{region} << {region_longname}')
        
            # --- figure ---
            # fixed parameters
            xs = years_5yrs
            _vmax = 2113
            if figure_type == 'main':
                figsize = (6, 6)
                fontsize_region_name = 22
                fontsize_ticks = 19.5
            elif figure_type == 'supplementary':
                figsize = (5, 6)
                fontsize_region_name = 19 
                fontsize_ticks = 17
            # ---
            fig = plt.figure(figsize=figsize)
            plt.subplots_adjust(left=0.15, bottom=0.052, right=0.96, top=0.945)
            ax = fig.add_subplot(111)
            ax.set_title(f'{region_id}. {region_longname}', fontsize=fontsize_region_name)
            for position in ['right', 'left']:  #, 'top', 'bottom']:
                ax.spines[position].set_visible(False)
            # y-axis
            #yticks = [0, 0.05, 0.25, 0.33, 0.5, 0.66, 0.75, 0.95, 1]
            yticks = [0, 0.33, 0.5, 0.66, 0.95, 1]
            yticklabels = yticks
            ax.set_yticks(yticks)
            ax.set_yticklabels(yticklabels, fontsize=fontsize_ticks)
            ax.set_ylim([0, 1])
            # x-axis
            xticks = [2010, 2030, 2050, 2070, 2090]
            xticklabels = xticks
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticklabels, fontsize=fontsize_ticks)
            ax.set_xlim([syear, eyear + 5])
            # horizontal lines at 0.05, 0.5, and 0.95
            #ax.axhline(y=0.05, color='#808080', linewidth=0.6, linestyle='-')  # gray
            ax.axhline(y=0.33, color='#808080', linewidth=0.6, linestyle='-')  # gray
            ax.axhline(y=0.50, color='#808080', linewidth=0.6, linestyle='-')  # gray
            ax.axhline(y=0.66, color='#808080', linewidth=0.6, linestyle='-')  # gray
            ax.axhline(y=0.95, color='#808080', linewidth=0.6, linestyle='-')  # gray
            # virtical lines at 2050
            ax.axvline(x=2050, color='#808080', linewidth=0.6, linestyle='-')  # gray
            # grid
            ax.grid(color='#e0e0e0')
            #
            ax_pos = ax.get_position()
            ax2 = fig.add_axes([ax_pos.x0, ax_pos.y0-0.2, ax_pos.x1-ax_pos.x0, 0.2])
            ax2.axis('off')
            ax2.set_ylim([-0.2, 0])
            ax2.set_xlim([syear, eyear + 1])

            cell_text = []
            dict_cdf = {}
            for rcp in rcps:

                # empty dictionaries...
                dict_countor = {counter_threshold: 0 for counter_threshold in counter_thresholds}

                # load bootstrap regional results
                full_samples = load_tfes_nc(region, rcp, season, target_chunk).flatten()  # (member x samples)
                print(f'({rcp}) full_samples: {full_samples.shape}  {full_samples.min()}-{full_samples.max()}')

                # calculate stats
                median = np.median(full_samples)
                robust_overall_05, robust_overall_95 = np.percentile(full_samples, 5), np.percentile(full_samples, 95)
                print(f'robust_overall_05={robust_overall_05}, robust_overall_95={robust_overall_95}')
                bs_medians = bootstrap_medians(full_samples)
                robust_median_05, robust_median_95 = np.percentile(bs_medians, 5), np.percentile(bs_medians, 95)
                print(f'robust_median_05={robust_median_05}, robust_median_95={robust_median_95}')
                # calc frac_tfe_member by 2050 and 2100
                for tfe in full_samples:
                    for counter_threshold in counter_thresholds:
                        if tfe <= counter_threshold:
                            dict_countor[counter_threshold] += 1
                _cell_text = []
                for counter_threshold in counter_thresholds:
                    frac_tfe_member = dict_countor[counter_threshold] / full_samples.shape[0]
                    _cell_text.append(np.round(frac_tfe_member, 2))
                    df2.loc[region, f'by{counter_threshold}.{rcp}.median'] = frac_tfe_member
                cell_text.append(_cell_text)

                cdf_main = estimate_frequency(full_samples)
                dict_cdf[rcp] = cdf_main

                # 2000 sub-groups for an uncertainty estimate
                rng.shuffle(full_samples)  # shuffle!!
                bs_cdfs = np.array([estimate_frequency(samples) for samples in full_samples.reshape(2000, 1000)])  # (2000, 19)
                bs_cdf_max = np.max(bs_cdfs, axis=0)
                bs_cdf_min = np.min(bs_cdfs, axis=0)
                for counter_threshold in counter_thresholds:
                    df2.loc[region, f'by{counter_threshold}.{rcp}.min'] = np.min(bs_cdfs[:, years_5yrs.index(counter_threshold)])
                    df2.loc[region, f'by{counter_threshold}.{rcp}.max'] = np.max(bs_cdfs[:, years_5yrs.index(counter_threshold)])
                # when CDF reach 95%
                _years, dummy = interpolated_intercepts(years_5yrs, cdf_main,   np.array([0.95]*len(years_5yrs)))
                if not _years.shape[0] == 0:
                    df3.loc[region, f'{rcp}.median'] = _years[0]
                _years, dummy = interpolated_intercepts(years_5yrs, bs_cdf_min, np.array([0.95]*len(years_5yrs)))
                if not _years.shape[0] == 0:
                    df3.loc[region, f'{rcp}.min'] = _years[0]
                _years, dummy = interpolated_intercepts(years_5yrs, bs_cdf_max, np.array([0.95]*len(years_5yrs)))
                if not _years.shape[0] == 0:
                    df3.loc[region, f'{rcp}.max'] = _years[0]

                # plot CDF uncertainty
                ax.fill_between(xs, bs_cdf_max, bs_cdf_min, color=dict_color[rcp], alpha=0.1)
                ax.plot(xs, bs_cdf_max, color=dict_color[rcp], linestyle='-', linewidth=0.2, markersize=0)
                ax.plot(xs, bs_cdf_min, color=dict_color[rcp], linestyle='-', linewidth=0.2, markersize=0)
                """
                # median error bars
                if rcp == 'rcp85': y_median_errorbar = -0.06
                elif rcp == 'rcp26': y_median_errorbar = -0.1
                ax2.plot([robust_median_05, robust_median_95], [y_median_errorbar, y_median_errorbar],
                         color=dict_color[rcp], linewidth=1.5, linestyle='--', marker='|')
                # overall error bars
                if rcp == 'rcp85': y_overall_errorbar = -0.08
                elif rcp == 'rcp26': y_overall_errorbar = -0.12
                ax2.plot([robust_overall_05, robust_overall_95], [y_overall_errorbar, y_overall_errorbar],
                         color=dict_color[rcp], linewidth=1.5, linestyle='-', marker='|')
                """

            for rcp in rcps:
                ax.plot(xs, dict_cdf[rcp], color=dict_color[rcp], linestyle='-', linewidth=2.5, markersize=0)

            """
            # table
            ax3 = fig.add_axes([ax_pos.x0+0.47, ax_pos.y1-0.75, 0.3, 0.15])
            ax3.axis('off')
            ax3.set_title('Likelihood')
            columns = ['by2050', 'by2100']
            #cell_text.reverse()
            table = ax3.table(cellText=cell_text,
                              rowLabels=rcps, colLabels=columns, 
                              cellColours=None, cellLoc='center',
                              edges='open',
                              loc='lower right')
            table.auto_set_font_size(False)
            table.set_fontsize(15)
            table.scale(1,2)
            ax3.axhline(y=1, color='#808080') 
            """
        
            for suffix in ['png', 'pdf']:
                figname = f'cdf_indevidual_split.{season}.chunk{target_chunk:03}.{region}.{suffix}'
                figdir = os.path.join(figure_directory, suffix)
                if not os.path.isdir(figdir): os.makedirs(figdir)
                figpath = os.path.join(figdir, figname)
                plt.savefig(figpath)
                print(f'savefig: {figpath}')
            plt.close()

        df.to_excel(writer, sheet_name='medianTFE')
        writer.save()
        writer.close()

        df2.to_excel(writer_2, sheet_name='CDatYear')
        df3.to_excel(writer_2, sheet_name='YEARat95%')
        writer_2.save()
        writer_2.close()

    print('Successfully DONE!! d(^o^)b')


if __name__ == '__main__':
    main()
