#!/usr/bin/env python
# To draw TFE maps based on bootstrap data 
# By Yusuke Satoh
# On 2021/4/9

import os
import sys
import itertools
import pathlib
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import colors, cm
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
from netCDF4 import Dataset
from make_referenceRegions import regenerate_shp
from utiltools import get_region_info
plt.rcParams['axes.linewidth'] = 0.5
plt.rcParams['font.family'] = 'Arial'
mpl.rcParams['hatch.linewidth'] = 0.055
#mpl.rcParams['hatch.linewidth'] = 0.05

TEST = False

experiment = 'basic'
#experiment = 'picontrol'
#experiment = 'co2exp'
#experiment = 'rcp26soc'

#ensemble_type = 'original'  # original 5ghm x 4gcm
ensemble_type = 'bootstrap'

#region_type = 'AR6_regions'
#region_type = 'HydroBASINS_lev1'
region_type = sys.argv[1]

sample_type = sys.argv[2]

trend_type = 'quadratic'

# sample size of bootstrap
#K = 1000  
K = 100000

trend_syear = sys.argv[3]

blocksize = 5

SHUFFLE = None

seasons = [sys.argv[4]]
#seasons = ['annual', 'DRY', 'WET']
#seasons = ['DRY', 'WET', 'annual']
#seasons = ['WET']
#seasons = ['DRY']

tChunks = [sys.argv[5]]
#tChunks = [5, 10]
#tChunks = [1, 4, 5, 6, 10, 15, 20]

#MEDIAN_HATCH = True
MEDIAN_HATCH = False

case = 'case_1'

if region_type == 'AR6_regions':
    dict_input_date = {                 # date of bootstrap
        #('Average', 'annual'): 20210701,
        #('Average', 'DRY'): 20210629,
        #('Average', 'WET'): 20210630,
        ('Mean', 'annual'): 20210707,
        ('Mean', 'DRY'   ): 20210707,
        ('Mean', 'WET'   ): 20210707,
        #('Mean', 'annual'): 20210726,  # with DryDridMask
        #('Mean', 'DRY'   ): 20210726,  # with DryDridMask
        #('Mean', 'WET'   ): 20210726,  # with DryDridMask
        }
elif region_type == 'HydroBASINS_lev1':
    dict_input_date = {                 # date of bootstrap
        ('Mean', 'annual'): 20220103,
        ('Mean', 'DRY'   ): 20220103,
        ('Mean', 'WET'   ): 20220103,
        }
else:
    raise ValueError(f'check region_type... {region_type}')

rcps = ['rcp26', 'rcp85']
#rcps = ['rcp26']
#rcps = ['rcp85']

socs = ['2005soc']
co2s = ['co2']
if experiment == 'basic':
    pass
elif experiment == 'picontrol':
    rcps = ['picontrol', 'rcp26', 'rcp85']
elif experiment =='rcp26soc':
    socs = ['2005soc', 'rcp26soc']
    rcps = ['rcp26']
elif experiment == 'co2exp':
    co2s = ['co2', '2005co2']

if TEST:
    tChunks = [5]
    seasons = ['DRY']
    rcps = ['rcp85']

missing_value = 9999

TFE_cmap = 'hot'
index = 'nDayTot'
drought_paras = 'Q80win15Len30tau4'
threshold_type = 'max'
soc = '2005soc'
stats = ['median']
#stats = ['median', 'Q05', 'Q95']
#stats = ['median', 'Q25', 'Q75']

#suffixes = ['pdf']
suffixes = ['png', 'pdf']
#suffixes = ['png', 'eps', 'pdf']
dpi_figure  = 300
dpi_savefig = 300


main_directory = '/data/rg001/sgec0017/data/figure_box/drought_tpcd_isimip2b'
mapmask_directory = '/data/rg001/sgec0017/data/mapmask'
fig_directory_main = os.path.join(main_directory, 'draw_tfe_map_from_bootstrap')
landseamask_path = os.path.join(mapmask_directory, 'ISIMIP2b_landseamask_generic.nc4')
landseamask = Dataset(landseamask_path)['LSM'][:][0].mask
if region_type == 'AR6_regions':
    region_shp_path = os.path.join(mapmask_directory, 'IPCC-AR6-WGI-reference-regions-v4_shapefile', 'shapefile_edited', 'land', 'reference-regions-AR6')
elif region_type == 'SREX_regions':
    region_shp_path = os.path.join(mapmask_directory, 'SREX.referenceRegions.editted', 'referenceRegions')
elif region_type == 'HydroBASINS_lev1':
    region_shp_path = os.path.join(mapmask_directory, 'HydroSHEDS', 'HydroBASINS', 'withLakes', 'Global_merge', 'hybas_lake____lev02_v1c_merge')
elif region_type == 'hydroregions':
    region_shp_path = os.path.join(mapmask_directory, 'Hydrobelt', 'hydroregions_shp', 'meybeck_et_al_2013_hydroregions')

SEED = 1
np.random.seed(SEED)
rng = np.random.default_rng()


def get_score(_stats, region_id, region_name, rcp, soc, co2, tChunk, season):
    # load bootstrapped data
    nc_file = f'resampled_tfe.{region_id:02}.{region_name}.nc4'
    input_directory = os.path.join(main_directory, 'bootstrap', region_type,
                                   f'{experiment}.{sample_type}.{trend_type}{trend_syear[2:]}99.{K}.block{blocksize}.{season}', 
                                   f'tChunk_{tChunk:02}', 
                                   str(dict_input_date[(sample_type,season)]), f'{rcp}_{soc}_{co2}')
    nc_path = os.path.join(input_directory, nc_file)

    if ensemble_type == 'original':  # tentatively, using bootstrap results.
        score = np.median(Dataset(nc_path)['resampled_tfe'][:,0])  # (20) --> single value
    elif ensemble_type == 'bootstrap':
        if _stats == 'median': score = np.median(Dataset(nc_path)['resampled_tfe'][:])  # (20, 100001) --> single value
        elif _stats == 'Q05': score = np.percentile(Dataset(nc_path)['resampled_tfe'][:],  5)
        elif _stats == 'Q95': score = np.percentile(Dataset(nc_path)['resampled_tfe'][:], 95)
        elif _stats == 'overall_Q95': score = np.percentile(Dataset(nc_path)['resampled_tfe'][:], 95)
        elif _stats == 'median_Q95':
            src = Dataset(nc_path)['resampled_tfe'][:]  # (20, 100001)
            src = src[:,:-1].flatten()                  # (20, 100000) --> (2000000)
            rng.shuffle(src)                            # shuffle
            src = np.median(src.reshape(2000, 1000), axis=1)  # (2000)
            score = np.percentile(src, 95)
    return score


def draw_tfe_map(src, robust_overall, robust_median, rcp, soc, co2, _stats, season, tChunk, fig_directory):  # Global map of TFE

    src = np.ma.masked_array(src, mask=landseamask)
    src = src[:293]  # cut out the target range
    src = np.ma.masked_greater(src, 2100-tChunk)
    print(f'src: {src.min()}-{src.max()}')
    # robust hatch
    robust_overall = np.ma.masked_array(robust_overall, mask=landseamask)
    robust_overall = robust_overall[:293]  # cut out the target range
    robust_median = np.ma.masked_array(robust_median, mask=landseamask)
    robust_median = robust_median[:293]  # cut out the target range
   
    # size
    figsize = (8, 3.6)
    fontsize = 10
    labelsize = 8

    fig1 = plt.figure(num=1, figsize=figsize, dpi=dpi_figure)
    fig1.subplots_adjust(left=0.005, bottom=0.02, right=0.995, top=0.98)
    ax = fig1.add_subplot(1,1,1, projection=ccrs.Robinson())
    ax.set_extent([-180, 180, -60, 90], ccrs.PlateCarree())
    ax.outline_patch.set_linewidth(0.3)
    ax.add_feature(cfeature.LAND, facecolor='#d4d4d4')
    ax_pos = ax.get_position()

    # TPCD map
    cmap = cm.get_cmap(TFE_cmap)
    vmin, vmax = 2010, 2100
    # boundaries = [vmin, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, vmax]
    boundaries = [vmin, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100, 2110]  # added extra to avoid white in the colorbar
    norm = colors.BoundaryNorm(boundaries, cmap.N)

    ##img_extent = [-180, 180, -60, 90]
    ##im = ax.imshow(src, origin='upper', vmin=int(vmin), vmax=int(vmax), norm=norm, cmap=cmap, transform=ccrs.PlateCarree(), extent=img_extent)
    lons, lats = np.arange(-179.75, 179.75+0.5, 0.5), np.arange(89.75, -56.25-0.5, -0.5)
    im = plt.contourf(lons, lats, src, levels=boundaries, norm=norm, cmap=cmap, transform=ccrs.PlateCarree())

    # robust_overall hatching   Caution!! Hatching is supported in the PostScript, PDF, SVG and Agg backends only.
    if MEDIAN_HATCH: hatch_type = r'\\\\'
    else:            hatch_type = r'///'
    im2 = plt.contourf(lons, lats, robust_overall, levels=3, hatches=[None, hatch_type], alpha=0, transform=ccrs.PlateCarree())
    for collection in im2.collections:
        collection.set_edgecolor('#d5d5d5')
        collection.set_linewidth(0.)
    # robust_median hatching
    if MEDIAN_HATCH:
        im3 = plt.contourf(lons, lats, robust_median, levels=3, hatches=[None, '///'], alpha=0, transform=ccrs.PlateCarree())
        for collection in im3.collections:
            collection.set_edgecolor('#d5d5d5')
            collection.set_linewidth(0.)

    # coastline
    if not region_type in ['HydroBASINS_lev1', 'hydroregions']:
        ax.add_feature(cfeature.COASTLINE, linewidth=0.3, edgecolor='#000000')

    # regional shape file
    shape_feature = ShapelyFeature(Reader(region_shp_path).geometries(), ccrs.PlateCarree(), edgecolor='#808080', linewidth=0.2, facecolor='none')
    ax.add_feature(shape_feature)
    # regional shape file (overall robust)
    if region_type == 'AR6_regions':
        region_shp_path2 = os.path.join(mapmask_directory, 'IPCC-AR6-WGI-reference-regions-v4_shapefile',
                                      'shapefile_edited', f'tfe{tChunk:02}overallrobust_regions_{season}_{rcp}',
                                      'reference-regions-AR6')
    elif region_type == 'HydroBASINS_lev1':
        region_shp_path2 = os.path.join(mapmask_directory, 'HydroSHEDS', 'HydroBASINS', 'withLakes', 'Global_merge', 
                                      'shapefile_edited', f'tfe{tChunk:02}overallrobust_regions_{season}_{rcp}',
                                      'hydrobasins_lev1')
    elif region_type == 'hydroregions':
        region_shp_path2 = os.path.join(mapmask_directory, 'Hydrobelt',
                                        'hydroregions_shp_edited', f'tfe{tChunk:02}overallrobust_regions_{season}_{rcp}',
                                        'hydroregions')
    else:
        raise ValueError
    shape_feature2 = ShapelyFeature(Reader(region_shp_path2).geometries(), ccrs.PlateCarree(), edgecolor='#000000', linewidth=0.8, facecolor='none')
    ax.add_feature(shape_feature2)

    # add a color bar
    cax = fig1.add_axes([ax_pos.x0 + 0.4, ax_pos.y0 + 0.08, 0.37, 0.03])
    cb = fig1.colorbar(im, cax=cax, orientation='horizontal')
    cb.ax.set_title(f'$TFE_{tChunk}$', fontsize=fontsize, pad=0.05)
    ticks_year = range(2020, 2100, 10)
    ticklabels_year = [2020, '', 2040, '', 2060, '', 2080, '']
    cb.set_ticks(ticks_year)

    cb.set_ticklabels(ticklabels_year)
    cb.ax.tick_params(labelsize=labelsize, width=0.2, direction='in')
    cb.outline.set_visible(False)

    if ensemble_type == 'original':
        figname = f'TFE.{season}.{tChunk:03}.{rcp}_{soc}_{co2}.{ensemble_type}.{_stats}'
        _info = f'{drought_paras}{season}.{tChunk:03}{threshold_type}.{rcp}_{soc}_{co2}.{_stats}.{ensemble_type}.{case}'
    elif ensemble_type == 'bootstrap':
        figname = f'TFE.{season}.{tChunk:03}.{rcp}_{soc}_{co2}.{ensemble_type}{_stats}.{trend_type}{trend_syear[2:]}99.K{K}.block{blocksize}.{case}.'
        _info = f'{drought_paras}{season}.{tChunk:03}{threshold_type}.{rcp}_{soc}_{co2}.{_stats}.{ensemble_type}.{trend_type}{trend_syear[2:]}99.K{K}.block{blocksize}.{case}'
    for suffix in suffixes:
        _fig_directory = os.path.join(fig_directory, f'medianhatch{str(MEDIAN_HATCH)[0]}', suffix)
        if not os.path.isdir(_fig_directory): os.makedirs(_fig_directory)
        figpath = os.path.join(_fig_directory, figname + suffix)
        plt.savefig(figpath, dpi=dpi_savefig)
        print(f'savefig: {figpath}')
        pathlib.Path(os.path.join(_fig_directory, f'originalfile_{_info}')).touch()

    plt.close(1)


# ----------------------------------------------------------------------------------------------------------------------
def main():

    region_map, dict_region = get_region_info(region_type, TEST=TEST)

    for _stats, season, tChunk in itertools.product(stats, seasons, tChunks):
        tChunk = int(tChunk)
        print(f'>>>>>> {season} {_stats} {tChunk}')
        info = f'{experiment}.{ensemble_type}{K}.{trend_type}{trend_syear[2:]}99.{blocksize}.{sample_type}.{season}'
        fig_directory = os.path.join(fig_directory_main, region_type, info)
        if not os.path.isdir(fig_directory): os.makedirs(fig_directory)

        summary_excel = f'tfe_summary.{info}.{_stats}.{tChunk}.xlsx'
        excel_path = os.path.join(fig_directory, summary_excel)
        writer = pd.ExcelWriter(excel_path)
        df = pd.DataFrame(index=[f'{region_id:02}.{region_name}' for region_id, region_name in dict_region.items()], 
                          columns=rcps+[f'medi-{rcp}' for rcp in rcps]+[f'all-{rcp}' for rcp in rcps])

        for rcp, soc, co2 in itertools.product(rcps, socs, co2s):
            
            print(f'\n------ {rcp} {soc} {co2} {_stats} {tChunk} {season}')
            # prepare empty field
            tfe_map            = np.full((360, 720), missing_value)  # generate empty with missing_value for TFE field
            robust_overall_map = np.full((360, 720), missing_value)  # generate empty with missing_value for robust_overall 1-0 field
            robust_median_map  = np.full((360, 720), missing_value)  # generate empty with missing_value for robust_median 1-0 field

            # handle bootstrapped tfe values
            df_regions_with_median_tfe_valid = pd.DataFrame(columns=['id', 'name'])  # a table of regions with color in Fig2a
            robust_overall_regions = []
            for region_id, region_name in dict_region.items():
                region = f'{region_id:02}.{region_name}'
                YY, XX = np.where(region_map==region_id)
                # median
                tfe_score = get_score(_stats, region_id, region_name, rcp, soc, co2, tChunk, season)
                tfe_score = int(tfe_score)
                tfe_map[YY, XX] = tfe_score  # plug in overall_median TFE for the map
                if tfe_score == 9999: tfe_score = ''
                df.loc[region, rcp] = tfe_score
                # check if median is robust (=if 95% of bootstrapped ensemble median has tfe) 
                robust_median_level = 1 if get_score('median_Q95', region_id, region_name, rcp, soc, co2, tChunk, season) <= 2099 else 0
                robust_median_map[YY, XX] = robust_median_level  # 1 or 0
                df.loc[region, f'medi-{rcp}'] = robust_median_level
                # check robustness with overall (if 95% of the large resample tfe is < 2099-tChunk)
                robust_overall_level = 1 if get_score('overall_Q95', region_id, region_name, rcp, soc, co2, tChunk, season) <= 2099 else 0
                robust_overall_map[YY, XX] = robust_overall_level  # 1 or 0
                df.loc[region, f'all-{rcp}'] = robust_overall_level
                if robust_overall_level == 1:
                    robust_overall_regions.append(region_name) 
                print(f'{rcp} {region_id:02}.{region_name:<4} >> {tfe_score:4} (median:{robust_median_level}, overall:{robust_overall_level})')

            # maskout regions where the meidan result is NOT robust, when MEDIAN_HATCH=False.
            if not MEDIAN_HATCH:
                tfe_map = np.where(robust_median_map==1, tfe_map, missing_value)
            # maskout missing_value
            tfe_map            = np.ma.masked_equal(tfe_map, missing_value)
            robust_median_map  = np.ma.masked_equal(robust_median_map, missing_value)
            robust_overall_map = np.ma.masked_equal(robust_overall_map, missing_value)

            # regenerate a shape file for robust regions
            shp_title = f'tfe{tChunk:02}overallrobust_regions_{season}_{rcp}'
            regenerate_shp(region_type, shp_title, robust_overall_regions)

            if not np.all(tfe_map.mask):
                draw_tfe_map(tfe_map, robust_overall_map, robust_median_map, rcp, soc, co2, _stats, season, tChunk, fig_directory)

        df.to_excel(writer, sheet_name=f'medianTFE({ensemble_type})')
        writer.save()
        writer.close()
        print(f'save: {excel_path}')

    print('Successfully DONE!! d(^o^)b')


if __name__ == '__main__':
    main()
