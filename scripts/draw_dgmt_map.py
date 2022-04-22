#!/usr/bin/env python
# To draw dGMT maps corresponding to tfe based on bootstrap data
# By Yusuke Satoh
# On 2021/6/29

import os
import sys
import pandas as pd
import numpy as np
import itertools
import datetime
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors, colorbar
from matplotlib import mathtext
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
from netCDF4 import Dataset
from utiltools import get_region_info
mathtext.FontConstantsBase = mathtext.ComputerModernFontConstants
plt.rcParams['font.family'] = 'Arial'
today = datetime.date.today().strftime('%Y%m%d')


TEST = False

experiment = 'basic'

#region_type = 'AR6_regions'
#region_type = 'HydroBASINS_lev1'
region_type = sys.argv[1]

sample_types = [sys.argv[2]]
#sample_types = ['Mean']
#sample_types = ['Average']
#sample_types = ['Extreme']

#ensemble_type = 'original'  # oroginal 20 ensemble member
ensemble_type = 'bootstrap'  # bootstrap   

trend_type = 'quadratic'

trend_syear = sys.argv[3]
#trend_syear = 2005

seasons = [sys.argv[4]]
#seasons = ['annual', 'DRY', 'WET']

target_chunk = int(sys.argv[5]) # tChunk

#K = 0
K = 100000

blocksize = 5

rcps = ['rcp26', 'rcp85']
#rcps = ['rcp26']
#rcps = ['rcp85']

if TEST:
    seasons = ['DRY']
    rcps = ['rcp85']

#colorscheme = 'YlOrRd'
#colorscheme = 'YlOrRd_r'
#colorscheme = 'hot'
#colorscheme = 'hot_r'
colorscheme = 'YlOrRd'
#colorscheme = 'copper'

main_directory = '/data/rg001/sgec0017/data/figure_box/drought_tpcd_isimip2b'
data_directory = '/data/rg001/sgec0017/data'
srcdirectory_tdgmt =  os.path.join(data_directory, 'figure_box', 'drought_tpcd_isimip2b', 'extract_temperature')
figure_directory_main = os.path.join(data_directory, 'figure_box', 'drought_tpcd_isimip2b', 'draw_dgmt_map')
mapmask_directory = os.path.join(data_directory, 'mapmask')
landseamask_path = os.path.join(mapmask_directory, 'ISIMIP2b_landseamask_generic.nc4')

landsea = Dataset(landseamask_path)['LSM'][:][0].filled(0)
landseamask = np.ma.make_mask(Dataset(landseamask_path).variables['LSM'][:][0] != 1.0)

if region_type == 'AR6_regions':

    dict_input_date = {  # date of bootstrap
        #('Average', 'annual'): 20210701,
        #('Average', 'DRY'): 20210625,
        #('Average', 'WET'): 20210630,
        ('Mean', 'annual'): 20210707,
        ('Mean', 'DRY'   ): 20210707,
        ('Mean', 'WET'   ): 20210707,
        #('Mean', 'annual'): 20210726,
        #('Mean', 'DRY'   ): 20210726,
        #('Mean', 'WET'   ): 20210726,
        }
    region_shp_path = os.path.join(
        mapmask_directory, 'IPCC-AR6-WGI-reference-regions-v4_shapefile', 
        'shapefile_edited', 'land', 'reference-regions-AR6'
        )
    regions = [
        '01.NWN',  '02.NEN', '03.WNA',  '04.CNA',  '05.ENA',
        '06.NCA',  '07.SCA', '08.CAR',  '09.NWS',  '10.NSA',
        '11.NES',  '12.SAM', '13.SWS',  '14.SES',  '15.SSA',
        '16.NEU',  '17.WCE', '18.EEU',  '19.MED',  '20.SAH',
        '21.WAF',  '22.CAF', '23.NEAF', '24.SEAF', '25.WSAF',
        '26.ESAF', '27.MDG', '28.RAR',  '29.WSB',  '30.ESB',
        '31.RFE',  '32.WCA', '33.ECA',  '34.TIB',  '35.EAS',
        '36.ARP',  '37.SAS', '38.SEA',  '39.NAU',  '40.CAU',
        '41.EAU',  '42.SAU', '43.NZ'
        ]
elif region_type == 'HydroBASINS_lev1':

    dict_input_date = {  # date of bootstrap
        ('Mean', 'annual'): 20211204,
        ('Mean', 'DRY'   ): 20211204,
        ('Mean', 'WET'   ): 20211204,
        }
    region_map, dict_regions = get_region_info(region_type, TEST=False)
    del region_map
    regions = [f'{region_id:02}.{region_name}' for region_id, region_name in dict_regions.items()]
else:
    raise ValueError(f'check region_type  {region_type}')
print(f'regions: {regions}')

ghms = ['cwatm', 'h08', 'lpjml', 'matsiro', 'watergap2']
gcms = ['hadgem2-es',  'ipsl-cm5a-lr', 'gfdl-esm2m', 'miroc5']
members = [f'{ghm}_{gcm}' for ghm in ghms for gcm in gcms]

dict_rcp = {
    'rcp26': '#0080ff',  # right blue
    'rcp85': '#ff6600',  # orange
    }

suffixes = ['png', 'pdf']

missing_value = 1e+20



# ----------------------------------------------------------------------------------------------------------------------
def load_regional_tfes(sample_type, region, rcp, season):
    nc_file = f'resampled_tfe.{region}.nc4'
    nc_directory = os.path.join(
        main_directory, 'bootstrap', region_type,
        f'{experiment}.{sample_type}.{trend_type}{trend_syear[2:]}99.{K}.block{blocksize}.{season}', 
        f'tChunk_{target_chunk:02}', 
        str(dict_input_date[(sample_type,season)]), f'{rcp}_2005soc_co2'
        )
    nc_path = os.path.join(nc_directory, nc_file)
    return Dataset(nc_path)['resampled_tfe'][:]


def load_tfe_regions_excel(sample_type, season):
    info = f'{experiment}.{ensemble_type}{K}.{trend_type}{trend_syear[2:]}99.{blocksize}.{sample_type}.{season}'
    excel_file = f'tfe_summary.{info}.median.{target_chunk}.xlsx'
    excel_directory = os.path.join(main_directory, 'draw_tfe_map_from_bootstrap', region_type, info) 
    excel_path = os.path.join(excel_directory, excel_file)
    print(f'loading... {excel_path}')
    return pd.read_excel(excel_path, index_col=0)


def if_too_late(tpcd, threshold_2084=False):
    if tpcd == 'nan':
        tpcd = np.nan
    elif np.isnan(tpcd):
        tpcd = np.nan
    else:
        if 1e+5 < tpcd:   tpcd = np.nan
        else: tpcd = int(tpcd)
    return tpcd


def getNearestValues(list_src, num):
    idx = np.abs(np.asarray(list_src) - num).argmin()
    return list_src[idx], list_src[idx+1] 


# ----------------------------------------------------------------------------------------------------------------------
def main(*argv):

    region_map, dict_region = get_region_info(region_type, TEST=TEST)

    # load TdGMT excel file
    tdgmt_excelname = 'isimip2b_temperature_timeseries__baseperiod_1850-1900.xlsx'
    tdgmt_excel_path = os.path.join(srcdirectory_tdgmt, tdgmt_excelname)
    df_dgmt_org = pd.read_excel(tdgmt_excel_path, sheet_name='dGMT', index_col=0)

    for sample_type, season in itertools.product(sample_types, seasons):

        # use regional info from Fig2.draw_tfe_map_from_bootstrap.py
        df_tfe_regions = load_tfe_regions_excel(sample_type, season)

        info = f'{experiment}.{ensemble_type}{K}.{trend_type}{trend_syear[2:]}99.{blocksize}.{sample_type}.{season}'
        figure_directory = os.path.join(figure_directory_main, region_type, info)
        if not os.path.isdir(figure_directory): os.makedirs(figure_directory)

        # output excel for logging
        excelfile = f'tpcd_dgmd.{sample_type}.{season}.xlsx'
        outpath = os.path.join(figure_directory, excelfile)
        writer = pd.ExcelWriter(outpath)

        for rcp in rcps:
            print(f'{sample_type}  {rcp}')
            
            # create empties
            # this refers overall ensemble Mean dGMT for a median TFE
            df_dgmt = pd.DataFrame(
                            index=regions,                                        
                            columns=['medianTFE',                  # emsemble
                                     'GMT_for_medianTFE',          # original GMT that cooresponds to the median TFE
                                     'median_of_eachGMT_for_all',  # estimate member specific GMT for each TFEs, then take their mean
                                     'mean_of_eachGMT_for_TFEs',   # estimate member specific GMT for each TFEs, then take their mean
                                     ]
                                    )

            dict_TFEregion_rate = {}
            dict_TFEregion_less2deg_rate = {}
            for region in regions:
                if not region in df_tfe_regions[df_tfe_regions[f'medi-{rcp}']==1].index:
                    print(f'>>> {region} pass...')
                    df_dgmt.loc[region, 'medianTFE'] = missing_value
                    df_dgmt.loc[region, 'GMT_for_medianTFE'] = missing_value
                    df_dgmt.loc[region, 'median_of_eachGMT_for_all'] = missing_value
                    df_dgmt.loc[region, 'mean_of_eachGMT_for_TFEs'] = missing_value
                    dict_TFEregion_rate[region] = missing_value
                    dict_TFEregion_less2deg_rate[region] = missing_value

                else:  # process specific regions with the robust median (=regions mapped in Fig2a)
                
                    # load regional TFE netCDF files here
                    tfes = load_regional_tfes(sample_type, region, rcp, season)  # (nmember, nK) from .nc
                    
                    # conv tfe to dgmt
                    dgmts = []
                    _tfes = []
                    tpcd_counter = 0
                    tpcd_less2deg_counter = 0
                    for (ighm, ghm), (igcm, gcm), k in itertools.product(enumerate(ghms), enumerate(gcms), range(K+1)):
                        column = f'{rcp}_{gcm}'
                        tfe = int(tfes[len(gcms)*ighm+igcm,k])
                        _tfes.append(tfe)
                        if tfe <= 2084:
                            dgmt = df_dgmt_org.loc[str(tfe), column]
                            dgmts.append(dgmt)
                        elif 2085 <= tfe <= 2099:
                            dgmt = df_dgmt_org.loc['2084', column]
                            dgmts.append(dgmt)
                        else:
                            dgmt = missing_value
                            dgmts.append(dgmt)
                        if tfe < 2100: tpcd_counter += 1
                        if dgmt <= 2: tpcd_less2deg_counter += 1
                    # overall median
                    medianTFE = np.median(_tfes)  

                    # find GMT for the median TFE
                    if medianTFE < 2099:  # this should be True because of the flug of df_tfe_regions, but just for the case...
                        n_indices = _tfes.count(medianTFE)

                        if n_indices == 1:  # the overall median is derived from one ensemble member.
                            GMT_for_medianTFE = dgmts[_tfes.index(medianTFE)]
                        elif n_indices == 0:  # the median is an average of two values (Caution: some ensemble members have one value...)
                            _tfes_sorted = sorted(_tfes)
                            nearest_value_small = getNearestValues(_tfes_sorted, medianTFE)[0]   # take the only smallest one
                            index_nearest_value_small = _tfes_sorted.index(nearest_value_small)  # get the index of the "first" smaller nearest value
                            n_nearest_value_small = _tfes_sorted.count(nearest_value_small)      # count how many the smaller neares value you have
                            nearest_value_big = _tfes_sorted[index_nearest_value_small+n_nearest_value_small]  # get the bigger neares value
                            _tfes = np.array(_tfes)
                            nearest_indices = np.where(_tfes==nearest_value_small)[0].tolist() + np.where(_tfes==nearest_value_big)[0].tolist()
                            GMT_for_medianTFE = np.mean([dgmts[nearest_index] for nearest_index in nearest_indices])
                        elif n_indices >= 2:
                            _tfes = np.array(_tfes)
                            median_indices = np.where(_tfes==medianTFE)[0].tolist()
                            GMT_for_medianTFE = np.mean([dgmts[median_index] for median_index in median_indices])
                        else:
                            raise ValueError
                    else:
                        GMT_for_medianTFE = None

                    # media TFE
                    df_dgmt.loc[region, 'medianTFE'] = medianTFE
                    # gmts
                    df_dgmt.loc[region, 'GMT_for_medianTFE'] = GMT_for_medianTFE
                    df_dgmt.loc[region, 'median_of_eachGMT_for_all'] = np.median(dgmts)   # <---- 50percentile in the CDF
                    df_dgmt.loc[region, 'mean_of_eachGMT_for_TFEs'] = np.ma.masked_equal(dgmts, 1e+20).mean()

                    if ensemble_type == 'original': 
                        N = len(ghms)*len(gcms)
                    elif ensemble_type == 'bootstrap':
                        N = len(ghms)*len(gcms)*(K+1)
                    dict_TFEregion_rate[region] = tpcd_counter / N
                    dict_TFEregion_less2deg_rate[region] = tpcd_less2deg_counter / N

                    print(f'>>> {region} processed  year_{int(medianTFE)} --> {GMT_for_medianTFE:.2}_deg')

            # draw global mean temperature rise maps
            figsize = (8, 3.6)
            dpi_figure = 300
            dpi_savefig = 300
            fontsize = 8
            labelsize = 8
            # color info
            bounds = [1, 1.5, 2, 2.5, 3, 3.5, 4]
            bounds_ticklabels = bounds
            colorlist = list(plt.get_cmap(colorscheme, len(bounds)+2)(range(len(bounds)+2)))
            cmap = colors.ListedColormap(colorlist[1:-1], "")
            norm = colors.BoundaryNorm(bounds, cmap.N)
            cmap.set_under(colorlist[0])
            cmap.set_over(colorlist[-1])
            # --- main
            for dGMT_type in ['GMT_for_medianTFE', 
                              #'median_of_eachGMT_for_all', 
                              #'mean_of_eachGMT_for_TFEs'
                              ]:
                print(f'--- {dGMT_type} ---')
                # initialize an empty field
                src = np.full((360, 720), missing_value)  

                # filling a global map with results
                for region in regions:
                    region_id = int(region.split('.')[0])
                    YY, XX = np.where(region_map==region_id)
                    gmt = df_dgmt.loc[region, dGMT_type]
                    src[YY,XX] = gmt
                    print(f'{region}: {gmt}')
                # masking
                src = np.ma.masked_equal(src, missing_value)
                src = np.ma.masked_array(src, mask=landseamask)
                src = src[:293]
                # draw a figure
                fig = plt.figure(num=1, figsize=figsize, dpi=dpi_figure)
                fig.subplots_adjust(left=0.005, bottom=0.02, right=0.995, top=0.98)
                ax = fig.add_subplot(1,1,1, projection=ccrs.Robinson())
                ax.set_extent([-180, 180, -60, 90], ccrs.PlateCarree())
                ax.outline_patch.set_linewidth(0.3)
                ax.add_feature(cfeature.LAND, facecolor='#d4d4d4')
                ax_pos = ax.get_position()
                lons, lats = np.arange(-179.75, 179.75+0.5, 0.5), np.arange(89.75, -56.25-0.5, -0.5)
                im = plt.contourf(lons, lats, src, levels=bounds, extend='both', norm=norm, cmap=cmap, transform=ccrs.PlateCarree())
                ax.add_feature(cfeature.COASTLINE, linewidth=0.3, edgecolor='#000000')
                # regional shape file
                shape_feature = ShapelyFeature(Reader(region_shp_path).geometries(), ccrs.PlateCarree(), 
                                               edgecolor='#808080', linewidth=0.3, facecolor='none')
                ax.add_feature(shape_feature)

                cax = fig.add_axes([ax_pos.x0+0.43, ax_pos.y0 + 0.08, 0.30, 0.03])
                cb = colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, ticks=bounds_ticklabels, 
                                           orientation='horizontal', extend='both')
                cb.ax.tick_params(labelsize=labelsize, direction='in')
                cb.outline.set_visible(False)
                cax.text(1.06, -0.98, '$[^{o}C]$', 
                         ha='left', va='bottom', fontsize=fontsize, 
                         horizontalalignment='left', verticalalignment='bottom',
                         transform=cax.transAxes
                         )

                for suffix in suffixes:
                    figname = f'map.{target_chunk:03}.{sample_type}.{rcp}.{dGMT_type}.{suffix}'
                    _fig_directory = os.path.join(figure_directory, suffix)
                    if not os.path.isdir(_fig_directory): os.makedirs(_fig_directory)
                    figpath = os.path.join(_fig_directory, figname)
                    plt.savefig(figpath, dpi=dpi_savefig)
                    print(f'\nsavefig: {figpath}')
                plt.close()

            df_dgmt.to_excel(writer, sheet_name=rcp)

        writer.save()
        writer.close()

        print(f'\nwrite out: {outpath}')


if __name__=='__main__':
    main(*sys.argv)
