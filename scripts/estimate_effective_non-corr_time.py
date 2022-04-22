#!/usr/bin/env python
# To
# By Yusuke Satoh
# On
import os
import sys
import itertools
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.basemap import Basemap
from matplotlib.ticker import PercentFormatter
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
from scipy import stats
from sklearn.utils import resample
from utiltools import get_region_info


#region_type = 'AR6_regions'
region_type = 'HydroBASINS_lev1'
input_date = 20220102

rcps = ['rcp85', 'rcp26']
gcms = ['hadgem2-es', 'ipsl-cm5a-lr', 'gfdl-esm2m', 'miroc5']  # full members
ghms = ['cwatm', 'h08', 'lpjml', 'matsiro', 'watergap2']  # full members
region_map, dict_region = get_region_info(region_type, TEST=False)

main_figure_directory = '/data/rg001/sgec0017/data/figure_box/drought_tpcd_isimip2b'
input_directory = os.path.join(
    main_figure_directory, 'TPCD.estimator.scns.v2_for_1stRevision',
    region_type, 'all_all', 'Q80win15Len30tau4DRY',
    'basic.abs.05.Mea.1865-2005.max.05.historical_histsoc_co2_1861-2005.1865-2099.plt1865-',
    str(input_date)
    )
mapmask_directory = '/data/rg001/sgec0017/data/mapmask'
region_shp_path = os.path.join(mapmask_directory, 'HydroSHEDS', 'HydroBASINS', 
                               'withLakes', 'Global_merge', 'hybas_lake____lev02_v1c_merge')
shape_feature = ShapelyFeature(Reader(region_shp_path).geometries(), ccrs.PlateCarree(), 
                               edgecolor='#808080', linewidth=0.2, facecolor='none')
figure_directory = os.path.join(
    main_figure_directory, 'estimate_effective_non-corr_time',
    region_type
    )
if not os.path.isdir(figure_directory): os.makedirs(figure_directory)
#bm = Basemap(projection='cyl', llcrnrlat=-56.5, urcrnrlat=84.5, llcrnrlon=-180., urcrnrlon=180., resolution='l')
dpi = 300


def calc_te(src):
    lag_sizes = range(1, 95 + 1)
    for lag_size in lag_sizes:
        r1 = np.corrcoef(src[lag_size:], src[:-lag_size])[0, 1]
        r1 = abs(r1)
        if r1 <= 0.3:
            break
    return r1, lag_size


def calc_r1_at_lag5(src):
    return abs(np.corrcoef(src[5:], src[:-5])[0, 1])

def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])
    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value


def main():
    
    for rcp in rcps:

        # load input
        excel_file = f'TPCDs.nDayTot.005_max.Q80win15Len30tau4DRY.all.Mean.{rcp}_2005soc_co2.xlsx'
        excel_path = os.path.join(input_directory, excel_file)
        df = pd.read_excel(excel_path, sheet_name='data', index_col=0)
        indices = list(df.index)

        # main
        r1s_org_full, tes_org_full = [], []
        r1s_anmly_full, tes_anmly_full = [], []
        dict_regional_r1s = {}
        for region_id, region_name in dict_region.items():
            #print(f'\n>>> process for {region_id} {region_name}')

            fig = plt.figure(figsize=(26, 16))
            fig.suptitle(f'{region_id:02} {region_name}')
            plt.subplots_adjust(left=0.03, bottom=0.05, right=0.99, top=0.95, wspace=0.075)

            gs = gridspec.GridSpec(3, 2)
            ax1 = plt.subplot(gs[0:2, 0])
            ax2 = plt.subplot(gs[2, 0])
            ax3 = plt.subplot(gs[0, 1])
            ax4 = plt.subplot(gs[1, 1])
            ax5 = plt.subplot(gs[2, 1])

            r1s_org_rgn, tes_org_rgn = [], []
            r1s_anmly_rgn, tes_anmly_rgn = [], []
            dict_regional_r1s[f'{region_id:02}.{region_name}'] = []

            for ghm, gcm in itertools.product(ghms, gcms):
                index_name = f'{ghm}_{gcm}_{region_id:02}.{region_name}'
                idx = indices.index(index_name)
                src_org = df.iloc[idx, 30:]  # pandas Sries
                years = src_org.index.astype('int').values
                src_org = src_org.astype('float').values
                # --- org
                r1, te = calc_te(src_org)
                r1s_org_full.append(r1); r1s_org_rgn.append(r1)
                tes_org_full.append(te); tes_org_rgn.append(te)
                label = f'({r1:.3f}) {gcm}_{ghm} '
                ax1.plot(years, src_org, label=label)
                # --- liner and anomaly
                a, b, c = np.polyfit(years, src_org, deg=2)
                src_liner = a*years**2 + b*years + c
                src_anomaly = src_org - src_liner
                r1, te = calc_te(src_anomaly)
                #print(f'{region_id:02} {region_name:18} {ghm:<9} {gcm:<12}  >>>  a: {a:>7.4f}, b: {b:>5.1f}, c: {c:>7.1f}, te: {te:>2}, autcorr: {r1:.3f}')
                r1s_anmly_full.append(r1); r1s_anmly_rgn.append(r1)
                tes_anmly_full.append(te); tes_anmly_rgn.append(te)
                dict_regional_r1s[f'{region_id:02}.{region_name}'].append(calc_r1_at_lag5(src_anomaly))
                label = f'{gcm}_{ghm}'
                ax3.plot(years, src_liner, label=label)
                label = f'({r1:.3f}) {gcm}_{ghm}'
                ax4.plot(years, src_anomaly, label=label)

            ax1.set_xlim([2000,2130])
            ax1.legend()
            ax3.set_xlim([2000,2130])
            ax3.legend()
            ax4.set_xlim([2000,2130])
            ax4.legend()

            ax2.scatter(r1s_org_rgn, tes_org_rgn)
            ax2.set_xlabel('r1 (lag=1 auto_corr')
            ax2.set_ylabel('te (effective non-corr time)')

            ax5.scatter(r1s_anmly_rgn, tes_anmly_rgn)
            ax5.set_xlabel('r1 (lag=1 auto_corr')
            ax5.set_ylabel('te (effective non-corr time)')

            figure_name = f'{rcp}_{region_id:02}.{region_name}.pdf'
            figure_path = os.path.join(figure_directory, figure_name)
            plt.savefig(figure_path)
            plt.close()
            print(f'savefig: {figure_path}\n')
            #sys.exit()

        #fig = plt.figure()
        #fig.suptitle('original time series')
        #ax = fig.add_subplot()
        #ax.scatter(r1s_org_full, tes_org_full)
        #ax.set_xlabel('r1 (lag=1 auto_corr')
        #ax.set_ylabel('te (effective non-corr time)')
        #figure_name = 'scatter_full_org.pdf'
        #figure_path = os.path.join(figure_directory, figure_name)
        #plt.savefig(figure_path)
        #plt.close()
        #print('savefig: {}'.format(figure_path))

        # anomaly scatter
        fig = plt.figure(dpi=dpi)
        fig.suptitle('processed time series (trend removed)')
        ax = fig.add_subplot()
        ax.scatter(r1s_anmly_full, tes_anmly_full)
        ax.set_xlabel('r1 (lag=1 auto_corr')
        ax.set_ylabel('te (effective non-corr time)')
        figure_name = f'{rcp}_scatter_full_anomaly.pdf'
        figure_path = os.path.join(figure_directory, figure_name)
        plt.savefig(figure_path, dpi=dpi)
        plt.close()
        print(f'savefig: {figure_path}\n')

        # anomaly histogram
        edges = np.arange(min(tes_anmly_full), max(tes_anmly_full)+2) - 0.5
        fig = plt.figure(dpi=dpi)
        fig.suptitle('processed time series (trend removed)')
        ax = fig.add_subplot()
        ax.grid(visible=True, which='major', axis='y', linewidth=0.2, color='gray')
        ax.hist(tes_anmly_full, bins=edges,
                cumulative=True, density=True,
                color='skyblue', edgecolor='black', linewidth=0.5,
                )
        ax.set_ylim([0, 1])
        ax.yaxis.set_major_formatter(PercentFormatter(1))
        ax.set_xlim([0,edges[-1]])
        ax.set_xlabel('te (effective non-corr time)')
        figure_name = f'{rcp}_histogram_full_anomaly.pdf'
        figure_path = os.path.join(figure_directory, figure_name)
        plt.savefig(figure_path, dpi=dpi)
        plt.close()
        print(f'savefig: {figure_path}\n')

        # r1 boxplot
        r1s = [dict_regional_r1s[f'{region_id:02}.{region_name}'] for region_id, region_name in dict_region.items()]
        regions = [f'{region_id}  {region_name}' for region_id, region_name in dict_region.items()]
        #figsize, left, right, bottom, top = (12, 13), 0.18, 0.98, 0.03, 0.97  # horizontal
        figsize, left, right, bottom, top = (12, 12), 0.05, 0.98, 0.165, 0.97  # vertical
        #is_vert = False; r1s = r1s[::-1]; regions = regions[::-1]
        is_vert = True
        # ---
        fig = plt.figure(figsize=figsize, dpi=dpi)
        plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top)
        ax = fig.add_subplot()
        for spine in ['right', 'top', 'bottom']:
            ax.spines[spine].set_visible(False)
        ax.boxplot(r1s, 
                   vert=is_vert, 
                   showfliers=False,
                   patch_artist=True,
                   boxprops=dict(facecolor='skyblue',  # color of box
                                 color='black', linewidth=0.2),  # color of box frame
                   medianprops=dict(color='black', linewidth=1.5),  # color of median line
                   whiskerprops=dict(color='gray', linewidth=0.5),  # color of winsker
                   capprops=dict(color='gray', linewidth=0.9),  # color of cap
                   )
        if is_vert:
            ax.grid(visible=True, which='major', axis='y', linewidth=0.6, color='black', alpha=0.3)
            #xticks = range(len(regions))
            #ax.set_xticks(xticks)
            ax.set_xticklabels(regions, rotation=75, ha='right')
            yticks = [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.]
            ax.set_yticks(yticks)
            ax.set_yticklabels(yticks, fontsize=12)
            ax.set_ylim([0,1.2])
        else:
            ax.grid(visible=True, which='major', axis='x', linewidth=0.2, color='black', alpha=0.3)
            yticks = range(len(regions))[::-1]
            ax.set_yticks(yticks)
            ax.set_yticklabels(regions)
            xticks = [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.]
            ax.set_xticks(xticks)
            ax.set_xticklabels(xticks)
            ax.set_xlim([0,1])
        figure_name = f'r1_at_lag5_boxplot_{rcp}_full_anomaly.pdf'
        figure_path = os.path.join(figure_directory, figure_name)
        plt.savefig(figure_path, dpi=dpi)
        plt.close()
        print(f'savefig: {figure_path}\n')

        # median r1 map
        missing_value = 1e+20
        median_r1_map = np.full(shape=(360, 720), fill_value=missing_value)
        for region_id, region_name in dict_region.items():
            median_r1 = np.median(dict_regional_r1s[f'{region_id:02}.{region_name}'])
            YY, XX = np.where(region_map==region_id)
            median_r1_map[YY,XX] = median_r1
            #print(f'@{region_id:02}: {median_r1}')
        #median_r1_map = median_r1_map[11:11+280, ...]
        median_r1_map = median_r1_map[:300, ...]
        median_r1_map = np.ma.masked_equal(median_r1_map, missing_value)
        # ---
        figsize, left, right, bottom, top = (6, 2.5), 0.001, 0.999, 0.01, 0.99
        vmin, vmax = 0, 0.6
        norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)
        cmap = mpl.cm.terrain
        extent = [-180., 180., -60., 90.]
        #proj_fig = ccrs.Robinson()
        proj_fig = ccrs.PlateCarree()
        proj_org = ccrs.PlateCarree()
        # ---
        fig = plt.figure(figsize=figsize)
        #fig.suptitle('median r1')
        plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top)
        #ax = fig.add_subplot(111)
        ax = fig.add_subplot(1,1,1, projection=proj_fig)
        ax.set_extent(extent, proj_fig)
        ax.outline_patch.set_linewidth(0.3)
        #ax.axis('off')
        ax.imshow(median_r1_map, 
                  origin='upper', extent=extent, transform=proj_org, 
                  #zorder=9,
                  vmin=vmin, vmax=vmax, cmap=cmap)
        #bm.imshow(median_r1_map[::-1], vmin=vmin, vmax=vmax, cmap=cmap)
        #bm.drawcoastlines(linewidth=0.05, color='#808080')
        ax.add_feature(cfeature.COASTLINE, linewidth=0.2, edgecolor='#000000')
        ax.add_feature(shape_feature)
        ax_pos = ax.get_position()
        #ax11 = fig.add_axes([ax_pos.x0+0.01, ax_pos.y0+0.05, 0.25, 0.02])
        #ax11 = fig.add_axes([ax_pos.x0+0.015, ax_pos.y0+0.075, 0.25, 0.02])
        ax11 = fig.add_axes([ax_pos.x0+0.025, ax_pos.y0+0.2, 0.25, 0.02])
        cb1 = mpl.colorbar.ColorbarBase(ax11, cmap=cmap, norm=norm, orientation='horizontal')
        cb1.outline.set_visible(False)
        cb1.set_label('median', color='black')
        #cb1_ticks = cb1.get_ticks()
        cb1_ticks = [0., 0.3, 0.5]
        cb1.set_ticks(cb1_ticks)
        cb1.set_ticklabels(['0'if v==0. else str(v) for v in cb1_ticks])
        figure_name = f'r1_median_at_lag5_map_{rcp}_full_anomaly.pdf'
        figure_path = os.path.join(figure_directory, figure_name)
        plt.savefig(figure_path, dpi=dpi)
        plt.close()
        print(f'savefig: {figure_path}\n')


if __name__ == '__main__':
    main()
