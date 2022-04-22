import os
import sys
import itertools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from netCDF4 import Dataset
from utiltools import get_region_info
plt.rcParams["font.size"] = 8


region_type = 'HydroBASINS_lev1'
input_date = 20220102
regions_selected = [
    8, 11, 13, 14, 16,
    17, 18, 19, 20, 22,
    23, 26, 29, 31, 32,
    35, 38, 51, 57
    ]

ghms = ['cwatm', 'h08', 'lpjml', 'matsiro', 'watergap2']  # full members
gcms = ['hadgem2-es',  'ipsl-cm5a-lr', 'gfdl-esm2m', 'miroc5']  # full members
rcps = ['rcp26', 'rcp85']

main_figure_directory = '/data/rg001/sgec0017/data/figure_box/drought_tpcd_isimip2b'
input_directory = os.path.join(
    main_figure_directory, 'TPCD.estimator.scns.v2_for_1stRevision',
    region_type, 'all_all', 'Q80win15Len30tau4DRY',
    'basic.abs.05.Mea.1865-2005.max.05.historical_histsoc_co2_1861-2005.1865-2099.plt1865-',
    str(input_date)
    )
figure_directory = os.path.join(
    main_figure_directory, 'find_the_best_fittingcurve.v2',
    region_type
    )
if not os.path.isdir(figure_directory): os.makedirs(figure_directory)

region_map, dict_region = get_region_info(region_type, TEST=False)
regions = [f'{region_id:02}.{region_name}' for region_id, region_name in dict_region.items()]
del region_map

syear, eyear = 2005, 2099
years = np.arange(syear, eyear+1)
xticks = [2010, 2030, 2050, 2070, 2090]
suffixs = ['png', 'pdf']


def load_df_from_xlsx(rcp):
    excel_file = f'TPCDs.nDayTot.005_max.Q80win15Len30tau4DRY.all.Mean.{rcp}_2005soc_co2.xlsx'
    excel_path = os.path.join(input_directory, excel_file)
    return pd.read_excel(excel_path, sheet_name='data', index_col=0)


def extract_target_timeseries(df, ghm, gcm, region):
    index = '{}_{}_{}'.format(ghm, gcm, region)
    indices = df.index
    columns = df.columns
    return df.iloc[list(indices).index(index), list(columns).index(syear):].astype('float').values


def main():

    for rcp in rcps:

        df = load_df_from_xlsx(rcp)
        print(df)

        for (ighm, ghm), (igcm, gcm) in itertools.product(enumerate(ghms), enumerate(gcms)):

            fig = plt.figure(figsize=(14, 10), dpi=300)
            fig.suptitle(f'{rcp}  {ghm}-{gcm}')
            #plt.subplots_adjust(left=0.05, right=0.985, bottom=0.05, top=0.88, wspace=0.15, hspace=0.25)
            plt.subplots_adjust(left=0.05, right=0.985, bottom=0.05, top=0.94, wspace=0.15, hspace=0.3)
            nrow = len(regions_selected)//5+1
            ncol = 5
            gs = gridspec.GridSpec(nrow, ncol)

            counter = 0
            for region in regions:
                
                if int(region[:2]) in regions_selected:
                    print(f'{rcp} {ighm}.{ghm} {igcm}.{gcm} @{region}')

                    irow, icol = divmod(counter, ncol)
                    ax = plt.subplot(gs[irow, icol])
                    ax.set_title(f'{region[:2]}  {region[3:]}')

                    # src
                    timeseries_org = extract_target_timeseries(df, ghm, gcm, region)

                    # original time series
                    ax.plot(years, timeseries_org, color='#000000', linewidth=0.8, linestyle='--', label='original')
                    # moving average
                    N = 31
                    window_size_half = int((N-1)/2)
                    weights = np.ones(N)/N
                    moving_average = np.convolve(timeseries_org, weights, 'valid')
                    ax.plot(years[window_size_half:-window_size_half], moving_average, color='#000000', linewidth=1.5, linestyle='-', label='moving avr ({}yr)'.format(N))
                    # fitting plots            
                    for deg in range(1,5+1):
                        coefficients = np.polyfit(years, timeseries_org, deg=deg)
                        timeseries_fitted = np.ones(years.shape)
                        for icf, cf in enumerate(coefficients):
                            timeseries_fitted += cf * years ** (len(coefficients)-icf-1)
                        ax.plot(years, timeseries_fitted, linewidth=1.2, linestyle='-', label='deg={}'.format(deg)) 
                    # ylabel
                    if icol == 0:
                        ax.set_ylabel('regional average FDD', fontsize=10)
                    # legend
                    if counter == len(regions_selected)-1:
                        ax.legend(bbox_to_anchor=(1.15, 1.), loc='upper left', borderaxespad=0, fontsize=10)
                    #if ighm == 0 and igcm == 0:
                    #    ax.legend(bbox_to_anchor=(0.01, 1.09), loc='lower left', borderaxespad=0, fontsize=8)

                    counter += 1

            for suffix in suffixs:
                fig_name = f'{rcp}.{ghm}.{gcm}.{suffix}'
                fig_path = os.path.join(figure_directory, fig_name)
                plt.savefig(fig_path, dpi=300)
                print(f'savefig: {fig_path}')
            plt.close()

            #sys.exit()


if __name__=='__main__':
    main()
