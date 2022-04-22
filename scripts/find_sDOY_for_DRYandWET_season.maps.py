#!/usr/bin/env python
# To
# By Yusuke Satoh
# On 20210623
import os
import itertools
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

#timing = 'beginning'
timing = 'middle'

INDEVIDUAL = False

seasons = ['DRY', 'WET']
ghms = ['cwatm', 'h08', 'lpjml','matsiro', 'watergap2']
gcms = ['hadgem2-es','ipsl-cm5a-lr','gfdl-esm2m','miroc5']
suffixes = ['png', 'pdf']

data_directory = '/data/rg001/sgec0017/data/figure_box/drought_tpcd_isimip2b/find_sDOY_for_DRYandWET_season'

def load_input(season, ghm, gcm):
    filename = 'startDOY_{}season_{}_{}_historical_histsoc_co2.nc4'.format(season, ghm, gcm)
    path_src = os.path.join(data_directory, filename)
    return Dataset(path_src)['sDOY'][:]


def draw_map(src, season, ghm='ensemble', gcm='ensemble', num=0):

    if timing == 'middle': 
        src = src + (91/2)
        src = np.where(src>365, src-365, src)
    src = np.ma.masked_greater_equal(src, 1e+19)
    print('range: {}-{}'.format(src.min(), src.max()))

    fig = plt.figure(figsize=(8, 5))
    plt.subplots_adjust(left=0.02, bottom=0.06, right=0.98, top=0.90)
    src = np.ma.masked_greater_equal(src, 1e+19)
    ax = fig.add_subplot(111)
    ax.axis('off')
    im = ax.imshow(src[:300], cmap=plt.get_cmap('twilight'))
    ax.set_title('{} season (3months)\n {} {}'.format(season, ghm, gcm))
    clb = plt.colorbar(im, orientation='horizontal', pad=0, aspect=45)
    clb.set_label('DOY of the {} of the season'.format(timing))
    clb.outline.set_visible(False)

    for suffix in suffixes:
        fig_directory = os.path.join(data_directory, 'season_map', season, suffix)
        if not os.path.isdir(fig_directory): os.makedirs(fig_directory)
        figname = '{}.{}.{}.{}.{}.{}'.format(season, timing, num, ghm, gcm, suffix)
        path_fig = os.path.join(fig_directory, figname)
        plt.savefig(path_fig)
        print('savefig: {}'.format(path_fig))
    plt.close()


def main():

    for season in seasons:
        src = np.array([load_input(season, ghm, gcm) for ghm in ghms for gcm in gcms])
        src = np.median(src, axis=0)
        draw_map(src, season)

    if INDEVIDUAL:
        for season, ghm, gcm in itertools.product(seasons, ghms, gcms):
            src = load_input(season, ghm, gcm)
            draw_map(src, season, ghm, gcm, 1)

if __name__ == '__main__':
    main()
