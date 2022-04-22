#!/usr/bin/python
# by Yusuke Satoh (NIES->KAIST)
import sys
import os
import itertools
import numpy as np
import matplotlib.pyplot as plt
import numpy.ma as ma
import pandas as pd
import matplotlib as mpl
import datetime
from numpy import fromfile, array, zeros, where
from numpy import divide, median, mean, subtract, percentile, var
from netCDF4 import Dataset 
from matplotlib import cm, colors, rcParams, gridspec
from utiltools import fillinMasked, extractMasked, flipud
from scipy import stats
#------------------------------------------------------------------------------------------------------
plt.rcParams['font.family'] = 'Arial'
hostname = os.uname()[1]

TEST = False

#SSNs = ['DJF','MAM','JJA','SON']
#SSNs = ['ALL']
#SSNs = ['DRY', 'WET']
SSNs = ['DRY', 'ALL', 'WET']
#SSNs = ['ALL', 'DRY', 'WET']
#SSNs = ['WET']
#SSNs = ['DRY']
#SSNs = ['ALL', 'WET']

#mmeTYPs = ['median','mean']  
mmeTYPs = ['median']  
#mmeTYPs = ['mean']  

#IDXs = ['nDayTot', 'dfcTot', 'nEvent', 'nOnset', 'avrDSL', 'maxDSL']
IDXs = ['nDayTot']

ghms = ['cwatm', 'h08', 'lpjml', 'matsiro', 'watergap2']
gcms = ['hadgem2-es', 'ipsl-cm5a-lr', 'gfdl-esm2m', 'miroc5']

# target simlation type
SIMs  = [
    ('rcp26', '2005soc', 'co2'),
    ('rcp85', '2005soc', 'co2'),
    #('picontrol', '2005soc', 'co2'),
    ]

#region_type = 'AR6_regions'
region_type = 'HydroBASINS_lev1'

# reference historical period of the drought detection
syear_ref, eyear_ref = 1861, 2005  # main
#syear_ref, eyear_ref = 1971, 2005
#syear_ref, eyear_ref = 1979, 2013  # test

# a whole target period of drought analysis
syear, eyear = syear_ref, 2099

# window of time slices
climate_window = 30

syear_hi = 1976
syear_mf = 2036
syear_ff = 2070
eyear_hi = syear_hi+climate_window-1
eyear_mf = syear_mf+climate_window-1
eyear_ff = syear_ff+climate_window-1
if syear_ref == 1979: syear_hi = 1979; eyear_hi = 2013  # test
futures = ['mid', 'late']

#agreeThrsh = [0.8, 0.7, 0.6]
agreeThrsh = [0.6]

projection = 'Basemap_PlateCarree'
#projection = 'Basemap_Robinson'
#projection = 'Cartopy_PlateCarree'
#projection = 'Cartopy_Robinson'

#suffixes = ['.png', '.eps', '.pdf']
suffixes = ['png', 'pdf']

# which figure do you want?
EnsmMap = True
AGREEMENT = True
S2NMap  = True
GHMMap  = False
UncertaintySource = False

# masks
DesertMASK = False
KStest = True
#KStest = False
ks = 0.05    # 95% level significance

# do you want to KS mask white? or gray?  i.e., do you want to set background color in gray?
KSmask_gray = True  # if True, gray. if False, white


# --- qvalType
#autQval = True                   # automatically specify, depending on sim type: 'Qvale_hist_{}'.format(sim)
autQval = False                  # set qvaltype, below
#qvalTYPE = 'Qvale_hist_nosoc'
#qvalTYPE = 'Qvale_hist_pressoc'
# --- Qvalue
#Qs = [90, 85, 80, 75, 70, 65]
#Qs = [90,80,70]
#Qs = [80,90,70]
#Qs = [80,90]
Qs = [80]
# --- window size for Qvalue
#WINs = [15, 10, 7, 5]
WINs = [15]
# --- minimum drought days
#LENs = [180,90,60,30,14,7]
#LENs = [30,14,60,90,180,7]
LENs = [30]
# --- Pool duration
TAUs = [4]

#coastal_line_color = '#808080'
#coastal_line_color = '#3c3c3c'
coastal_line_color = '#2c2c2c'

#coastal_line_width = 0.1
coastal_line_width = 0.2
#coastal_line_width = 0.3
#coastal_line_width = 0.5


if TEST:
    #ghms = ['matsiro']
    #gcms = ['hadgem2-es']
    ghms = ['cwatm', 'matsiro']
    gcms = ['hadgem2-es', 'ipsl-cm5a-lr']
    KStest = False
    SIMs = [('rcp85', '2005soc')]


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
today = datetime.date.today().strftime('%Y%m%d')

### set input data directory: 
if hostname == 'scs01' or hostname == 'sun220s.nies.go.jp':
    if   hostname == 'scs01':              data_directory = '/data/rg001/sgec0017/data'
    elif hostname == 'sun220s.nies.go.jp': data_directory = '/sraid02/satoh/data'
    elif hostname == 'sun90h.nies.go.jp':  data_directory = '/sraid02/satoh/data'
    drghtDIR = os.path.join(data_directory, 'isimip2b.drought')
    if KStest: fig_directory_main = os.path.join(data_directory, 'figure_box', 'drought_tpcd_isimip2b', 'basicMaps.RCPxFuture_withKStest',    region_type, projection)
    else:      fig_directory_main = os.path.join(data_directory, 'figure_box', 'drought_tpcd_isimip2b', 'basicMaps.RCPxFuture_withoutKStest', region_type, projection)
    if TEST: fig_directory_main = os.path.join(data_directory, 'figure_box', 'drought_tpcd_isimip2b', 'basicMaps.RCPxFuture_TEST', today)
    mapmaskDIR = os.path.join(data_directory, 'mapmask')
lndmskPath = os.path.join(mapmaskDIR, 'ISIMIP2b_landseamask_generic.nc4')
grlmskPath = os.path.join(mapmaskDIR, 'GAUL/flt/gaul2014_05deg.flt')      # GreenLand is 98
dstmskPath = os.path.join(mapmaskDIR, 'grid_index.nc4')
grdaraPath = os.path.join(mapmaskDIR, 'grd_ara.hlf')
if region_type == 'AR6_regions':
    ipcc_shp_path = os.path.join(mapmaskDIR, 'IPCC-AR6-WGI-reference-regions-v4_shapefile', 
                                 'shapefile_edited', '10regions_for_timeseries', 'reference-regions-AR6')
    shp_name = 'reference-regions-AR6'

gridarea = np.fromfile(grdaraPath, 'float32').byteswap().reshape(360,720)


dictUnit = {
    'nDayTot': '%',
    'dfcTot':  'm3',
    'nEvent':  'times',
    'nOnset':  'times',
    'avrDSL':  'days',
    'maxDSL':  'days',
    }

dictVlim = {
    'nDayTot': {'absolute': [0,1], 'abschange': [0,2.5], 'percentchange': [-10,200]},
    'dfcTot':  {'absolute': [None, None], 'abschange': [None, None], 'percentchange': [None, None]},
    'nEvent':  {'absolute': [None, None], 'abschange': [None, None], 'percentchange': [None, None]},
    'nOnset':  {'absolute': [None, None], 'abschange': [None, None], 'percentchange': [None, None]},
    'avrDSL':  {'absolute': [None, None], 'abschange': [None, None], 'percentchange': [None, None]},
    'maxDSL':  {'absolute': [None, None], 'abschange': [None, None], 'percentchange': [None, None]},
    }

dictSeason = {
    'ALL': 365,
    'DJF':  90,
    'MAM':  92,
    'JJA':  92,
    'SON':  91,
    'DRY':  91,
    'WET':  91,
    }

years = range(syear, eyear+1)
nyears = len(years)
nyears_hist = len(range(syear_hi, eyear_hi+1))
nyears_future = climate_window
histStr, histEnd = years.index(syear_hi), years.index(eyear_hi)  # 1976-2005 (30)
nfStr, nfEnd     = years.index(syear_mf), years.index(eyear_mf)  # 2036-2065 (30)
ffStr, ffEnd     = years.index(syear_ff), years.index(eyear_ff)  # 2070-2099 (30)
print('\n--------------------------------------------')
print('full            : {}-{} (n={})'.format(syear, eyear, nyears))
print('historical (ref): {}-{} (n={})'.format(syear_ref, eyear_ref, len(range(syear_ref, eyear_ref+1))))
print('hist            : {}-{} (n={}, #{}-{})'.format(syear_hi, eyear_hi, nyears_hist,    years.index(syear_hi), years.index(eyear_hi)))
print('mid             : {}-{} (n={}, #{}-{})'.format(syear_mf, eyear_mf, climate_window, years.index(syear_mf), years.index(eyear_mf)))
print('late            : {}-{} (n={}, #{}-{})'.format(syear_ff, eyear_ff, climate_window, years.index(syear_ff), years.index(eyear_ff)))
print('--------------------------------------------')
print('projection = {}\n'.format(projection))

nGHM, nGCM = len(ghms), len(gcms)
nMME = nGHM*nGCM

grl_mask = ma.masked_equal(fromfile(grlmskPath,'float32').reshape(360,720), 98).mask
dst_mask = ma.masked_equal(Dataset(dstmskPath).variables['grid_index'][:], 11).mask
lnd_mask = ma.masked_not_equal(Dataset(lndmskPath).variables['LSM'][:][0], 1.0).mask
lnd_mask = ma.mask_or(lnd_mask, grl_mask)
if DesertMASK: lnd_mask = ma.mask_or(lnd_mask, dst_mask)
nlndgrid = 360*720 - lnd_mask.sum()
total_gridarea = ma.masked_array(gridarea, mask=Dataset(lndmskPath)['LSM'][0].mask).sum()

if projection == 'Basemap_PlateCarree':
    from mpl_toolkits.basemap import Basemap
    bm1 = Basemap(projection='cyl',llcrnrlat=-56.5,urcrnrlat=84.5,llcrnrlon=-180.,urcrnrlon=180.,resolution='l') # correct
elif projection == 'Basemap_Robinson':  # this does not work... coordinate won't be changed...
    from mpl_toolkits.basemap import Basemap
    bm1 = Basemap(projection='robin', lon_0=0, lat_0=0, resolution='l') # correct
    lons, lats = np.arange(-179.75, 179.75+0.5, 0.5), np.arange(-89.75, 89.75+0.5, 0.5)
    size = 10000.
    nx = int((bm1.xmax-bm1.xmin)/size)+1
    ny = int((bm1.ymax-bm1.ymin)/size)+1
elif 'Cartopy' in projection:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    from cartopy.io.shapereader import Reader
    from cartopy.feature import ShapelyFeature
    original_projection = ccrs.PlateCarree()
    if projection == 'Cartopy_PlateCarree': axis_projection = ccrs.PlateCarree()
    elif projection == 'Cartopy_Robinson':  axis_projection = ccrs.Robinson()
    lons, lats = np.arange(-179.75, 179.75+0.5, 0.5), np.arange(89.75, -56.25-0.5, -0.5)

dpi = 300
missing_value = 1e+20


# ----------------------------------------------------------------------------------------------------------------------
def read_nc(index, soc, scn, co2, qvalType, win, Q, Len, tau, ghm, gcm, season):
    # Read only Land Data

    srcFile = 'Q{:02}win{:02}_Len{:03}tau{}_{}_{}.nc4'.format(Q,win,Len,tau,season,index)
    srcDir = os.path.join(drghtDIR, ghm, gcm,
                          '{}.{}_{}'.format(qvalType, syear_ref, eyear_ref),
                          'droughtStats.{}_{}_{}.{}_{}'.format(scn, soc, co2, syear, eyear)
                          )
    srcPath = os.path.join(srcDir,srcFile)

    if not os.path.isfile(srcPath): 
        print('Error!! The file NOT exist... Check!! {}'.format(srcPath))
        sys.exit()
    else:
        aSrc = extractMasked(Dataset(srcPath).variables[index][:], lnd_mask)
        print('Read: {} {}'.format(srcPath, aSrc.shape))
        return aSrc


# ----------------------------------------------------------------------------------------------------------------------
def change_arearate(target_type, src, intensity=None):
    """
    To estimate the rate of certain change (decrease/increase)
    :param target_type: (str) decrease, increase, significant increase
    :param src: (array) percentchange or absolutechange
    :return: gridrate
    """

    if target_type == 'decrease':               target_mask = ma.make_mask(0<=src)
    elif target_type == 'increase':             target_mask = ma.make_mask(src<=0)
    elif target_type == 'significant_increase': target_mask = ma.make_mask(src<=intensity)
    else: print('Something is wrong. Check target_type!!')

    target_gridarea = ma.masked_array(gridarea, mask=target_mask).sum()
    print('change_area_rate: ({}) {} / {} = {}'.format(target_type, target_gridarea, total_gridarea, target_gridarea / total_gridarea))

    return target_gridarea / total_gridarea


# ----------------------------------------------------------------------------------------------------------------------
def main(*args):
    strTIME = datetime.datetime.now()
    print('START detect.drought.py')

    for win, Len, tau, Q, season, mmeType, index in itertools.product(WINs, LENs, TAUs, Qs, SSNs, mmeTYPs, IDXs):

        strTime = datetime.datetime.now()
        print('\nJob started !!!  {} with Q{} win{} for Len{} tau{}'.format(season, Q, win, Len, tau))

        qvalType = 'Qvalue_historical_histsoc_co2'
        #if autQval: qvalType = 'Qvalue_historical_{}'.format(soc)
        #else:
        #    if soc == 'nosoc': qvalType = 'Qvalue_historical_nosoc'
        #    else:              qvalType = 'Qvalue_historical_histsoc'
        #print('autQval : {}'.format(autQval))
        #print('qvalType: {}\n'.format(qvalType))

        # Reading data
        aSRC = array([[[read_nc(index,soc,scn,co2,qvalType,win,Q,Len,tau,ghm,gcm,season) for ghm in ghms] 
                                                                                         for gcm in gcms] 
                                                                                         for (scn, soc, co2) in SIMs])
        aSRC = ma.masked_equal(aSRC,1e+20)             # orignally, missing_value in each data is 1e+20
        print('aSRC.shape {}'.format(aSRC.shape))      #(nSCN,nGCM,nGHM,nYear,nLand))

        # Cutout specific periods
        aHIS = aSRC[:,:,:,histStr:histEnd+1,:]         #      (nSIM,nGCM,nGHM,30,nLand)
        aFTR = array([aSRC[:,:,:,nfStr:nfEnd+1,:],
                      aSRC[:,:,:,ffStr:ffEnd+1,:]])    # (nPRD,nSIM,nGCM,nGHM,30,nLand)
        del aSRC

        ### Kolmogorov-Smirnov test and dictKSmask
        ### if p-value is less than ks, change in the grid is significant.
        print('KS testing...')
        dictKSmask = {scn: {'mid': {}, 'late': {}} for scn,soc,co2 in SIMs}
        if KStest:
            for s, (scn, soc, co2) in enumerate(SIMs):
                for i, f in enumerate(futures):
                    dictKSmask[scn][f]['all'] = ma.make_mask(
                                                    fillinMasked(
                                                        array([ 0 if stats.ks_2samp(aHis.reshape(-1),aFut.reshape(-1)).pvalue < ks else 1
                                                                for aHis, aFut in zip(aHIS[s].T, aFTR[i,s].T)]),
                                                        lnd_mask)
                                                    == 1) #(nY,nX)
                    print('{} {} all {}  {}'.format(scn, f, dictKSmask[scn][f]['all'].shape, 
                                                    (nlndgrid-dictKSmask[scn][f]['all'].sum())/float(nlndgrid)))
                    for j, ghm in enumerate(ghms):
                        dictKSmask[scn][f][ghm] = ma.make_mask(
                                                    fillinMasked(
                                                        array([ 0 if stats.ks_2samp(aHis.reshape(-1),aFut.reshape(-1)).pvalue < ks else 1
                                                                for aHis, aFut in zip(aHIS[s,:,j,...].T, aFTR[i,s,:,j,...].T)]),
                                                        lnd_mask)
                                                    == 1) #(nY,nX)
                        print('{} {} {} {}  {}'.format(scn, f, ghm, dictKSmask[scn][f][ghm].shape, 
                                                       (nlndgrid-dictKSmask[scn][f][ghm].sum())/float(nlndgrid)))
        else:
            print('KStest is False. Skep it.')
            for i, f in enumerate(futures):
                for s, (scn, soc, co2) in enumerate(SIMs):
                    dictKSmask[scn][f]['all'] = zeros((360,720),'bool')
                    for j, ghm in enumerate(ghms):
                        dictKSmask[scn][f][ghm] = zeros((360,720),'bool')

        print('genarating values...')
        # Make ensemble value for each combination of GHM&GCM
        if index == 'nDayTot':
            clmHIS = aHIS.sum(axis=3)                            # (     nSIM,nGCM,nGHM,nLand)
            clmFTR = aFTR.sum(axis=4)                            # (nPRD,nSIM,nGCM,nGHM,nLand)
            # unit conv [day] >> [%]
            clmHIS = array([ divide(clmHIS[iscn],             nyears_hist*dictSeason[season]) for iscn, (scn, soc, co2) in enumerate(SIMs)]) * 100
            clmFTR = array([[divide(clmFTR[ifuture, iscn], climate_window*dictSeason[season]) for iscn, (scn, soc, co2) in enumerate(SIMs)]
                                                                                              for ifuture, future  in enumerate(futures)]) * 100
        else:
            clmHIS = aHIS.mean(axis=3)                           # (     nSIM,nGCM,nGHM,nLand)
            clmFTR = aFTR.mean(axis=4)                           # (nPRD,nSIM,nGCM,nGHM,nLand)
        pc_clmFTR= (clmFTR-clmHIS)/clmHIS*100                    # (nPRD,nSIM,nGCM,nGHM,nLand)          include inf
        ##pc_clmFTR[np.isinf(pc_clmFTR)] = 1e+20
        ##pc_clmFTR = np.ma.masked_equal(pc_clmFTR, 1e+20)         # inf is masked out

        if mmeType == 'median':
            allMMEhis = median(clmHIS, axis=(1,2))               # (     nSIM,         nLand) >> Fig
            gcmMMEhis = median(clmHIS, axis=(1))                 # (     nSIM,    nGHM,nLand) >> Fig
            allMMEftr = median(clmFTR, axis=(2,3))               # (nPRD,nSIM,         nLand) >> Fig(sub)
            gcmMMEftr = median(clmFTR, axis=(2))                 # (nPRD,nSIM,    nGHM,nLand) >> Fig(sub)
        elif mmeType == 'mean':
            allMMEhis = mean(clmHIS, axis=(1,2))                 # (     nSIM,         nLand) >> Fig
            gcmMMEhis = mean(clmHIS, axis=(1))                   # (     nSIM,    nGHM,nLand) >> Fig
            allMMEftr = mean(clmFTR, axis=(2,3))                 # (nPRD,nSIM,         nLand) >> Fig(sub)
            gcmMMEftr = mean(clmFTR, axis=(2))                   # (nPRD,nSIM,    nGHM,nLand) >> Fig(sub)
        else:
            print('Warning. Check mmeType.'); sys.exit()

        # Change (chng)
        chng_allMME = allMMEftr - allMMEhis                      # (nPRD,nSIM,         nLand) >> Fig
        chng_gcmMME = gcmMMEftr - gcmMMEhis                      # (nPRD,nSIM,    nGHM,nLand) >> Fig
        # Percentage change (pc) [%]
        pc_allMME = chng_allMME / allMMEhis * 100             # (nPRD,nSIM,         nLand) >> Fig [%]   (include inf)
        pc_gcmMME = chng_gcmMME / gcmMMEhis * 100             # (nPRD,nSIM,    nGHM,nLand) >> Fig [%]   (include inf)

        ### Spread(spr) among ensemble samples
        if mmeType == 'median':   ### get inter quartile range (IQR)
            spr_allMMEhis = subtract(*percentile(clmHIS,    [75,25], axis=(1,2)))  #      (nSIM,nLand)
            spr_pc_allMME = subtract(*percentile(pc_clmFTR, [75,25], axis=(2,3)))  # (nPRD,nSIM,nLand)  (include inf??)
        elif mmeType == 'mean':   # get standard deviation (std)
            spr_allMMEhis = clmHIS.std(axis=(1,2))                                 #      (nSIM,nLand)
            spr_pc_allMME = pc_clmFTR.std(axis=(2,3))                              # (nPRD,nSIM,nLand)  (include inf??)

        """ 
        This part is to investigate why there are regions where S2N is empty.
        Still not sure why... leave this issue for a while...

        fig = plt.figure()

        ax1 = fig.add_subplot(421)
        ax1.set_title('allMMEftr[0,0]')
        im1 = ax1.imshow(np.ma.masked_equal(fillinMasked(allMMEftr[0,0], lnd_mask), 1e+20))
        plt.colorbar(im1)
       
        ax2 = fig.add_subplot(422)
        ax2.set_title('allMMEhis[0]')
        im2 = ax2.imshow(np.ma.masked_equal(fillinMasked(allMMEhis[0], lnd_mask), 1e+20))
        plt.colorbar(im2)
       
        ax3 = fig.add_subplot(423)
        ax3.set_title('pc_allMME[0,0]')
        im3 = ax3.imshow(np.ma.masked_equal(fillinMasked(pc_allMME[0,0], lnd_mask), 1e+20))
        plt.colorbar(im3)
       
        ax5 = fig.add_subplot(425)
        ax5.set_title('Q25')
        im5 = ax5.imshow(np.ma.masked_equal(fillinMasked(percentile(pc_clmFTR, 25, axis=(2,3))[0,0], lnd_mask), 1e+20))
        plt.colorbar(im5)

        ax6 = fig.add_subplot(426)
        ax6.set_title('Q75')
        im6 = ax6.imshow(np.ma.masked_equal(fillinMasked(percentile(pc_clmFTR, 75, axis=(2,3))[0,0], lnd_mask), 1e+20))
        plt.colorbar(im6)

        ax7 = fig.add_subplot(427)
        ax7.set_title('spr_pc_allMME[0,0]')
        im7 = ax7.imshow(np.ma.masked_equal(fillinMasked(spr_pc_allMME[0,0], lnd_mask), 1e+20))
        plt.colorbar(im5)

        plt.show()
        sys.exit()
        """

        # Uncertainty  (Signal to noise)  TODO!!
        if mmeType == 'median':   ### Singal to noise (s2n)   (ref. Prudhomme et al. 2014)
            s2n_allMMEhis = allMMEhis / spr_allMMEhis            #      (nSIM,nLand) >> Fig             (include inf&nan)
            s2n_pc_allMME = np.abs(pc_allMME) / spr_pc_allMME            # (nPRD,nSIM,nLand) >> Fig     (include inf&nan)
        elif mmeType == 'mean':   ### Coefficient of variation (Cov)
            s2n_allMMEhis = spr_allMMEhis / allMMEhis            #      (nSIM,nLand) >> Fig             (include inf&nan)
            s2n_pc_allMME = np.abs(spr_pc_allMME) / pc_allMME            # (nPRD,nSIM,nLand) >> Fig     (include inf&nan)

        """
        [Boulange Julien et al 2018 Environ. Res. Lett.]
        The SN ratios deﬁned by the IPCC Third Assessment Report(IPCC 2013), were calculated individually for each grid cell, 
        and for the entire sequence of monthly values within a given 30 year period:
            𝑆𝑁 = | Δ / 𝜎 |
        where 
          SN is the signal to noise ratio of a speciﬁc variable, 
          Δ is the difference between the mean value of a variable for the experimental and control periods,
          𝜎 is the standard deviation over the control period of the variable (Hawkins and Sutton 2011). 
        The SN ratio represents the strength of the change signal compared to natural variability(noise), 
        and the signal stands out against the noise when and where this ratio is large(IPCC 2013).
        """

        # Uncertainty comparison: GHM vs GCM  TODO!!
        # Ratio of GCM variation to total variation   (ref. Schewe et al. 2014)
        # GCM variation was computed across all gcms for each GHM individually and then averaged over all ghms
        ratGCMvar_clmHIS    = var(clmHIS,axis=1).mean(axis=1) / var(clmHIS,axis=(1,2))         #      (nSIM,nLand) >> Fig (include nan?) 
        ratGCMvar_pc_allMME = var(pc_clmFTR,axis=2).mean(axis=2) / var(pc_clmFTR,axis=(2,3))   # (nPRD,nSIM,nLand) >> Fig (include nan?) 
        del clmHIS, clmFTR
        del spr_allMMEhis, spr_pc_allMME

        # Agreement on the sign of change (increase/decrease)  (0-1)
        # This is only for allMME, for now
        # all member
        flug_allMME  =  where(pc_allMME>0,1,0) + where(pc_allMME<0,-1,0)                       # (          nPRD,nSIM,nLand)
        flug_eachMME = (where(pc_clmFTR>0,1,0) + where(pc_clmFTR<0,-1,0)).transpose(2,3,0,1,4) # (nGCM,nGHM,nPRD,nSIM,nLand)
        agreementALL =  where(flug_eachMME==flug_allMME,1,0).sum(axis=(0,1)) / float(nMME)     # (          nPRD,nSIM,nLand) >> Fig
        # GCM ansemble for each GHM
        flug_allMME  = where(pc_gcmMME>0,1,0) + where(pc_gcmMME<0,-1,0)                        #      (nPRD,nSIM,nGHM,nLand)
        flug_eachMME = (where(pc_clmFTR>0,1,0) + where(pc_clmFTR<0,-1,0)).transpose(2,0,1,3,4) # (nGCM,nPRD,nSIM,nGHM,nLand)
        agreementGCM =  where(flug_eachMME==flug_allMME,1,0).sum(axis=0) / float(nGCM)         #      (nPRD,nSIM,nGHM,nLand) >> Fig
        del flug_allMME, flug_eachMME

        """
        ### Convert nLand >> nY, nX                                       (Note: missing_value is -999)
        allMMEhis = fillinMasked(allMMEhis, lnd_mask)                          #      (nSIM,     nY,nX)
        allMMEftr = fillinMasked(allMMEftr, lnd_mask)                          # (nPRD,nSIM,     nY,nX)
        gcmMMEhis = fillinMasked(gcmMMEhis, lnd_mask)                          #      (nSIM,nGHM,nY,nX)
        gcmMMEftr = fillinMasked(gcmMMEftr, lnd_mask)                          # (nPRD,nSIM,nGHM,nY,nX)
        chng_allMME = fillinMasked(chng_allMME, lnd_mask)                      # (nPRD,nSIM,     nY,nX)
        chng_gcmMME = fillinMasked(chng_gcmMME, lnd_mask)                      # (nPRD,nSIM,nGHM,nY,nX)
        pc_allMME = fillinMasked(pc_allMME, lnd_mask)                          # (nPRD,nSIM,     nY,nX)
        pc_gcmMME = fillinMasked(pc_gcmMME, lnd_mask)                          # (nPRD,nSIM,nGHM,nY,nX)
        s2n_allMMEhis = fillinMasked(s2n_allMMEhis, lnd_mask)                  #      (nSIM,     nY,nX)
        s2n_pc_allMME = fillinMasked(s2n_pc_allMME, lnd_mask)                  # (nPRD,nSIM,     nY,nX)
        ratGCMvar_clmHIS = fillinMasked(ratGCMvar_clmHIS, lnd_mask)            #      (nSIM,     nY,nX)
        ratGCMvar_pc_allMME = fillinMasked(ratGCMvar_pc_allMME, lnd_mask)      # (nPRC,nSIM,     nY,nX)
        agreementALL = fillinMasked(agreementALL, lnd_mask)                    # (nPRD,nSIM,     nY,nX)
        agreementGCM = fillinMasked(agreementGCM, lnd_mask)                    # (nPRD,nSIM,nGHM,nY,nX)
        """


        # --------------------------------------------------------------------------------------------------------------
        # Make figure
        # --------------------------------------------------------------------------------------------------------------
        print("\nOK, let's start figure making...")
        #for prd,   srcType,                (mmeSrc,   agreementall),          s2nSrc,        (ghmSRC,   agreementgcm),uncertainty_source_rate in [
        #   ['hist',     'absolute',     (allMMEhis,           None),   s2n_allMMEhis,     (gcmMMEhis,           None),      ratGCMvar_clmHIS],
        #   ['mid',      'absolute',  (allMMEftr[0],           None),            None,  (gcmMMEftr[0],           None),                  None],
        #   ['late',     'absolute',  (allMMEftr[1],           None),            None,  (gcmMMEftr[1],           None),                  None],
        #   ['mid',     'abschange',(chng_allMME[0],agreementALL[0]),s2n_pc_allMME[0],(chng_gcmMME[0],agreementGCM[0]),                  None],
        #   ['late',    'abschange',(chng_allMME[1],agreementALL[1]),s2n_pc_allMME[1],(chng_gcmMME[1],agreementGCM[1]),                  None],
        #   ['mid', 'percentchange',  (pc_allMME[0],agreementALL[0]),s2n_pc_allMME[0],  (pc_gcmMME[0],agreementGCM[0]),ratGCMvar_pc_allMME[0]],
        #   ['late','percentchange',  (pc_allMME[1],agreementALL[1]),s2n_pc_allMME[1],  (pc_gcmMME[1],agreementGCM[1]),ratGCMvar_pc_allMME[1]],
        #   ]:
        for srcType,          prd,      (mmeSrc,      agreementall), s2nSrc,        (ghmSRC,      agreementgcm), uncertainty_source_rate in [
            ['absolute',      'hist',   (allMMEhis,   None),         None,          (gcmMMEhis,   None),         ratGCMvar_clmHIS],
            ['absolute',      'future', (allMMEftr,   None),         None,          (gcmMMEftr,   None),         None],
            ['abschange',     'future', (chng_allMME, agreementALL), s2n_pc_allMME, (chng_gcmMME, agreementGCM), None],
            ['percentchange', 'future', (pc_allMME,   agreementALL), s2n_pc_allMME, (pc_gcmMME,   agreementGCM), ratGCMvar_pc_allMME],
            ]:
            print('\n====================\n{} {}\n===================='.format(prd, srcType))

            mmeSrc = ma.masked_equal(ma.masked_equal(fillinMasked(mmeSrc, lnd_mask), 1e+20), 0) # his:(nSIM,     nY,nX),future:(nPRD,nSIM,     nY,nX)
            ghmSRC = ma.masked_equal(ma.masked_equal(fillinMasked(ghmSRC, lnd_mask), 1e+20), 0) # his:(nSIM,nGHM,nY,nX),future:(nPRD,nSIM,nGHM,nY,nX)
            if S2NMap and s2nSrc is not None:
                s2nSrc = ma.masked_equal(fillinMasked(s2nSrc, lnd_mask), 1e+20)                 # his:(nSIM,     nY,nX),future:(nPRD,nSIM,     nY,nX)
            if agreementall is not None: agreementall = ma.masked_equal(fillinMasked(agreementall, lnd_mask), 1e+20)  # (nPRD,nSIM,     nY,nX)
            if agreementgcm is not None: agreementgcm = ma.masked_equal(fillinMasked(agreementgcm, lnd_mask), 1e+20)  # (nPRD,nSIM,nGHM,nY,nX)
            if uncertainty_source_rate is not None:
                uncertainty_source_rate = ma.masked_equal(fillinMasked(uncertainty_source_rate, lnd_mask), 1e+20)  # his:(nSIM,nY,nX), future:(nPRC,nSIM,nY,nX)

            mmeSrc[np.isnan(mmeSrc)] = missing_value
            mmeSrc[np.isinf(mmeSrc)] = missing_value
            mmeSrc = np.ma.masked_equal(mmeSrc, missing_value)
            ghmSRC[np.isnan(ghmSRC)] = missing_value
            ghmSRC[np.isinf(ghmSRC)] = missing_value
            ghmSRC = np.ma.masked_equal(ghmSRC, missing_value)
            print('mmeSrc: {}-{} {}'.format(mmeSrc.min(), mmeSrc.max(), mmeSrc.shape))
            print('ghmSRC: {}-{} {}'.format(ghmSRC.min(), ghmSRC.max(), ghmSRC.shape))

            if EnsmMap:
                #=============#
                # 1. Ensemble #
                #=============#
                print('<<<<<<<<<<>>>>>>>>>>\n    Ensemble map    \n<<<<<<<<<<>>>>>>>>>>')

                if 'Basemap' in projection:
                    if projection == 'Basemap_PlateCarree':
                        figsize = (8, 3.8)
                        left, right, bottom, top, hspace, wspace =0.007, 0.995, 0.01, 0.975, 0.01, 0.01
                    elif projection == 'Basemap_Robinson':
                        figsize = (8, 4)
                        left, right, bottom, top, hspace, wspace =0.01, 0.995, 0.02, 0.95, 0.005, 0.002
                    dx, dx2, dx3, width = 0.028, 0.015, 0.001, 0.1
                    fig = plt.figure(num=1, figsize=figsize)
                    gs = gridspec.GridSpec(2, 2)  # (rows,cols)
                    gs.update(left=left, right=right, bottom=bottom, top=top, hspace=hspace, wspace=wspace)
                    ax1 = plt.subplot(gs[0,0])  # for mid-26
                    ax2 = plt.subplot(gs[0,1])  # for mid-85
                    ax3 = plt.subplot(gs[1,0])  # for late-26
                    ax4 = plt.subplot(gs[1,1])  # for late-8
                elif 'Cartopy' in projection:
                    dx, dx2, dx3, width = 0, -0.012, -0.025, 0.12
                    fig = plt.figure(num=1, figsize=(8.8, 4.3))
                    gs = gridspec.GridSpec(2, 2)  # (rows,cols)
                    gs.update(left=0.04, right=0.995, bottom=0.02, top=0.95, hspace=0.015, wspace=0.0005)
                    ax1 = plt.subplot(gs[0,0], projection=axis_projection)  # for mid-26
                    ax2 = plt.subplot(gs[0,1], projection=axis_projection)  # for mid-85
                    ax3 = plt.subplot(gs[1,0], projection=axis_projection)  # for late-26
                    ax4 = plt.subplot(gs[1,1], projection=axis_projection)  # for late-8
                axes = [ax1, ax2, ax3, ax4]

                if 'change' in srcType: changes = ['decrease', 'increase', 'significant_increase']
                else:                   changes = ['decrease', 'increase']
                df = pd.DataFrame(index=changes, columns=['{}-{}'.format(scn, future) for (scn, soc, co2) in SIMs for future in futures])

                #istep = 0
                for istep, ((ifuture, future), (iscn, (scn, soc, co2))) in enumerate(itertools.product(enumerate(futures), enumerate(SIMs))):  # [mid, late], [rcp26, rcp85]
                    #if not (ifuture == 0 and iscn == 0): istep = istep + 1
                    print(ifuture, iscn, istep)
                    tag = ['a', 'b', 'c', 'd'][istep]
                    print('ax{}..'.format(istep+1))
                    ax = axes[istep]
                    plt.sca(ax)
                    if 'Basemap' in projection:
                        ax.axis('off')
                        if KSmask_gray:
                            bm1.fillcontinents(color='#d2d2d2', zorder=-1)
                    elif 'Cartopy' in projection:
                        ax.set_extent([-180, 180, -60, 90], original_projection)
                        ax.outline_patch.set_linewidth(0.1)
                        ax.outline_patch.set_edgecolor('#b0b0b0')
                        ax.add_feature(cfeature.LAND, facecolor='#fdfdfd')
                    ax_pos = ax.get_position()
                    ax.text(0.5, 1.0, scn.upper(), ha="center", va="bottom", fontsize=7, transform=ax.transAxes)

                    if srcType == 'percentchange' or srcType == 'abschange':

                        aSrc = ma.masked_array(ma.masked_equal(mmeSrc[ifuture, iscn],0), mask=dictKSmask[scn][future]['all'])
                        print('{} {}_{}_{}: {}-{}'.format(future, scn, soc, co2, aSrc.min(), aSrc.max()))
                        if srcType == 'percentchange': 
                            bounds = [-200, -100, -50, -10, 0, 10, 50, 100, 200]
                            #bounds = [-100, -50, -10, 0, 10, 50, 100, 200]
                        elif srcType == 'abschange': bounds = [-20, -10, -5, 0, 5, 10, 20, 30]
                        print('bounds: {}'.format(bounds))

                        if len(agreeThrsh) == 3:
                            mask1  = ma.make_mask(agreementall[ifuture, iscn]<agreeThrsh[0])
                            mask21 = ma.make_mask(agreementall[ifuture, iscn]>=agreeThrsh[0])
                            mask22 = ma.make_mask(agreementall[ifuture, iscn]<agreeThrsh[1])
                            mask2  = ma.mask_or(mask21, mask22)
                            mask31 = ma.make_mask(agreementall[ifuture, iscn]>=agreeThrsh[1])
                            mask32 = ma.make_mask(agreementall[ifuture, iscn]<agreeThrsh[2])
                            mask3  = ma.mask_or(mask31, mask32)
                            mask4  = ma.make_mask(agreementall[ifuture, iscn]>=agreeThrsh[2])
                            signal1 = ma.masked_array(aSrc, mask=mask1); print('signal1: {}-{}'.format(signal1.min(), signal1.max()))
                            signal2 = ma.masked_array(aSrc, mask=mask2); print('signal2: {}-{}'.format(signal2.min(), signal2.max()))
                            signal3 = ma.masked_array(aSrc, mask=mask3); print('signal3: {}-{}'.format(signal3.min(), signal3.max()))
                            signal4 = ma.masked_array(aSrc, mask=mask4); print('signal4: {}-{}'.format(signal4.min(), signal4.max()))
                            colors1 = divide([[ 23, 23,171],[ 39,161,242],[127,  0,255],[  0,132,132],[240,210,  0],[230,120,  0],[170,  0,  0]],255.)
                            colors2 = divide([[102,109,206],[125,202,253],[175, 94,255],[ 84,172,172],[244,224, 84],[238,164, 84],[198, 84, 84]],255.)
                            colors3 = divide([[195,195,240],[164,219,255],[215,175,255],[168,213,213],[249,239,168],[246,209,168],[226,168,168]],255.)
                            colors4 = divide([[210,210,210],[210,210,210],[210,210,210],[210,210,210],[210,210,210],[210,210,210],[210,210,210]],255.)  # gray d2d2d2
                            cmap1 = mpl.colors.ListedColormap(colors1)
                            cmap2 = mpl.colors.ListedColormap(colors2)
                            cmap3 = mpl.colors.ListedColormap(colors3)
                            cmap4 = mpl.colors.ListedColormap(colors4)
                            norm1 = mpl.colors.BoundaryNorm(bounds, cmap1.N)
                            norm2 = mpl.colors.BoundaryNorm(bounds, cmap2.N)
                            norm3 = mpl.colors.BoundaryNorm(bounds, cmap3.N)
                            norm4 = mpl.colors.BoundaryNorm(bounds, cmap4.N)
                            if projection == 'Basemap_PlaceCarree':
                                ims3 = bm1.imshow(flipud(signal3[11:293]), norm=norm3, cmap=cmap3, interpolation='nearest')
                                ims2 = bm1.imshow(flipud(signal2[11:293]), norm=norm2, cmap=cmap2, interpolation='nearest')
                                ims1 = bm1.imshow(flipud(signal1[11:293]), norm=norm1, cmap=cmap1, interpolation='nearest')
                                ims4 = bm1.imshow(flipud(signal4[11:293]), norm=norm4, cmap=cmap4, interpolation='nearest')
                                bm1.drawcoastlines(linewidth=coastal_line_width, color=coastal_line_color)
                            elif projection == 'Basemap_Robin':
                                ims3 = bm1.imshow(bm1.fransform_scalar(flipud(signal3[11:293])), norm=norm3, cmap=cmap3, interpolation='nearest')
                                ims2 = bm1.imshow(bm1.fransform_scalar(flipud(signal2[11:293])), norm=norm2, cmap=cmap2, interpolation='nearest')
                                ims1 = bm1.imshow(bm1.fransform_scalar(flipud(signal1[11:293])), norm=norm1, cmap=cmap1, interpolation='nearest')
                                ims4 = bm1.imshow(bm1.fransform_scalar(flipud(signal4[11:293])), norm=norm4, cmap=cmap4, interpolation='nearest')
                                bm1.drawcoastlines(linewidth=coastal_line_width, color=coastal_line_color)
                            elif 'Cartopy' in projection:
                                ##img_extent = [-180, 180, -60, 90]
                                ##im = ax.imshow(src, origin='upper', vmin=int(vmin), vmax=int(vmax), norm=norm, cmap=cmap, transform=ccrs.PlateCarree(), extent=img_extent)
                                #im3 = plt.contourf(lons, lats, signal3[:293,:], levels=bounds, norm=norm3, cmap=cmap3, transform=original_projection)
                                #im2 = plt.contourf(lons, lats, signal2[:293,:], levels=bounds, norm=norm2, cmap=cmap2, transform=original_projection)
                                #im1 = plt.contourf(lons, lats, signal1[:293,:], levels=bounds, norm=norm1, cmap=cmap1, transform=original_projection)
                                #im4 = plt.contourf(lons, lats, signal4[:293,:], levels=bounds, norm=norm4, cmap=cmap4, transform=original_projection)
                                im3 = plt.contourf(lons, lats, signal3[:293,:], levels=bounds, cmap=cmap3, transform=original_projection)
                                im2 = plt.contourf(lons, lats, signal2[:293,:], levels=bounds, cmap=cmap2, transform=original_projection)
                                im1 = plt.contourf(lons, lats, signal1[:293,:], levels=bounds, cmap=cmap1, transform=original_projection)
                                im4 = plt.contourf(lons, lats, signal4[:293,:], levels=bounds, cmap=cmap4, transform=original_projection)
                                ax.add_feature(cfeature.COASTLINE, linewidth=coastal_line_width, edgecolor=coastal_line_color)
                            if tag == 'c':  # 2D colorbar for Kaye et al.-plot:
                                ax14 = fig.add_axes([ax_pos.x0+dx, ax_pos.y0+0.140, 0.10, 0.02])
                                ax13 = fig.add_axes([ax_pos.x0+dx, ax_pos.y0+0.120, 0.10, 0.02])
                                ax12 = fig.add_axes([ax_pos.x0+dx, ax_pos.y0+0.100, 0.10, 0.02])
                                ax11 = fig.add_axes([ax_pos.x0+dx, ax_pos.y0+0.080, 0.10, 0.02])
                                cmap = [cmap1, cmap2, cmap3, cmap4]
                                for iax, axs in enumerate([ax12, ax13, ax14]):
                                    norm = mpl.colors.BoundaryNorm(bounds, cmap[iax+1].N)
                                    cb = mpl.colorbar.ColorbarBase(axs, cmap=cmap[iax+1], norm=norm, boundaries=bounds, spacing='uniform', orientation='horizontal')
                                    cb.set_ticks(bounds)
                                    cb.set_ticklabels([])
                                    cb.ax.tick_params(labelsize=5, width=0.25, direction='in')
                                    cb.outline.set_visible(False)
                                norm = mpl.colors.BoundaryNorm(bounds, cmap[0].N)
                                cb1 = mpl.colorbar.ColorbarBase(ax11, cmap=cmap[0], norm=norm, boundaries=bounds, spacing='uniform', orientation='horizontal')
                                cb1.ax.set_xticks(bounds)
                                cb1.ax.set_xticklabels([str(int(ibound)) for ibound in bounds], rotation=-45, fontsize=4.5)
                                cb1.ax.tick_params(labelsize=5, width=0.25, direction='in')
                                # labels and ticks
                                if srcType == 'percentchange': cb1.set_label('relative change [%]',                      fontsize=4.5, labelpad=-0.9)
                                elif srcType == 'abschange':   cb1.set_label('{} [{}]'.format(srcType, dictUnit[index]), fontsize=4.5, labelpad=-0.6)
                                cb1.outline.set_visible(False)
                                fig.text(ax_pos.x0+0.015, ax_pos.y0+0.140, str(int(agreeThrsh[2]*1e2)), va='center', ha='center', fontsize=4.5)
                                fig.text(ax_pos.x0+0.015, ax_pos.y0+0.120, str(int(agreeThrsh[1]*1e2)), va='center', ha='center', fontsize=4.5)
                                fig.text(ax_pos.x0+0.015, ax_pos.y0+0.100, str(int(agreeThrsh[0]*1e2)), va='center', ha='center', fontsize=4.5)
                                fig.text(ax_pos.x0+0.001, ax_pos.y0+0.130, 'agreement [%]',             va='center', ha='center', fontsize=4.5, rotation='vertical')
                        elif len(agreeThrsh) == 2:
                            mask1  = ma.make_mask(agreementall[ifuture, iscn]<agreeThrsh[0])   # mask <80%
                            mask21 = ma.make_mask(agreementall[ifuture, iscn]>=agreeThrsh[0])
                            mask22 = ma.make_mask(agreementall[ifuture, iscn]<agreeThrsh[1])
                            mask2  = ma.mask_or(mask21, mask22)                                 # mask <60% 80%<
                            mask3  = ma.make_mask(agreementall[ifuture, iscn]>=agreeThrsh[1])   # mask 60%<
                            signal1 = ma.masked_array(aSrc, mask=mask1); print('signal1: {}-{}'.format(signal1.min(), signal1.max()))
                            signal2 = ma.masked_array(aSrc, mask=mask2); print('signal2: {}-{}'.format(signal2.min(), signal2.max()))
                            signal3 = ma.masked_array(aSrc, mask=mask3); print('signal3: {}-{}'.format(signal3.min(), signal3.max()))
                            ### Map   (reference: Kaye et al., Fig. 7b)
                            colors1 = divide([[ 23, 23,171],[ 39,161,242],[127,  0,255],[  0,132,132],[240,210,  0],[230,120,  0],[170,  0,  0]],255.)
                            colors2 = divide([[102,109,206],[125,202,253],[175, 94,255],[ 84,172,172],[244,224, 84],[238,164, 84],[198, 84, 84]],255.)
                            colors3 = divide([[210,210,210],[210,210,210],[210,210,210],[210,210,210],[210,210,210],[210,210,210],[210,210,210]],255.)  # gray
                            cmap1 = mpl.colors.ListedColormap(colors1)
                            cmap2 = mpl.colors.ListedColormap(colors2)
                            cmap3 = mpl.colors.ListedColormap(colors3)  # gray
                            norm1 = mpl.colors.BoundaryNorm(bounds, cmap1.N)
                            norm2 = mpl.colors.BoundaryNorm(bounds, cmap2.N)
                            norm3 = mpl.colors.BoundaryNorm(bounds, cmap3.N)
                            if projection == 'Basemap_PlateCarree':
                                ims1 = bm1.imshow(flipud(signal1[11:293]), norm=norm1, cmap=cmap1, interpolation='nearest')
                                ims2 = bm1.imshow(flipud(signal2[11:293]), norm=norm2, cmap=cmap2, interpolation='nearest')
                                ims3 = bm1.imshow(flipud(signal3[11:293]), norm=norm3, cmap=cmap3, interpolation='nearest')
                                bm1.drawcoastlines(linewidth=coastal_line_width, color=coastal_line_color)
                            elif projection == 'Basemap_Robinson':
                                ims1 = bm1.imshow(bm1.transform_scalar(flipud(signal1), lons, lats, nx, ny), norm=norm1, cmap=cmap1, interpolation='nearest')
                                ims2 = bm1.imshow(bm1.transform_scalar(flipud(signal2), lons, lats, nx, ny), norm=norm2, cmap=cmap2, interpolation='nearest')
                                ims3 = bm1.imshow(bm1.transform_scalar(flipud(signal3), lons, lats, nx, ny), norm=norm3, cmap=cmap3, interpolation='nearest')
                                bm1.drawcoastlines(linewidth=coastal_line_width, color=coastal_line_color)
                            elif 'Cartopy' in projection:
                                ###img_extent = [-180, 180, -60, 90]
                                ###im = ax.imshow(src, origin='upper', vmin=int(vmin), vmax=int(vmax), norm=norm, cmap=cmap, transform=ccrs.PlateCarree(), extent=img_extent)
                                ###im3 = ax.contourf(lons, lats, signal3[:293], levels=bounds, norm=norm3, cmap=cmap3, transform=original_projection)
                                ###im2 = ax.contourf(lons, lats, signal2[:293], levels=bounds, norm=norm2, cmap=cmap2, transform=original_projection)
                                ###im1 = ax.contourf(lons, lats, signal1[:293], levels=bounds, norm=norm1, cmap=cmap1, transform=original_projection)
                                im3 = ax.contourf(lons, lats, signal3[:293], levels=bounds, cmap=cmap3, transform=original_projection)
                                im2 = ax.contourf(lons, lats, signal2[:293], levels=bounds, cmap=cmap2, transform=original_projection)
                                im1 = ax.contourf(lons, lats, signal1[:293], levels=bounds, cmap=cmap1, transform=original_projection)
                                ax.add_feature(cfeature.COASTLINE, linewidth=coastal_line_width, edgecolor=coastal_line_color)
                            if tag == 'a' or tag == 'c':  # 2D colorbar for Kaye et al.-plot:
                                ax13 = fig.add_axes([ax_pos.x0+dx, ax_pos.y0+0.160, width, 0.03])
                                ax12 = fig.add_axes([ax_pos.x0+dx, ax_pos.y0+0.130, width, 0.03])
                                ax11 = fig.add_axes([ax_pos.x0+dx, ax_pos.y0+0.100, width, 0.03])
                                #ax14 = fig.add_axes([ax_pos.x0+0.028, ax_pos.y0+0.190, 0.10, 0.03])  # dummy...
                                cmap = [cmap1, cmap2, cmap3]
                                for iax, axs in enumerate([ax12, ax13]):
                                    norm = mpl.colors.BoundaryNorm(bounds, cmap[iax+1].N)
                                    cb = mpl.colorbar.ColorbarBase(axs, cmap=cmap[iax+1], norm=norm, boundaries=bounds, spacing='uniform', orientation='horizontal')
                                    cb.set_ticks(bounds)
                                    cb.set_ticklabels([])
                                    cb.ax.tick_params(labelsize=5, width=0.25, direction='in')
                                    cb.outline.set_visible(False)
                                norm = mpl.colors.BoundaryNorm(bounds, cmap[0].N)
                                cb1 = mpl.colorbar.ColorbarBase(ax11, cmap=cmap[0], norm=norm, boundaries=bounds, spacing='uniform', orientation='horizontal')
                                cb1.set_ticks(bounds)
                                #cb1.ax.set_xticks(bounds)
                                cb1.ax.set_xticklabels([str(int(ibound)) for ibound in bounds], rotation=-45, fontsize=6)
                                cb1.ax.tick_params(labelsize=5, width=0.25, direction='in')
                                if srcType == 'percentchange': cb1.set_label('Relative change [%]',                      fontsize=6.5, labelpad=-0.9)
                                elif srcType == 'abschange':   cb1.set_label('{} [{}]'.format(srcType, dictUnit[index]), fontsize=6.5, labelpad=-0.6)
                                cb1.outline.set_visible(False)
                                fig.text(ax_pos.x0+dx2, ax_pos.y0+0.155, str(int(agreeThrsh[1]*1e2)), ha='center', va='center', fontsize=6)
                                fig.text(ax_pos.x0+dx2, ax_pos.y0+0.125, str(int(agreeThrsh[0]*1e2)), ha='center', va='center', fontsize=6)
                                fig.text(ax_pos.x0+dx2, ax_pos.y0+0.095, str(int(1e2)),               ha='center', va='center', fontsize=6)
                                fig.text(ax_pos.x0+dx3, ax_pos.y0+0.135, 'Agreement [%]',             ha='center', va='center', fontsize=6.5, rotation='vertical')
                        elif len(agreeThrsh) == 1:
                            dx, width = 0.005, 0.123
                            mask1  = ma.make_mask(agreementall[ifuture, iscn]<agreeThrsh[0])   # mask <60%
                            mask2  = ma.make_mask(agreementall[ifuture, iscn]>=agreeThrsh[0])   # mask 60%<
                            signal1 = ma.masked_array(aSrc, mask=mask1); print('signal1: {}-{}'.format(signal1.min(), signal1.max()))
                            signal2 = ma.masked_array(aSrc, mask=mask2); print('signal2: {}-{}'.format(signal2.min(), signal2.max()))
                            ### Map   (reference: Kaye et al., Fig. 7b)
                            #colors1 = divide([[ 68,116,180],[143,190,218],[221,241,246],[254,254,190],[252,225,143],[250,140, 87],[213, 45, 37]],255.)
                            colors1 = divide([[ 68,116,180],[125,165,225],[179,207,229],[221,241,246],[254,254,190],[252,225,143],[250,140, 87],[213, 45, 37]],255.)
                            colors2 = divide([[210,210,210],[210,210,210],[210,210,210],[210,210,210],[210,210,210],[210,210,210],[210,210,210],[210,210,210]],255.)  # gray
                            cmap1 = mpl.colors.ListedColormap(colors1)
                            cmap2 = mpl.colors.ListedColormap(colors2)  # gray
                            norm1 = mpl.colors.BoundaryNorm(bounds, cmap1.N)
                            norm2 = mpl.colors.BoundaryNorm(bounds, cmap2.N)
                            if projection == 'Basemap_PlateCarree':
                                ims1 = bm1.imshow(flipud(signal1[11:293]), norm=norm1, cmap=cmap1, interpolation='nearest')
                                ims2 = bm1.imshow(flipud(signal2[11:293]), norm=norm2, cmap=cmap2, interpolation='nearest')
                                bm1.drawcoastlines(linewidth=coastal_line_width, color=coastal_line_color)
                            elif projection == 'Basemap_Robinson':
                                ims1 = bm1.imshow(bm1.transform_scalar(flipud(signal1), lons, lats, nx, ny), norm=norm1, cmap=cmap1, interpolation='nearest')
                                ims2 = bm1.imshow(bm1.transform_scalar(flipud(signal2), lons, lats, nx, ny), norm=norm2, cmap=cmap2, interpolation='nearest')
                                bm1.drawcoastlines(linewidth=coastal_line_width, color=coastal_line_color)
                            elif 'Cartopy' in projection:
                                ###img_extent = [-180, 180, -60, 90]
                                ###im = ax.imshow(src, origin='upper', vmin=int(vmin), vmax=int(vmax), norm=norm, cmap=cmap, transform=ccrs.PlateCarree(), extent=img_extent)
                                ###im3 = ax.contourf(lons, lats, signal3[:293], levels=bounds, norm=norm3, cmap=cmap3, transform=original_projection)
                                ###im2 = ax.contourf(lons, lats, signal2[:293], levels=bounds, norm=norm2, cmap=cmap2, transform=original_projection)
                                ###im1 = ax.contourf(lons, lats, signal1[:293], levels=bounds, norm=norm1, cmap=cmap1, transform=original_projection)
                                im2 = ax.contourf(lons, lats, signal2[:293], levels=bounds, cmap=cmap2, transform=original_projection)
                                im1 = ax.contourf(lons, lats, signal1[:293], levels=bounds, cmap=cmap1, transform=original_projection)
                                ax.add_feature(cfeature.COASTLINE, linewidth=coastal_line_width, edgecolor=coastal_line_color)
                            if tag == 'a' or tag == 'c':  # 2D colorbar for Kaye et al.-plot:
                                cmap = [cmap1, cmap2]  #####, cmap3]
                                norm = mpl.colors.BoundaryNorm(bounds, cmap[0].N)
                                ax11 = fig.add_axes([ax_pos.x0+dx, ax_pos.y0+0.130, width, 0.03])  # for horizontal
                                cb1 = mpl.colorbar.ColorbarBase(ax11, cmap=cmap[0], norm=norm, boundaries=bounds, spacing='uniform', orientation='horizontal')
                                #ax11 = fig.add_axes([ax_pos.x0+3*dx, ax_pos.y0+0.03, 0.03, width])  # for vertical
                                #cb1 = mpl.colorbar.ColorbarBase(ax11, cmap=cmap[0], norm=norm, boundaries=bounds, spacing='uniform', orientation='vertical')
                                cb1.set_ticks(bounds)
                                #cb1.ax.set_xticklabels([str(int(ibound)) for ibound in bounds], fontsize=6.5, rotation=-45)
                                cb1.ax.set_xticklabels([str(int(ibound)) for ibound in bounds], fontsize=6.5, rotation=-50)
                                cb1.ax.tick_params(labelsize=6.5, width=0.25, direction='in')
                                #if srcType == 'percentchange': cb1.set_label('Relative change [%]',                      fontsize=7.5, labelpad=-1.1)
                                if srcType == 'percentchange': cb1.set_label('Relative change [%]',                      fontsize=7.5, labelpad=2.5)
                                elif srcType == 'abschange':   cb1.set_label('{} [{}]'.format(srcType, dictUnit[index]), fontsize=7.5, labelpad=-0.6)
                                cb1.outline.set_visible(False)
                                #####fig.text(ax_pos.x0+dx2, ax_pos.y0+0.155, str(int(agreeThrsh[1]*1e2)), ha='center', va='center', fontsize=6)
                                #####fig.text(ax_pos.x0+dx2, ax_pos.y0+0.125, str(int(agreeThrsh[0]*1e2)), ha='center', va='center', fontsize=6)
                                #####fig.text(ax_pos.x0+dx2, ax_pos.y0+0.095, str(int(1e2)),               ha='center', va='center', fontsize=6)
                                #####fig.text(ax_pos.x0+dx3, ax_pos.y0+0.135, 'Agreement [%]',             ha='center', va='center', fontsize=6.5, rotation='vertical')

                        """
                        if tag == 'd': 
                            ax.text(0.38, 0.01, 
                                    '{}\n (Simulation: {}, {}_{}  (Period: {}-historical  (Drought: Q{}, Len{}'.format(srcType, scn, scn, soc, prd, Q, Len),
                                    va='bottom', ha='left', fontsize=5, transform=ax.transAxes)
                        """

                    else:   # absolute

                        if prd == 'hist' and srcType == 'absolute': 
                            aSrc = mmeSrc[iscn]
                        elif prd == 'future' and srcType == 'absolute':
                            aSrc = mmeSrc[ifuture, iscn]
                        bounds = [0, 10, 20, 30, 40, 50]
                        cmap = cm.rainbow
                        norm1 = mpl.colors.BoundaryNorm(bounds, cmap.N)

                        if projection == 'Basemap_PlateCarree':
                            ims1 = bm1.imshow(flipud(aSrc[11:293]), norm=norm1, cmap=cmap, interpolation="nearest")
                            bm1.drawcoastlines(linewidth=coastal_line_width, color=coastal_line_color)
                        elif projection == 'Basemap_Robinson':
                            ims1 = bm1.imshow(bm1.tfansform_scalar(flipud(aSrc[11:293])), norm=norm1, cmap=cmap, interpolation="nearest")
                            bm1.drawcoastlines(linewidth=coastal_line_width, color=coastal_line_color)
                        elif 'Cartopy' in projection:
                            ##img_extent = [-180, 180, -60, 90]
                            ##im = ax.imshow(src, origin='upper', vmin=int(vmin), vmax=int(vmax), norm=norm, cmap=cmap, transform=original_projection, extent=img_extent)
                            #ims1 = plt.contourf(lons, lats, aSrc[:293,:], levels=bounds, norm=norm1, cmap=cmap, transform=original_projection)
                            ims1 = plt.contourf(lons, lats, aSrc[:293,:], levels=bounds, cmap=cmap, transform=original_projection)
                            ax.add_feature(cfeature.COASTLINE, linewidth=coastal_line_width, edgecolor=coastal_line_color)

                        if tag == 'c':
                            ax11 = fig.add_axes([ax_pos.x0+0.02, ax_pos.y0+0.080,  0.10, 0.03])
                            cb1 = mpl.colorbar.ColorbarBase(ax11, cmap=cmap, norm=norm1, orientation='horizontal')
                            cb1.set_ticks(bounds)
                            cb1.set_ticklabels([str(int(ibound)) for ibound in bounds])
                            cb1.ax.tick_params(labelsize=5, width=0.25, direction='in')
                            cb1.outline.set_visible(False)

                            if srcType == 'absolute':
                                cb1.set_label('[{}]'.format(dictUnit[index]), fontsize=5, labelpad=-0.6)

                        """
                        if tag == 'd':
                            if srcType == 'abschange':
                                ax.text(0.38, 0.01, 
                                         '{}\n{}, {}_{}  {}-historical Q{}, Len{}'.format(srcType, scn, scn, soc, prd, Q, Len),
                                         va='bottom', ha='left', fontsize=5, transform=ax.transAxes)
                            elif srcType == 'absolute':
                                ax.text(0.38, 0.01,
                                        '{}\n{}, {}_{}  {}   Q{}, Len{}'.format(srcType, scn, scn, soc, prd, Q, Len),
                                         va='bottom', ha='left', fontsize=5, transform=ax.transAxes)
                        """

                    #if tag == 'a' or tag == 'b':
                    #if tag == 'b':
                    if region_type == 'AR6_regions' and tag == 'b' and season == 'DRY':
                        if 'Basemap' in projection:
                            bm1.readshapefile(ipcc_shp_path, shp_name, linewidth=0.6, color='k')
                        elif 'Cartopy' in projection:
                            shape_feature = ShapelyFeature(Reader(ipcc_shp_path).geometries(), axis_projection, edgecolor='#000000', linewidth=0.6, facecolor='none')
                            ax.add_feature(shape_feature)

                    # information of change
                    if 'change' in srcType:
                        if   len(agreeThrsh) == 3: low_agreement_mask = ma.make_mask(agreementall[ifuture, iscn]<agreeThrsh[2])
                        elif len(agreeThrsh) == 2: low_agreement_mask = ma.make_mask(agreementall[ifuture, iscn]<agreeThrsh[1])
                        elif len(agreeThrsh) == 1: low_agreement_mask = ma.make_mask(agreementall[ifuture, iscn]<agreeThrsh[0])
                        _aSrc = ma.masked_array(aSrc, mask=low_agreement_mask)
                        decrease_rate = change_arearate('decrease', _aSrc)
                        increase_rate = change_arearate('increase', _aSrc)
                        df.loc['decrease', '{}-{}'.format(scn, future)] = decrease_rate
                        df.loc['increase', '{}-{}'.format(scn, future)] = increase_rate
                        if srcType == 'percentchange':
                            significant_level = 100
                            significant_increase_rate = change_arearate('significant_increase', _aSrc, significant_level)
                            df.loc['significant_increase', '{}-{}'.format(scn, future)] = significant_increase_rate

                for suffix in suffixes:
                    outFile = 'globMap.Ensemble{}.{}.{}.{}.Q{:02}win{:02}_Len{:03}tau{}_{}.{}'.format(mmeType, srcType, prd, season, 
                                                                                                      Q, win, Len, tau, index, suffix)
                    outDir = os.path.join(fig_directory_main, season, suffix)
                    #outDir = os.path.join(fig_directory_main, '{}.{}_{}'.format(qvalType,syear_ref,eyear_ref), 
                    #                      '{}-{}.{}-{}.{}-{}'.format(syear_hi, eyear_hi, syear_mf, eyear_mf, syear_ff, eyear_ff))
                    if not os.path.isdir(outDir): os.makedirs(outDir)
                    outPath = os.path.join(outDir, outFile)
                    if os.path.exists(outPath):
                        print("File exists, will be overwritten. ")
                        #raw_input("Print figure? Press key to continue...")
                    plt.savefig(outPath, dpi=dpi)
                    print('savefig: {}\n'.format(outPath))
                plt.close(1)

                change_excel_name = 'globMap.Ensemble{}.{}.{}.Q{:02}win{:02}_Len{:03}tau{}_{}.xlsx'.format(mmeType, srcType, 
                                                                                                           season, Q, win, Len, tau, index)
                change_excel_path = os.path.join(outDir, change_excel_name)
                writer = pd.ExcelWriter(change_excel_path)
                df.to_excel(writer, 'change_area_rate')
                writer.close()
                print('ExcelWriter: {}\n'.format(change_excel_path))


            if AGREEMENT and agreementall is not None:
                #==============#
                # 2. Agreement #
                #==============#

                print('<<<<<<<<<<>>>>>>>>>>\n    Agreement map    \n<<<<<<<<<<>>>>>>>>>>')
                if 'Basemap' in projection:
                    if projection == 'Basemap_PlateCarree':
                        dx, width = 0.028, 0.03
                        figsize = (8, 3.8) 
                        left, right, bottom, top, hspace, wspace = 0.007, 0.995, 0.01, 0.975, 0.01, 0.01
                    elif projection == 'Basemap_Robinson':
                        dx, dx2, dx3, width = 0.028, 0.015, 0.001, 0.1
                        figsize = (8, 4)
                        left, right, bottom, top, hspace, wspace =0.01, 0.995, 0.02, 0.95, 0.005, 0.002
                    fig = plt.figure(num=2, figsize=figsize)
                    gs = gridspec.GridSpec(2, 2)  # (rows,cols)
                    gs.update(left=left, right=right, bottom=bottom, top=top, hspace=hspace, wspace=wspace)
                    ax1 = plt.subplot(gs[0,0])  # for mid-26
                    ax2 = plt.subplot(gs[0,1])  # for mid-85
                    ax3 = plt.subplot(gs[1,0])  # for late-26
                    ax4 = plt.subplot(gs[1,1])  # for late-8
                elif 'Cartopy' in projection:
                    dx, dx2, dx3, width = 0, -0.012, -0.025, 0.12
                    fig = plt.figure(num=2, figsize=(8.8, 4.3))
                    gs = gridspec.GridSpec(2, 2)  # (rows,cols)
                    gs.update(left=0.04, right=0.995, bottom=0.02, top=0.95, hspace=0.015, wspace=0.0005)
                    ax1 = plt.subplot(gs[0,0], projection=axis_projection)  # for mid-26
                    ax2 = plt.subplot(gs[0,1], projection=axis_projection)  # for mid-85
                    ax3 = plt.subplot(gs[1,0], projection=axis_projection)  # for late-26
                    ax4 = plt.subplot(gs[1,1], projection=axis_projection)  # for late-8
                axes = [ax1, ax2, ax3, ax4]

                istep = 0
                for istep, ((ifuture, future), (iscn, (scn, soc, co2))) in enumerate(itertools.product(enumerate(futures), enumerate(SIMs))):  # [mid,late],[rcp26,rcp85]
                    ax = axes[istep]
                    tag = ['a', 'b', 'c', 'd'][istep]

                    print('check... {}-{}'.format(agreementall[ifuture,iscn,:293,:].min(), agreementall[ifuture,iscn,:293,:].max()))

                    print('ax{}..'.format(istep+1))
                    plt.sca(ax)
                    if 'Basemap' in projection:
                        ax.axis('off')
                    elif 'Cartopy' in projection:
                        ax.set_extent([-180, 180, -60, 90], original_projection)
                        ax.outline_patch.set_linewidth(0.1)
                        ax.outline_patch.set_edgecolor('#b0b0b0')
                        ax.add_feature(cfeature.LAND, facecolor='#fdfdfd')

                    bounds = [0, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
                    colors1 = divide([[235,235,235], [210,210,210], [255,240,168], [253,160,68], [216,18,30], [129,0,38]], 255.)
                    labelName = 'agreement'
                    cmap = colors.ListedColormap(colors1)
                    norm = colors.BoundaryNorm(bounds, cmap.N)

                    if projection == 'Basemap_PlateCarree':
                        ims = bm1.imshow(flipud(agreementall[ifuture,iscn,11:293,:]), norm=norm, cmap=cmap, interpolation="nearest")
                        bm1.drawcoastlines(linewidth=coastal_line_width, color=coastal_line_color)
                    elif projection == 'Basemap_Robinson':
                        ims = bm1.imshow(bm1.transform_scalar(flipud(agreementall[ifuture,iscn]), lons, lats, nx, ny), norm=norm, cmap=cmap, interpolation="nearest")
                        bm1.drawcoastlines(linewidth=coastal_line_width, color=coastal_line_color)
                    elif 'Cartopy' in projection:
                        ##img_extent = [-180, 180, -60, 90]
                        ##im = ax.imshow(src, origin='upper', vmin=int(vmin), vmax=int(vmax), norm=norm, cmap=cmap, transform=original_projection, extent=img_extent)
                        #ims = ax.contourf(lons, lats, s2nSrc[ifuture, iscn,:293], levels=bounds, norm=norm, cmap=cmap, transform=original_projection)
                        ims = ax.contourf(lons, lats, agreementall[ifuture, iscn,:293], levels=bounds, cmap=cmap, transform=original_projection)
                        ax.add_feature(cfeature.COASTLINE, linewidth=coastal_line_width, edgecolor=coastal_line_color)

                    ax.text(0., 1.,  '({})'.format(tag),         va="bottom", ha="left",   fontsize=8, transform=ax.transAxes)
                    ax.text(0.5, 1., '{} {}'.format(scn,future), va="bottom", ha="center", fontsize=8, transform=ax.transAxes)

                    if tag == 'c':
                        ax_pos = ax.get_position()
                        ax11 = fig.add_axes([ax_pos.x0+0.008, ax_pos.y0+0.1,  0.112, 0.03])
                        cb1 = mpl.colorbar.ColorbarBase(ax11,
                                                        cmap=cmap,
                                                        norm=norm,
                                                        boundaries=bounds,
                                                        #spacing='proportional',
                                                        spacing='uniform',
                                                        orientation='horizontal')
                        #for t in cb1.ax11.get_yticklabels(): t.set_fontsize(6)
                        cb1.set_ticks(bounds)
                        cb1.set_ticklabels([str(i) for i in bounds])
                        cb1.ax.tick_params(labelsize=5, width=0.25, direction='in')
                        cb1.set_label(labelName, fontsize=7.5)#, labelpad=-0.6)
                        cb1.outline.set_visible(False)

                for suffix in suffixes:
                    outFile = 'globMap.agreement.Ensemble{}.{}.{}.{}.Q{:02}win{:02}_Len{:03}tau{}_{}.{}'.format(mmeType, srcType, prd, season,
                                                                                                                Q, win, Len, tau, index, suffix)
                    outDir = os.path.join(fig_directory_main, season, suffix)
                    #outDir = os.path.join(fig_directory_main, '{}.{}_{}'.format(qvalType,syear_ref,eyear_ref),
                    #                      '{}-{}.{}-{}.{}-{}'.format(syear_hi, eyear_hi, syear_mf, eyear_mf, syear_ff, eyear_ff),
                    #                      )
                    if not os.path.isdir(outDir): os.makedirs(outDir)
                    outPath = os.path.join(outDir, outFile)

                    if os.path.exists(outPath):
                        print("File exists, will be overwritten. ")
                        #raw_input("Print figure? Press key to continue...")
                    plt.savefig(outPath, dpi=dpi)
                    print('savefig: {}\n'.format(outPath))
                plt.close(2)


            if S2NMap and s2nSrc is not None:
                #========#
                # 3. S2N #
                #========#
                
                print('<<<<<<<<<<>>>>>>>>>>\n       S2N map       \n<<<<<<<<<<>>>>>>>>>>')
                if 'Basemap' in projection:
                    if projection == 'Basemap_PlateCarree':
                        dx, width = 0.028, 0.03
                        figsize = (8, 3.8)
                        left, right, bottom, top, hspace, wspace =0.007, 0.995, 0.01, 0.975, 0.01, 0.01
                    elif projection == 'Basemap_Robinson':
                        dx, dx2, dx3, width = 0.028, 0.015, 0.001, 0.1
                        figsize = (8, 4)
                        left, right, bottom, top, hspace, wspace =0.01, 0.995, 0.02, 0.95, 0.005, 0.002
                    fig = plt.figure(num=3, figsize=figsize)
                    gs = gridspec.GridSpec(2, 2)  # (rows,cols)
                    gs.update(left=left, right=right, bottom=bottom, top=top, hspace=hspace, wspace=wspace)
                    ax1 = plt.subplot(gs[0,0])  # for mid-26
                    ax2 = plt.subplot(gs[0,1])  # for mid-85
                    ax3 = plt.subplot(gs[1,0])  # for late-26
                    ax4 = plt.subplot(gs[1,1])  # for late-8
                elif 'Cartopy' in projection:
                    dx, dx2, dx3, width = 0, -0.012, -0.025, 0.12
                    fig = plt.figure(num=3, figsize=(8.8, 4.3))
                    gs = gridspec.GridSpec(2, 2)  # (rows,cols)
                    gs.update(left=0.04, right=0.995, bottom=0.02, top=0.95, hspace=0.015, wspace=0.0005)
                    ax1 = plt.subplot(gs[0,0], projection=axis_projection)  # for mid-26
                    ax2 = plt.subplot(gs[0,1], projection=axis_projection)  # for mid-85
                    ax3 = plt.subplot(gs[1,0], projection=axis_projection)  # for late-26
                    ax4 = plt.subplot(gs[1,1], projection=axis_projection)  # for late-8
                axes = [ax1, ax2, ax3, ax4]

                istep = 0
                for istep, ((ifuture, future), (iscn, (scn, soc, co2))) in enumerate(itertools.product(enumerate(futures), enumerate(SIMs))):  # [mid,late],[rcp26,rcp85]
                    #if not (ifuture == 0 and iscn == 0): istep = istep + 1
                    ax = axes[istep]
                    tag = ['a', 'b', 'c', 'd'][istep]

                    print('check... {}-{}'.format(s2nSrc[ifuture,iscn,:293,:].min(), s2nSrc[ifuture,iscn,:293,:].max()))

                    print('ax{}..'.format(istep+1))
                    plt.sca(ax)
                    if 'Basemap' in projection:
                        ax.axis('off')
                    elif 'Cartopy' in projection:
                        ax.set_extent([-180, 180, -60, 90], ccrs.PlateCarree())
                        ax.outline_patch.set_linewidth(0.1)
                        ax.outline_patch.set_edgecolor('#b0b0b0')
                        ax.add_feature(cfeature.LAND, facecolor='#fdfdfd')

                    if   mmeType == 'median':
                        bounds = [0,0.5,1,1.5,2,2.5,3,3.5]
                        #colors1 = divide([[28,125,199], [0,0,180], [238,172,172], [228,121,121], [239,0,0], [198,0,0], [158,0,0]], 255.)
                        colors1 = divide([[153,204,255], [0,0,180], [238,172,172], [228,121,121], [239,0,0], [198,0,0], [158,0,0]], 255.)
                        labelName = 'signal to noise ratio [-]'
                    elif mmeType == 'mean':
                        bounds = [0,0.05,0.1,0.5,1,1.5,2,2.5]
                        colors1 = divide([[28,125,199], [0,0,180], [238,172,172], [228,121,121], [239,0,0], [198,0,0], [158,0,0]], 255.)
                        labelName = 'coefficient of variation [-]'
                    cmap = colors.ListedColormap(colors1)
                    norm = colors.BoundaryNorm(bounds, cmap.N)

                    if projection == 'Basemap_PlateCarree':
                        ims = bm1.imshow(flipud(s2nSrc[ifuture,iscn,11:293,:]), norm=norm, cmap=cmap, interpolation="nearest")
                        bm1.drawcoastlines(linewidth=coastal_line_width, color=coastal_line_color)
                    elif projection == 'Basemap_Robinson':
                        ims = bm1.imshow(bm1.transform_scalar(flipud(s2nSrc[ifuture,iscn]), lons, lats, nx, ny), norm=norm, cmap=cmap, interpolation="nearest")
                        bm1.drawcoastlines(linewidth=coastal_line_width, color=coastal_line_color)
                    elif 'Cartopy' in projection:
                        ##img_extent = [-180, 180, -60, 90]
                        ##im = ax.imshow(src, origin='upper', vmin=int(vmin), vmax=int(vmax), norm=norm, cmap=cmap, transform=original_projection, extent=img_extent)
                        #ims = ax.contourf(lons, lats, s2nSrc[ifuture, iscn,:293], levels=bounds, norm=norm, cmap=cmap, transform=original_projection)
                        ims = ax.contourf(lons, lats, s2nSrc[ifuture, iscn,:293], levels=bounds, cmap=cmap, transform=original_projection)
                        ax.add_feature(cfeature.COASTLINE, linewidth=coastal_line_width, edgecolor=coastal_line_color)

                    ax.text(0., 1.,  '({})'.format(tag),         va="bottom", ha="left",   fontsize=8, transform=ax.transAxes)
                    ax.text(0.5, 1., '{} {}'.format(scn,future), va="bottom", ha="center", fontsize=8, transform=ax.transAxes)

                    if tag == 'c':
                        ax_pos = ax.get_position()
                        ax11 = fig.add_axes([ax_pos.x0+0.008, ax_pos.y0+0.1,  0.115, 0.03])
                        cb1 = mpl.colorbar.ColorbarBase(ax11,
                                                        cmap=cmap,
                                                        norm=norm,
                                                        boundaries=bounds,
                                                        #spacing='proportional',
                                                        spacing='uniform',
                                                        orientation='horizontal')
                        #for t in cb1.ax11.get_yticklabels(): t.set_fontsize(6)
                        cb1.set_ticks(bounds)
                        cb1.set_ticklabels([str(i) for i in bounds])
                        cb1.ax.tick_params(labelsize=5, width=0.25, direction='in')
                        cb1.set_label(labelName, fontsize=7.5)#, labelpad=-0.6)
                        cb1.outline.set_visible(False)

                        """
                        # add information
                        if srcType == 'abschange' or srcType == 'percentchange':
                            ax.text(0.38, 0.01, 
                                    '{}\n{}, {},  {}-historical,  Q{}, Len{}'.format(srcType, scn, scn, prd, Q, Len),
                                    va='bottom', ha='left', fontsize=5, transform=ax.transAxes)
                        elif srcType == 'absolute':
                            ax.text(0.38, 0.01, 
                                    '{}\n{}, {},  {},  Q{}, Len{}'.format(srcType, scn, scn, prd, Q, Len),
                                    va='bottom', ha='left', fontsize=5, transform=ax.transAxes)
                        """

                for suffix in suffixes:
                    outFile = 'globMap.S2N.Ensemble{}.{}.{}.{}.Q{:02}win{:02}_Len{:03}tau{}_{}.{}'.format(mmeType, srcType, prd, season,
                                                                                                          Q, win, Len, tau, index, suffix)
                    outDir = os.path.join(fig_directory_main, season, suffix)
                    #outDir = os.path.join(fig_directory_main, '{}.{}_{}'.format(qvalType,syear_ref,eyear_ref), 
                    #                      '{}-{}.{}-{}.{}-{}'.format(syear_hi, eyear_hi, syear_mf, eyear_mf, syear_ff, eyear_ff),
                    #                      )
                    if not os.path.isdir(outDir): os.makedirs(outDir)
                    outPath = os.path.join(outDir, outFile)

                    if os.path.exists(outPath): 
                        print("File exists, will be overwritten. ")
                        #raw_input("Print figure? Press key to continue...")
                    plt.savefig(outPath, dpi=dpi)
                    print('savefig: {}\n'.format(outPath))
                plt.close(3)


            if GHMMap:
                #==============#
                # 4. GHM table #
                #==============#
                print('<<<<<<<<<<>>>>>>>>>>\n   GHM table   \n<<<<<<<<<<>>>>>>>>>>')
                fig = plt.figure(num=4, figsize=(len(ghms)+1, 4))
                #hold(True)
        
                print('GridSpec...')
                gs = gridspec.GridSpec(len(ghms)+1, 4)  # (rows,cols)
                gs.update(left=0.03, right=0.99, bottom=0.02, top=0.96, hspace=0.005, wspace=0.0025)
                axes = []
                for ighm in range(len(ghms)):
                    if 'Basemap' in projection:
                        for i in range(4):
                            axes.append(plt.subplot(gs[ighm,i]))
                    elif 'Cartopy' in projection:
                        for i in range(4):
                            axes.append(plt.subplot(gs[ighm,i], projection=axis_projection))

                char_list = [chr(ichr) for ichr in range(97, 97+4*len(ghms))]

                for g, ghm in enumerate(ghms):
                    #for s, ssn in enumerate(SSNs):
                    for istep, ((ifuture, future), (iscn, (scn, soc, co2))) in enumerate(itertools.product(enumerate(futures), enumerate(SIMs))):  #[mid, late],[rcp26,rcp85]
                        gs = g*4+istep
                        print('{} {} {} (ax{})...'.format(ghm, future, scn, istep))

                        if prd == 'hist' and srcType == 'absolute': 
                            aSrc = ghmSRC[iscn,g][11:293]
                        elif prd == 'future':
                            aSrc = ghmSRC[ifuture,iscn,g][11:293]

                        ax = axes[gs]
                        plt.sca(ax)
                        if 'Basemap' in projection:
                            ax.axis('off')
                        elif 'Cartopy' in projection:
                            ax.set_extent([-180, 180, -60, 90], original_projection)
                            ax.outline_patch.set_linewidth(0.1)
                            ax.outline_patch.set_edgecolor('#b0b0b0')
                            ax.add_feature(cfeature.LAND, facecolor='#fdfdfd')
                        ax_pos = ax.get_position()

                        if srcType == 'percentchange' or srcType == 'abschange':

                            ks_mask = dictKSmask[scn][future][ghm][11:293]
                            aSrc = ma.masked_array(ma.masked_equal(aSrc,0), mask=ks_mask)

                            agreement = agreementgcm[ifuture,iscn,g][11:293]
                            agreement = ma.masked_equal(agreement, -999)

                            mask1  = ma.make_mask(agreement<agreeThrsh[0])
                            mask21 = ma.make_mask(agreement>=agreeThrsh[0])
                            mask22 = ma.make_mask(agreement<agreeThrsh[1])
                            mask2  = ma.mask_or(mask21, mask22)
                            mask31 = ma.make_mask(agreement>=agreeThrsh[1])
                            mask32 = ma.make_mask(agreement<agreeThrsh[2])
                            mask3  = ma.mask_or(mask31, mask32)
                            mask4  = ma.make_mask(agreement>=agreeThrsh[2])
                            signal1 = ma.masked_array(aSrc, mask=mask1)
                            signal2 = ma.masked_array(aSrc, mask=mask2)
                            signal3 = ma.masked_array(aSrc, mask=mask3)
                            signal4 = ma.masked_array(aSrc, mask=mask4)

                            ### Map   (reference: Kaye et al., Fig. 7b)
                            if srcType == 'percentchange': bounds = [-100, -50, -10, 0, 10, 50, 100, 200]
                            elif srcType == 'abschange': bounds = [-20, -10, -5, 0, 5, 10, 20, 30]
                            ### blue >> green >> yellow >> orange >> red (5)
                            colors1 = divide([[ 23, 23,171],[ 39,161,242],[127,  0,255],[  0,132,132],[240,210,  0],[230,120,  0],[170,  0,  0]],255.)
                            colors2 = divide([[102,109,206],[125,202,253],[175, 94,255],[ 84,172,172],[244,224, 84],[238,164, 84],[198, 84, 84]],255.)
                            colors3 = divide([[195,195,240],[164,219,255],[215,175,255],[168,213,213],[249,239,168],[246,209,168],[226,168,168]],255.)
                            colors4 = divide([[210,210,210],[210,210,210],[210,210,210],[210,210,210],[210,210,210],[210,210,210],[210,210,210]],255.)  # gray
                            cmap1 = mpl.colors.ListedColormap(colors1)
                            cmap2 = mpl.colors.ListedColormap(colors2)
                            cmap3 = mpl.colors.ListedColormap(colors3)
                            cmap4 = mpl.colors.ListedColormap(colors4)
                            norm1 = mpl.colors.BoundaryNorm(bounds, cmap1.N)
                            norm2 = mpl.colors.BoundaryNorm(bounds, cmap2.N)
                            norm3 = mpl.colors.BoundaryNorm(bounds, cmap3.N)
                            norm4 = mpl.colors.BoundaryNorm(bounds, cmap4.N)

                            if 'Basemap' in projection:
                                ims3 = bm2.imshow(flipud(signal3), norm=norm3, cmap=cmap3, interpolation="nearest")
                                ims2 = bm2.imshow(flipud(signal2), norm=norm2, cmap=cmap2, interpolation="nearest")
                                ims1 = bm2.imshow(flipud(signal1), norm=norm1, cmap=cmap1, interpolation="nearest")
                                ims4 = bm2.imshow(flipud(signal4), norm=norm4, cmap=cmap4, interpolation="nearest")
                                bm2.drawcoastlines(linewidth=0.05, color=coastal_line_color)
                            elif 'Cartopy' in projection:
                                ##img_extent = [-180, 180, -60, 90]
                                ##im = ax.imshow(src, origin='upper', vmin=int(vmin), vmax=int(vmax), norm=norm, cmap=cmap, transform=ccrs.PlateCarree(), extent=img_extent)
                                im1 = plt.contourf(lons, lats, flipud(signal1[:293]), levels=bounds, norm=norm1, cmap=cmap1, transform=original_projection)
                                im2 = plt.contourf(lons, lats, flipud(signal2[:293]), levels=bounds, norm=norm2, cmap=cmap2, transform=original_projection)
                                im3 = plt.contourf(lons, lats, flipud(signal3[:293]), levels=bounds, norm=norm3, cmap=cmap3, transform=original_projection)
                                im4 = plt.contourf(lons, lats, flipud(signal4[:293]), levels=bounds, norm=norm4, cmap=cmap4, transform=original_projection)
                                ax.add_feature(cfeature.COASTLINE, linewidth=coastal_line_width, edgecolor=coastal_line_color)

                        else:   # historical absolute or change
                            #if srcType == 'abschange':
                            #    ks_mask = dictKSmask[prd][ssn][ghm][11:293]
                            #    aSrc = ma.masked_array(ma.masked_equal(aSrc,0), mask=ks_mask)
                            #    bounds = [-60,-30,0,30,60]
                            #    cmap = cm.bwr
                            if srcType == 'absolute':
                                bounds = [0,10,20,30,40,50,]
                                cmap = cm.jet
                            norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

                            if 'Basemap' in projection:
                                ims1 = bm2.imshow(flipud(aSrc), norm=norm, cmap=cmap, interpolation='nearest')
                                bm2.drawcoastlines(linewidth=0.05, color=coastal_line_color)
                            elif 'Cartopy' in projection:
                                #ims1 = plt.contourf(lons, lats, flipud(aSrc[:293,:]), levels=bounds, norm=norm, cmap=cmap, transform=original_projection)
                                ims1 = plt.contourf(lons, lats, flipud(aSrc[:293,:]), levels=bounds, cmap=cmap, transform=original_projection)
                                ax.add_feature(cfeature.COASTLINE, linewidth=coastal_line_width, edgecolor=coastal_line_color)

                        ax.set_xticks([])
                        ax.set_yticks([])
                        if ghm == ghms[0]: ax.text(0.5, 1.05, '{}-{}'.format(scn,future), va='bottom', ha='center', fontsize=8, transform=ax.transAxes)
                        if future == futures[0] and scn == SIMs[0][0]: ax.text(-0.05, 0.5, ghm, va='center', ha='center', fontsize=8, transform=ax.transAxes, rotation='vertical',)

                # Colorbar
                if srcType == 'percentchange' or srcType == 'abschange':

                    ### 2D colorbar for Kaye et al.-plot:
                    ax_pos = axes[-1].get_position()
                    ax11 = fig.add_axes([ax_pos.x0+0.10, ax_pos.y0-0.090, 0.12, 0.02])
                    ax12 = fig.add_axes([ax_pos.x0+0.10, ax_pos.y0-0.070, 0.12, 0.02])
                    ax13 = fig.add_axes([ax_pos.x0+0.10, ax_pos.y0-0.050, 0.12, 0.02])
                    ax14 = fig.add_axes([ax_pos.x0+0.10, ax_pos.y0-0.030, 0.12, 0.02])
                    cmap = [cmap1, cmap2, cmap3, cmap4]
                    
                    for i, axs in enumerate([ax12, ax13, ax14]): 
                        norm = mpl.colors.BoundaryNorm(bounds, cmap[i+1].N)
                        cb = mpl.colorbar.ColorbarBase(axs, cmap=cmap[i+1], norm=norm,
                                                       boundaries=bounds,
                                                       #extend='neither',
                                                       #extendfrac='auto',
                                                       #spacing='proportional',
                                                       spacing='uniform',
                                                       orientation='horizontal')
                        cb.set_ticks(bounds)
                        cb.set_ticklabels([])
                        cb.ax.tick_params(labelsize=5, width=0.25, direction='in')
                        cb.outline.set_visible(False)
                    norm = mpl.colors.BoundaryNorm(bounds, cmap[0].N)
                    cb1 = mpl.colorbar.ColorbarBase(ax11, cmap=cmap[0], norm=norm,
                                                    boundaries=bounds,
                                                    #extend='neither',
                                                    #extendfrac='auto',
                                                    #spacing='proportional',
                                                    spacing='uniform',
                                                    orientation='horizontal')
                    #for t in cb1.ax11.get_yticklabels(): t.set_fontsize(6)
                    cb1.ax.set_xticks(bounds)
                    cb1.ax.set_xticklabels([str(int(ibound)) for ibound in bounds], rotation=-45, fontsize=4.5)
                    cb1.ax.tick_params(labelsize=5, width=0.25, direction='in')
                    cb1.outline.set_visible(False)
                    if srcType == 'percentchange':
                        cb1.set_label('relative change [%]', fontsize=4.5)  #, labelpad=-0.9)
                    elif srcType == 'abschange':
                        cb1.set_label('{} [{}]'.format(srcType, dictUnit[index]), fontsize=4.5)  # labelpad=-0.6)
                    fig.text(ax_pos.x0+0.09,  ax_pos.y0-0.032, str(int(agreeThrsh[2]*1e2)), va='center', ha='center', fontsize=4.5)
                    fig.text(ax_pos.x0+0.09,  ax_pos.y0-0.052, str(int(agreeThrsh[1]*1e2)), va='center', ha='center', fontsize=4.5)
                    fig.text(ax_pos.x0+0.09,  ax_pos.y0-0.072, str(int(agreeThrsh[0]*1e2)), va='center', ha='center', fontsize=4.5)
                    fig.text(ax_pos.x0+0.075, ax_pos.y0-0.065, 'agreement [%]',             va='center', ha='center',  rotation='vertical', fontsize=4.5)

                else:
                    ax_pos = axes[-1].get_position()
                    ax11 = fig.add_axes([ax_pos.x0+0.06, ax_pos.y0-0.080,  0.15, 0.03])
                    cb1 = mpl.colorbar.ColorbarBase(ax11, cmap=cmap, norm=norm, orientation='horizontal')
                    cb1.set_ticks(bounds)
                    cb1.set_ticklabels([str(int(ibound)) for ibound in bounds])
                    cb1.ax.tick_params(labelsize=5, width=0.25, direction='in')
                    cb1.outline.set_visible(False)
                    #if srcType == 'abschange':
                    #    cb1.set_label('{} [{}]'.format(srcType, dictUnit[index]), fontsize=5, labelpad=-0.6)
                    if srcType == 'absolute':
                        cb1.set_label('[{}]'.format(dictUnit[index]), fontsize=5, labelpad=-0.6)

                for suffix in suffixes:
                    outFile = 'globMap.GHM.{}.{}.{}.Q{:02}win{:02}_Len{:03}tau{}_{}.{}'.format(prd, srcType, season, Q, win, Len, tau, index, suffix)
                    outDir = os.path.join(fig_directory_main, season, suffix)
                    #outDir = os.path.join(fig_directory_main, '{}.{}_{}'.format(qvalType,syear_ref,eyear_ref), 
                    #                      '{}-{}.{}-{}.{}-{}'.format(syear_hi, eyear_hi, syear_mf, eyear_mf, syear_ff, eyear_ff),
                    #                      )
                    if not os.path.isdir(outDir): os.makedirs(outDir)
                    outPath = os.path.join(outDir,outFile)

                    if os.path.exists(outPath): 
                        print("File exists, will be overwritten.")
                        #raw_input("Print figure? Press key to continue...")
                    plt.savefig(outPath, dpi=dpi)
                    print('savefig: {}\n'.format(outPath))
                plt.close(4)

            if UncertaintySource and uncertainty_source_rate is not None:
                #=======================#
                # 5. Uncertainty source #
                #=======================#
                
                if prd == 'hist':
                    aSrc = uncertainty_source_rate[:,11:293,:] * 100
                elif prd == 'future':
                    aSrc = uncertainty_source_rate[:,:,11:293,:] * 100

                bounds = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
                norm = colors.Normalize(vmin=0, vmax=100)

                if 'Basemap' in projectioan:
                    if projectioan == 'Basemap_PlateCarree':
                        figsize = (8, 3.8)
                        left, right, bottom, top, hspace, wspace =0.007, 0.995, 0.01, 0.975, 0.01, 0.01
                    elif projection == 'Basemap_Robinson':
                        dx, dx2, dx3, width = 0.028, 0.015, 0.001, 0.1
                        figsize = (8, 4)
                        left, right, bottom, top, hspace, wspace = 0.01, 0.995, 0.02, 0.95, 0.005, 0.002
                    fig = plt.figure(num=5, figsize=figsize)
                    gs = gridspec.GridSpec(2, 2)  # (rows,cols)
                    gs.update(left=left, right=right, bottom=bottom, top=top, hspace=hspace, wspace=wspace)
                    ax1 = plt.subplot(gs[0,0])
                    ax2 = plt.subplot(gs[0,1])
                    ax3 = plt.subplot(gs[1,0])
                    ax4 = plt.subplot(gs[1,1])
                elif 'Cartopy' in projection:
                    fig = plt.figure(num=5, figsize=(8, 3.5))
                    gs = gridspec.GridSpec(2, 2)  # (rows,cols)
                    gs.update(left=0.004, right=0.99, bottom=0.02, top=0.95, hspace=0.005, wspace=0.001)
                    ax1 = plt.subplot(gs[0,0], projection=axis_projection)
                    ax2 = plt.subplot(gs[0,1], projection=axis_projection)
                    ax3 = plt.subplot(gs[1,0], projection=axis_projection)
                    ax4 = plt.subplot(gs[1,1], projection=axis_projection)
                axes = [ax1, ax2, ax3, ax4]

                #for i, (ssn, ax, tag) in enumerate(zip(SSNs, [ax1,ax2,ax3,ax4], ['a', 'b', 'c', 'd'])):
                istep = 0
                for (ifuture, future), (iscn, (scn, soc, co2)) in itertools.product(enumerate(futures), enumerate(SIMs)):  # [mid, late], [rcp26, rcp85]
                    if (ifuture == 0 and iscn == 0): istep = istep + 1
                    ax = axes[istep]
                    tag = ['a', 'b', 'c', 'd'][istep]

                    print('ax{}..'.format(istep+1))
                    plt.sca(ax)
                    #ax.axis('off')
                    ax.set_extent([-180, 180, -60, 90], ccrs.PlateCarree())
                    ax.outline_patch.set_linewidth(0.1)
                    ax.outline_patch.set_edgecolor('#b0b0b0')
                    ax.add_feature(cfeature.LAND, facecolor='#fdfdfd')

                    if prd == 'hist':
                        _aSrc = aSrc[iscn]
                    elif prd == 'future':
                        _aSrc = aSrc[ifuture,iscn]
                    #ax.set_title('Ratio of GCM variation to total variation', fontsize=8)
                    if projection == 'Basemap_PlateCarree':
                        ims1 = bm1.imshow(flipud(_aSrc), norm=norm, cmap=cm.RdYlBu, interpolation="nearest")
                        bm1.drawcoastlines(linewidth=coastal_line_width, color=coastal_line_color)
                    elif projection == 'Basemap_Robinson':
                        ims1 = bm1.imshow(bm1.transform_scalar(flipud(_aSrc), lons, lats, nx, ny), norm=norm, cmap=cm.RdYlBu, interpolation="nearest")
                        bm1.drawcoastlines(linewidth=coastal_line_width, color=coastal_line_color)
                    elif 'Cartopy' in projection:
                        #im = plt.contourf(lons, lats, _aSrc[:293], levels=bounds, norm=norm, cmap=cm.RdYlBu, transform=original_projection)
                        im = plt.contourf(lons, lats, _aSrc[:293], levels=bounds, cmap=cm.RdYlBu, transform=original_projection)
                        ax.add_feature(cfeature.COASTLINE, linewidth=coastal_line_width, edgecolor=coastal_line_color)
                    ax.text(0., 1., '({})'.format(tag),         va="bottom", ha="left",   fontsize=8, transform=ax.transAxes)
                    ax.text(0.5, 1., '{} {}'.format(scn,future), va="bottom", ha="center", fontsize=8, transform=ax.transAxes)

                    if tag == 'd':
                        ax_pos = ax.get_position()
                        ax2 = fig.add_axes([ax_pos.x0+0.2, ax_pos.y0+0.055,  0.225, 0.02])
                        cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=cm.RdYlBu, norm=norm, orientation='horizontal')
                        cb1.set_ticks(bounds)
                        cb1.set_xticklabels([str(int(ibound)) for ibound in bounds])
                        cb1.ax.tick_params(labelsize=5, width=0.25, direction='in')
                        cb1.outline.set_visible(False)
                        cb1.set_label('[%]', fontsize=5)#, labelpad=-0.6)

                for suffix in suffixes:
                    outFile = 'globMap.UNCSRC.{}.{}.{}.Q{:02}win{:02}_Len{:03}tau{}_{}.{}'.format(prd, srcType, season, Q, win, Len, tau, index, suffix)
                    outDir = os.path.join(fig_directory_main, season, suffix)
                    #outDir = os.path.join(fig_directory_main, '{}.{}_{}'.format(qvalType,syear_ref,eyear_ref), 
                    #                      '{}-{}.{}-{}.{}-{}'.format(syear_hi, eyear_hi, syear_mf, eyear_mf, syear_ff, eyear_ff),
                    #                      )
                    if not os.path.isdir(outDir): os.makedirs(outDir)
                    outPath = os.path.join(outDir,outFile)

                    if os.path.exists(outPath): 
                        print("File exists, will be overwritten.")
                        #raw_input("Print figure? Press key to continue...")
                    plt.savefig(outPath, dpi=dpi)
                    print('savefig: {}\n'.format(outPath))
                plt.close(5)

        endTime = datetime.datetime.now()
        diffTime = endTime - strTime
        print('end @{}'.format(endTime.strftime("%Y-%m-%d %H:%M:%S")))
        print('took {} min in total.'.format(int(diffTime.seconds/60)))

    return

if __name__=='__main__':
    main(*sys.argv)


