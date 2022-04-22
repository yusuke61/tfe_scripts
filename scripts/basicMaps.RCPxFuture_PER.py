#!/usr/bin/python
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
from mpl_toolkits.basemap import Basemap
from utiltools import fillinMasked, extractMasked, flipud
from scipy import stats
#------------------------------------------------------------------------------------------------------
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = 'Arial'
hostname = os.uname()[1]
TEST = False

#mmeTYPs = ['median','mean']  
mmeTYPs = ['median']  
#mmeTYPs = ['mean']  

#variables = ['pr', 'evap', 'qtot', 'soilmoist']
#variables = ['pr', 'evap', 'swe']#, 'qtot']
variables = ['pr', 'evap']#, 'qtot']
#variables = ['swe']

_ghms = ['cwatm', 'h08', 'lpjml', 'matsiro', 'watergap2']

gcms = ['hadgem2-es', 'ipsl-cm5a-lr', 'gfdl-esm2m', 'miroc5']


# target simlation type
SIMs  = [
    ('rcp26', '2005soc'),
    ('rcp85', '2005soc'),
    #('rcp26', 'nosoc'),
    #('rcp85', 'nosoc'),
    ]

# a whole target period of drought analysis
syear, eyear = 1861, 2099
syear_hist, eyear_hist = 1861, 2005
climate_window = 30
futures = ['mid', 'far']
syear_nf, syear_ff = 2036, 2070
eyear_nf, eyear_ff = syear_nf+climate_window-1, syear_ff+climate_window-1

EnsmMap = True
S2NMap  = False
GHMMap  = False
UncertaintySource = False
DesertMASK = False

#KStest = True
KStest = False
ks = 0.05    # 95% level significance

#agreeThrsh = [0.8, 0.7, 0.6]
#agreeThrsh = [0.8, 0.6]
agreeThrsh = [0.6]

dict_variable = {
    'pr':   'a',
    'evap': 'b',
    'swe':  'c',
    }

projection = 'Basemap_PlateCarree'

#coastal_line_color = '#808080'
coastal_line_color = '#3c3c3c'

coastal_line_width = 0.1

suffixes = ['png', 'pdf', 'eps']

if TEST:
    #variables = ['pr', 'evap']
    variables = ['evap']
    _ghms = ['matsiro']
    gcms = ['hadgem2-es']
    #ghms = ['cwatm', 'matsiro']
    #gcms = ['hadgem2-es', 'ipsl-cm5a-lr']
    #KStest = False
    #SIMs = [('rcp85', '2005soc')]


#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
today = datetime.date.today().strftime('%Y%m%d')

### set input data directory: 
if 'scs' in hostname  or hostname == 'sun220s.nies.go.jp':
    if   'scs' in hostname:                data_directory = '/data/rg001/sgec0017/data'
    elif hostname == 'sun220s.nies.go.jp': data_directory = '/sraid02/satoh/data'
    elif hostname == 'sun90h.nies.go.jp':  data_directory = '/sraid02/satoh/data'
    drghtDIR = os.path.join(data_directory, 'isimip2b.drought')
    if KStest: fig_directory_main = os.path.join(data_directory, 'figure_box', 'drought_tpcd_isimip2b', 'basicMaps.RCPxFuture_PER_withKStest')  #, today)
    else:      fig_directory_main = os.path.join(data_directory, 'figure_box', 'drought_tpcd_isimip2b', 'basicMaps.RCPxFuture_PER_withoutKStest')  #, today)
    if TEST: fig_directory_main = os.path.join(data_directory, 'figure_box', 'drought_tpcd_isimip2b', 'basicMaps.RCPxFuture_PER_TEST')
    mapmaskDIR = os.path.join(data_directory, 'mapmask')
lndmskPath = os.path.join(mapmaskDIR, 'ISIMIP2b_landseamask_generic.nc4')
grlmskPath = os.path.join(mapmaskDIR, 'GAUL/flt/gaul2014_05deg.flt')      # GreenLand is 98
dstmskPath = os.path.join(mapmaskDIR, 'grid_index.nc4')
grdaraPath = os.path.join(mapmaskDIR, 'grd_ara.hlf')

gridarea = np.fromfile(grdaraPath, 'float32').byteswap().reshape(360,720)


dictUnit = {
    'pr': 'kg/m2/s',
    'evap': 'kg/m2/s',
    'qtot': 'kg/m2/s',
    'soilmoist': 'kg/m2',
    }

dictVlim = {
    'pr':        {'absolute': [None, None], 'abschange': [None, None], 'percentchange': [None, None]},
    'evap':      {'absolute': [None, None], 'abschange': [None, None], 'percentchange': [None, None]},
    'qtot':      {'absolute': [None, None], 'abschange': [None, None], 'percentchange': [None, None]},
    'soilmoist': {'absolute': [None, None], 'abschange': [None, None], 'percentchange': [None, None]},
    }

dictSeason = {
    'ALL': 365,
    'DJF':  90,
    'MAM':  92,
    'JJA':  92,
    'SON':  91,
    }


years = range(syear, eyear+1)
years_hist = range(syear_hist, eyear_hist+1)
nyears_hist = len(years_hist)
nyears_future = climate_window
histStr, histEnd = years.index(syear_hist), years.index(eyear_hist)  # 1971-2004 (34)
nfStr, nfEnd = years.index(syear_nf), years.index(eyear_nf)          # 2036-2065 (30)
ffStr, ffEnd = years.index(syear_ff), years.index(eyear_ff)          # 2070-2099 (30)

grl_mask = ma.masked_equal(fromfile(grlmskPath,'float32').reshape(360,720), 98).mask
dst_mask = ma.masked_equal(Dataset(dstmskPath).variables['grid_index'][:], 11).mask
lnd_mask = ma.masked_not_equal(Dataset(lndmskPath).variables['LSM'][:][0], 1.0).mask
lnd_mask = ma.mask_or(lnd_mask, grl_mask)
if DesertMASK: lnd_mask = ma.mask_or(lnd_mask, dst_mask)
nlndgrid = 360*720 - lnd_mask.sum()
total_gridarea = ma.masked_array(gridarea, mask=Dataset(lndmskPath)['LSM'][0].mask).sum()

#resolution = 'i'
resolution = 'l'
bm1 = Basemap(projection='cyl',llcrnrlat=-56.5,urcrnrlat=84.5,llcrnrlon=-180.,urcrnrlon=180.,resolution=resolution) # correct
#bm2 = Basemap(projection='cyl',llcrnrlat=-90.0,urcrnrlat=90.0,llcrnrlon=-180.,urcrnrlon=180.,resolution=resolution) # correct
bm2 = Basemap(projection='cyl',llcrnrlat=-56.5,urcrnrlat=84.5,llcrnrlon=-180.,urcrnrlon=180.,resolution=resolution) # correct

dpi = 300


# ----------------------------------------------------------------------------------------------------------------------
def read_nc_core(variable, soc, scn, ghm, gcm, period):
    # Read only Land Data

    if variable == 'pr':
        input_directory = os.path.join(data_directory, 'isimip2b.standardized_drought', 'preprocess_0', 'pr')
        if period == 'historical': srcFile = '{}_historical_pr_monthly_1861_2005.nc4'.format(gcm)
        else:                      srcFile = '{}_{}_pr_monthly_2006_2099.nc4'.format(gcm, scn)
    else:
        input_directory_main = os.path.join(data_directory, 'isimip2b', 'out', 'nc', 'water_global')
        if period == 'historical': 
            srcFile = '{}_{}_ewembi_historical_histsoc_co2_{}_global_monthly_1861_2005.nc4'.format(ghm, gcm, variable)
        else:                      
            period = 'future'
            srcFile = '{}_{}_ewembi_{}_{}_co2_{}_global_monthly_2006_2099.nc4'.format(ghm, gcm, scn, soc, variable)
        input_directory = os.path.join(input_directory_main, ghm, gcm, period)

    srcPath = os.path.join(input_directory, srcFile)

    if not os.path.isfile(srcPath): 
        print('Error!! The file NOT exist... Check!! {}'.format(srcPath))
        sys.exit()
    else:
        #print('{} >> original shape {}'.format(srcPath, Dataset(srcPath).variables[variable][:].shape))
        aSrc = extractMasked(Dataset(srcPath).variables[variable][:].reshape(-1,12,360,720), lnd_mask)  # (nyear, 12, nland)
        aSrc = aSrc.mean(axis=1)                           # (nyear, nland)
        print('Read: {} {}'.format(srcPath, aSrc.shape))   # (nyear, nland)
        return aSrc


# ----------------------------------------------------------------------------------------------------------------------
def read_nc(variable, soc, scn, ghm, gcm):
    src =  np.concatenate([read_nc_core(variable, soc, scn, ghm, gcm, period) for period in ['historical', scn]], axis=0)  # (nYEAR, nland) nYEAR = 2099-1861+1 = 239
    print(src.shape, src.min(), src.max())
    return src


# ----------------------------------------------------------------------------------------------------------------------
def change_arearate(target_type, src, intensity=None):
    """
    To estimate the rate of certain change (decrease/increase)
    :param target_type: (str) decrease, increase, significant increase
    :param src: (array) percentchange or absolutechange
    :return: gridrate
    """

    if target_type == 'decrease':
        target_mask = ma.make_mask(0<=src)
    elif target_type == 'increase':
        target_mask = ma.make_mask(src<=0)
    elif target_type == 'significant_increase':
        target_mask = ma.make_mask(src<=intensity)
    else:
        print('Something is wrong. Check target_type!!')

    target_gridarea = ma.masked_array(gridarea, mask=target_mask).sum()
    print('change_arearate: ({}) {} / {} = {}'.format(target_type, target_gridarea, total_gridarea, target_gridarea / total_gridarea))

    return target_gridarea / total_gridarea


# ----------------------------------------------------------------------------------------------------------------------
def main(*args):
    strTIME = datetime.datetime.now()
    print('START detect.drought.py')

    for mmeType, variable in itertools.product(mmeTYPs, variables):

        strTime = datetime.datetime.now()
        print('\nJob started !!!  {} {}'.format(mmeType, variable))

        if variable == 'pr':
            ghms = ['gcm_output']
        else:
            ghms = _ghms.copy()
        nGHM, nGCM = len(ghms), len(gcms)
        nMME = nGHM*nGCM

        # Reading data
        print('reading data...')
        aSRC = array([[[read_nc(variable, soc, scn, ghm, gcm) for ghm in ghms] for gcm in gcms] for (scn, soc) in SIMs])
        aSRC = ma.masked_equal(aSRC,1e+20)                                          # orignally, missing_value in each data is 1e+20
        print('aSRC.shape {}'.format(aSRC.shape))                                   #(nSCN,nGCM,nGHM,nYear,nLand))

        # Cutout specific periods
        aHIS = aSRC[:,:,:,histStr:histEnd+1,:]                                      #      (nSIM,nGCM,nGHM,30,nLand)
        aFTR = array([aSRC[:,:,:,nfStr:nfEnd+1,:], aSRC[:,:,:,ffStr:ffEnd+1,:]])    # (nPRD,nSIM,nGCM,nGHM,30,nLand)
        #print('aHIS {}, aFTR {}'.format(aHIS.shape, aFTR.shape))
        del aSRC

        ### Kolmogorov-Smirnov test and dictKSmask
        ### if p-value is less than ks, change in the grid is significant.
        print('KS testing...')
        dictKSmask = {scn: {'mid': {}, 'far': {}} for scn,soc in SIMs}
        if KStest:
            for s, (scn, soc) in enumerate(SIMs):
                for i, f in enumerate(futures):
                    dictKSmask[scn][f]['all'] = ma.make_mask(
                                                    fillinMasked(
                                                        array([ 0 if stats.ks_2samp(aHis.reshape(-1),aFut.reshape(-1)).pvalue < ks else 1
                                                                for aHis, aFut in zip(aHIS[s].T, aFTR[i,s].T)]),
                                                        lnd_mask)
                                                    == 1) #(nY,nX)
                    print('{} {} all {}  {}'.format(scn, f, dictKSmask[scn][f]['all'].shape, (nlndgrid-dictKSmask[scn][f]['all'].sum())/float(nlndgrid)))
                    for j, ghm in enumerate(ghms):
                        dictKSmask[scn][f][ghm] = ma.make_mask(
                                                    fillinMasked(
                                                        array([ 0 if stats.ks_2samp(aHis.reshape(-1),aFut.reshape(-1)).pvalue < ks else 1
                                                                for aHis, aFut in zip(aHIS[s,:,j,...].T, aFTR[i,s,:,j,...].T)]),
                                                        lnd_mask)
                                                    == 1) #(nY,nX)
                        print('{} {} {} {}  {}'.format(scn, f, ghm, dictKSmask[scn][f][ghm].shape, (nlndgrid-dictKSmask[scn][f][ghm].sum())/float(nlndgrid)))
        else:
            print('KStest is False. Skep it.')
            for i, f in enumerate(futures):
                for s, (scn, soc) in enumerate(SIMs):
                    dictKSmask[scn][f]['all'] = zeros((360,720),'bool')
                    for j, ghm in enumerate(ghms):
                        dictKSmask[scn][f][ghm] = zeros((360,720),'bool')

        print('genarating values...')
        # Make ensemble value for each combination of GHM&GCM
        clmHIS = aHIS.mean(axis=3)                           # (     nSIM,nGCM,nGHM,nLand)
        clmFTR = aFTR.mean(axis=4)                           # (nPRD,nSIM,nGCM,nGHM,nLand)
        pc_clmFTR= (clmFTR-clmHIS)/clmHIS*100                    # (nPRD,nSIM,nGCM,nGHM,nLand)

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
        pc_allMME = chng_allMME / allMMEhis * 100             # (nPRD,nSIM,         nLand) >> Fig [%]
        pc_gcmMME = chng_gcmMME / gcmMMEhis * 100             # (nPRD,nSIM,    nGHM,nLand) >> Fig [%]

        ### Spread(spr) among ensemble samples
        if mmeType == 'median':   ### get inter quartile range (IQR)
            spr_allMMEhis = subtract(*percentile(clmHIS,    [75,25], axis=(1,2)))  #      (nSIM,nLand)
            spr_pc_allMME = subtract(*percentile(pc_clmFTR, [75,25], axis=(2,3)))  # (nPRD,nSIM,nLand)
        elif mmeType == 'mean':   # get standard deviation (std)
            spr_allMMEhis = clmHIS.std(axis=(1,2))                                 #      (nSIM,nLand)
            spr_pc_allMME = pc_clmFTR.std(axis=(2,3))                              # (nPRD,nSIM,nLand)

        ### Uncertainty  (Signal to noise)  TODO!!
        if mmeType == 'median':   ### Singal to noise (s2n)   (ref. Prudhomme et al. 2014)
            s2n_allMMEhis = allMMEhis / spr_allMMEhis            #      (nSIM,nLand) >> Fig
            s2n_pc_allMME = pc_allMME / spr_pc_allMME            # (nPRD,nSIM,nLand) >> Fig
        elif mmeType == 'mean':   ### Coefficient of variation (Cov)
            s2n_allMMEhis = spr_allMMEhis / allMMEhis            #      (nSIM,nLand) >> Fig
            s2n_pc_allMME = spr_pc_allMME / pc_allMME            # (nPRD,nSIM,nLand) >> Fig

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
        ratGCMvar_clmHIS    = var(clmHIS,axis=1).mean(axis=1) / var(clmHIS,axis=(1,2))         #      (nSIM,nLand) >> Fig
        ratGCMvar_pc_allMME = var(pc_clmFTR,axis=2).mean(axis=2) / var(pc_clmFTR,axis=(2,3))   # (nPRD,nSIM,nLand) >> Fig
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
        for prd,               srcType,      (mmeSrc,    agreementall),        s2nSrc,      (ghmSRC, agreementgcm), uncertainty_source_rate in [
            #['hist',        'absolute',   (allMMEhis,         None), s2n_allMMEhis,   (gcmMMEhis,         None),        ratGCMvar_clmHIS],
            #['future',      'absolute',   (allMMEftr,         None),          None,   (gcmMMEftr,         None),                    None],
            #['future',     'abschange', (chng_allMME, agreementALL), s2n_pc_allMME, (chng_gcmMME, agreementGCM),                    None],
            ['future', 'percentchange',   (pc_allMME, agreementALL), s2n_pc_allMME,   (pc_gcmMME, agreementGCM),     ratGCMvar_pc_allMME],
            ]:
            print('\n==================\n{} {}\n=================='.format(prd, srcType))

            mmeSrc = ma.masked_equal(ma.masked_equal(fillinMasked(mmeSrc, lnd_mask), 1e+20), 0)  # his:(nSIM,     nY,nX), future:(nPRD,nSIM,     nY,nX)
            ghmSRC = ma.masked_equal(ma.masked_equal(fillinMasked(ghmSRC, lnd_mask), 1e+20), 0)  # his:(nSIM,nGHM,nY,nX), future:(nPRD,nSIM,nGHM,nY,nX)
            if S2NMap and s2nSrc is not None:
                s2nSrc = ma.masked_equal(fillinMasked(s2nSrc, lnd_mask), 1e+20)                  # his:(nSIM,     nY,nX), future:(nPRD,nSIM,     nY,nX)
            if agreementall is not None: agreementall = ma.masked_equal(fillinMasked(agreementall, lnd_mask), 1e+20)  # (nPRD,nSIM,     nY,nX)
            if agreementgcm is not None: agreementgcm = ma.masked_equal(fillinMasked(agreementgcm, lnd_mask), 1e+20)  # (nPRD,nSIM,nGHM,nY,nX)
            if uncertainty_source_rate is not None: uncertainty_source_rate = ma.masked_equal(fillinMasked(uncertainty_source_rate,lnd_mask),1e+20)  # his:(nSIM,nY,nX),future:(nPRC,nSIM,nY,nX)
            print('mmeSrc.shape', mmeSrc.shape)
            print('ghmSRC.shape', ghmSRC.shape)

            vmin1, vmax1 = percentile(mmeSrc.compressed(), [10,80])
            vmin2, vmax2 = percentile(ghmSRC.compressed(), [10,80]) 

            if EnsmMap:
                #=============#
                # 1. Ensemble #
                #=============#
                print('Just a moment. preparing a figure object... (num=1)')
                fig = plt.figure(num=1, figsize=(8, 3.8))
                #plt.hold(True)
                print('GridSpec...')
                gs = gridspec.GridSpec(2, 2)  # (rows,cols)
                gs.update(left=0.007, right=0.995, bottom=0.01, top=0.975, hspace=0.01, wspace=0.01)

                ax1 = plt.subplot(gs[0,0])  # for mid-26
                ax2 = plt.subplot(gs[0,1])  # for mid-85
                ax3 = plt.subplot(gs[1,0])  # for far-26
                ax4 = plt.subplot(gs[1,1])  # for far-85
                axes = [ax1, ax2, ax3, ax4]

                if 'change' in srcType: changes = ['decrease', 'increase', 'significant_increase']
                else:                   changes = ['decrease', 'increase']
                df = pd.DataFrame(index=changes, columns=['{}-{}'.format(scn, future) for (scn, soc) in SIMs for future in futures])

                for istep, ((ifuture, future), (iscn, (scn, soc))) in enumerate(itertools.product(enumerate(futures), enumerate(SIMs))): # [mid, far]  # [rcp26, rcp85]
                        print(ifuture, iscn, istep)
                        ax = axes[istep]
                        tag = dict_variable[variable]+['1', '2', '3', '4'][istep]

                        print('ax{}..'.format(istep+1))
                        plt.sca(ax)
                        ax.axis('off')
                        ax_pos = ax.get_position()
                        #norm1 = mpl.colors.Normalize(vmin=vmin1, vmax=vmax1)

                        if srcType == 'percentchange' or srcType == 'abschange':

                            aSrc = ma.masked_array(ma.masked_equal(mmeSrc[ifuture, iscn],0), mask=dictKSmask[scn][future]['all'])

                            if len(agreeThrsh) == 3:

                                mask1  = ma.make_mask(agreementall[ifuture, iscn]<agreeThrsh[0])
                                mask21 = ma.make_mask(agreementall[ifuture, iscn]>=agreeThrsh[0])
                                mask22 = ma.make_mask(agreementall[ifuture, iscn]<agreeThrsh[1])
                                mask2  = ma.mask_or(mask21, mask22)
                                mask31 = ma.make_mask(agreementall[ifuture, iscn]>=agreeThrsh[1])
                                mask32 = ma.make_mask(agreementall[ifuture, iscn]<agreeThrsh[2])
                                mask3  = ma.mask_or(mask31, mask32)
                                mask4  = ma.make_mask(agreementall[ifuture, iscn]>=agreeThrsh[2])
                                signal1 = ma.masked_array(aSrc, mask=mask1)
                                signal2 = ma.masked_array(aSrc, mask=mask2)
                                signal3 = ma.masked_array(aSrc, mask=mask3)
                                signal4 = ma.masked_array(aSrc, mask=mask4)

                                ### Map   (reference: Kaye et al., Fig. 7b)
                                if srcType == 'percentchange': 
                                    bounds = [-20, -10, -5, 0, 5, 10, 20]    # Caution! bounds[0] = vmin, bounds[-1]=vmax
                                elif srcType == 'abschange': bounds = [-20, -10, -5, 0, 5, 10, 20]
                                colors1 = divide([[170,  0,  0],[230,120,  0],[240,210,  0],[127,  0,255],[ 39,161,242],[ 23, 23,171]],255.)
                                colors2 = divide([[198, 84, 84],[238,164, 84],[244,224, 84],[175, 94,255],[125,202,253],[102,109,206]],255.)
                                colors3 = divide([[226,168,168],[246,209,168],[249,239,168],[215,175,255],[164,219,255],[195,195,240]],255.)
                                colors4 = divide([[210,210,210],[210,210,210],[210,210,210],[210,210,210],[210,210,210],[210,210,210]],255.)  # gray
                                cmap1 = mpl.colors.ListedColormap(colors1)
                                cmap2 = mpl.colors.ListedColormap(colors2)
                                cmap3 = mpl.colors.ListedColormap(colors3)
                                cmap4 = mpl.colors.ListedColormap(colors4)
                                norm1 = mpl.colors.BoundaryNorm(bounds, cmap1.N)
                                norm2 = mpl.colors.BoundaryNorm(bounds, cmap2.N)
                                norm3 = mpl.colors.BoundaryNorm(bounds, cmap3.N)
                                norm4 = mpl.colors.BoundaryNorm(bounds, cmap4.N)
                                ims3 = bm1.imshow(flipud(signal3[11:293]), norm=norm3, cmap=cmap3, interpolation='nearest')
                                ims2 = bm1.imshow(flipud(signal2[11:293]), norm=norm2, cmap=cmap2, interpolation='nearest')
                                ims1 = bm1.imshow(flipud(signal1[11:293]), norm=norm1, cmap=cmap1, interpolation='nearest')
                                ims4 = bm1.imshow(flipud(signal4[11:293]), norm=norm4, cmap=cmap4, interpolation='nearest')
                                bm1.drawcoastlines(linewidth=0.005, color='#808080')
                                  
                                if tag == dict_variable[variable]+'3':
                                    ### 2D colorbar for Kaye et al.-plot:
                                    ax11 = fig.add_axes([ax_pos.x0+0.028, ax_pos.y0+0.080, 0.10, 0.02])
                                    ax12 = fig.add_axes([ax_pos.x0+0.028, ax_pos.y0+0.100, 0.10, 0.02])
                                    ax13 = fig.add_axes([ax_pos.x0+0.028, ax_pos.y0+0.120, 0.10, 0.02])
                                    ax14 = fig.add_axes([ax_pos.x0+0.028, ax_pos.y0+0.140, 0.10, 0.02])
                                    cmap = [cmap1, cmap2, cmap3, cmap4]

                                    for iax, axs in enumerate([ax12, ax13, ax14]):
                                        norm = mpl.colors.BoundaryNorm(bounds, cmap[iax+1].N)
                                        cb = mpl.colorbar.ColorbarBase(axs, cmap=cmap[iax+1], norm=norm,
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
                                    cb1.ax.set_xticks(bounds)
                                    cb1.ax.set_xticklabels([str(int(ibound)) for ibound in bounds], rotation=-45, fontsize=6)
                                    cb1.ax.tick_params(labelsize=5, width=0.25, direction='in')
                                    if srcType == 'percentchange':
                                        cb1.set_label('relative change [%]', fontsize=6, labelpad=-0.9)
                                    elif srcType == 'abschange':
                                        cb1.set_label('{} [{}]'.format(srcType, dictUnit[variable]), fontsize=6, labelpad=-0.6)
                                    cb1.outline.set_visible(False)
                                    fig.text(ax_pos.x0+0.015, ax_pos.y0+0.140, str(int(agreeThrsh[2]*1e2)), va='center', ha='center', fontsize=6)
                                    fig.text(ax_pos.x0+0.015, ax_pos.y0+0.120, str(int(agreeThrsh[1]*1e2)), va='center', ha='center', fontsize=6)
                                    fig.text(ax_pos.x0+0.015, ax_pos.y0+0.100, str(int(agreeThrsh[0]*1e2)), va='center', ha='center', fontsize=6)
                                    fig.text(ax_pos.x0+0.001, ax_pos.y0+0.130, 'agreement [%]',             va='center', ha='center', fontsize=6, rotation='vertical')

                            elif len(agreeThrsh) == 2:

                                mask1  = ma.make_mask(agreementall[ifuture, iscn]<agreeThrsh[0])
                                mask21 = ma.make_mask(agreementall[ifuture, iscn]>=agreeThrsh[0])
                                mask22 = ma.make_mask(agreementall[ifuture, iscn]<agreeThrsh[1])
                                mask2  = ma.mask_or(mask21, mask22)
                                mask3  = ma.make_mask(agreementall[ifuture, iscn]>=agreeThrsh[1])
                                signal1 = ma.masked_array(aSrc, mask=mask1)
                                signal2 = ma.masked_array(aSrc, mask=mask2)
                                signal3 = ma.masked_array(aSrc, mask=mask3)

                                ### Map   (reference: Kaye et al., Fig. 7b)
                                if srcType == 'percentchange': 
                                    bounds = [-20, -10, -5, 0, 5, 10, 20]    # Caution! bounds[0] = vmin, bounds[-1]=vmax
                                elif srcType == 'abschange':
                                    bounds = [-20, -10, -5, 0, 5, 10, 20]
                                colors1 = divide([[170,  0,  0],[230,120,  0],[240,210,  0],[127,  0,255],[ 39,161,242],[ 23, 23,171]],255.)
                                colors2 = divide([[198, 84, 84],[238,164, 84],[244,224, 84],[175, 94,255],[125,202,253],[102,109,206]],255.)
                                colors3 = divide([[210,210,210],[210,210,210],[210,210,210],[210,210,210],[210,210,210],[210,210,210]],255.)  # gray
                                cmap1 = mpl.colors.ListedColormap(colors1)
                                cmap2 = mpl.colors.ListedColormap(colors2)
                                cmap3 = mpl.colors.ListedColormap(colors3)
                                norm1 = mpl.colors.BoundaryNorm(bounds, cmap1.N)
                                norm2 = mpl.colors.BoundaryNorm(bounds, cmap2.N)
                                norm3 = mpl.colors.BoundaryNorm(bounds, cmap3.N)
                                ims2 = bm1.imshow(flipud(signal2[11:293]), norm=norm2, cmap=cmap2, interpolation='nearest')
                                ims1 = bm1.imshow(flipud(signal1[11:293]), norm=norm1, cmap=cmap1, interpolation='nearest')
                                ims3 = bm1.imshow(flipud(signal3[11:293]), norm=norm3, cmap=cmap3, interpolation='nearest')
                                bm1.drawcoastlines(linewidth=0.005, color='#808080')
                                  
                                if tag == dict_variable[variable]+'3':
                                    ### 2D colorbar for Kaye et al.-plot:
                                    ax11 = fig.add_axes([ax_pos.x0+0.028, ax_pos.y0+0.080, 0.10, 0.02])
                                    ax12 = fig.add_axes([ax_pos.x0+0.028, ax_pos.y0+0.100, 0.10, 0.02])
                                    ax13 = fig.add_axes([ax_pos.x0+0.028, ax_pos.y0+0.120, 0.10, 0.02])
                                    cmap = [cmap1, cmap2, cmap3]

                                    for iax, axs in enumerate([ax12, ax13]):
                                        norm = mpl.colors.BoundaryNorm(bounds, cmap[iax+1].N)
                                        cb = mpl.colorbar.ColorbarBase(axs, cmap=cmap[iax+1], norm=norm,
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
                                    cb1.ax.set_xticks(bounds)
                                    cb1.ax.set_xticklabels([str(int(ibound)) for ibound in bounds], rotation=-45, fontsize=6)
                                    cb1.ax.tick_params(labelsize=5, width=0.25, direction='in')
                                    if srcType == 'percentchange':
                                        cb1.set_label('relative change [%]', fontsize=6, labelpad=-0.9)
                                    elif srcType == 'abschange':
                                        cb1.set_label('{} [{}]'.format(srcType, dictUnit[variable]), fontsize=6, labelpad=-0.6)
                                    cb1.outline.set_visible(False)
                                    fig.text(ax_pos.x0+0.015, ax_pos.y0+0.120, str(int(agreeThrsh[1]*1e2)), va='center', ha='center', fontsize=6)
                                    fig.text(ax_pos.x0+0.015, ax_pos.y0+0.100, str(int(agreeThrsh[0]*1e2)), va='center', ha='center', fontsize=6)
                                    fig.text(ax_pos.x0+0.001, ax_pos.y0+0.130, 'agreement [%]',             va='center', ha='center', fontsize=6, rotation='vertical')

                            elif len(agreeThrsh) == 1:
                                mask1  = ma.make_mask(agreementall[ifuture, iscn]<agreeThrsh[0])   # mask <60%
                                mask2  = ma.make_mask(agreementall[ifuture, iscn]>=agreeThrsh[0])   # mask 60%<
                                signal1 = ma.masked_array(aSrc, mask=mask1); print('signal1: {}-{}'.format(signal1.min(), signal1.max()))
                                signal2 = ma.masked_array(aSrc, mask=mask2); print('signal2: {}-{}'.format(signal2.min(), signal2.max()))
                                ### Map   (reference: Kaye et al., Fig. 7b)
                                if srcType == 'percentchange': 
                                    #bounds = [-20, -10, -5, 0, 5, 10, 20]    # Caution! bounds[0] = vmin, bounds[-1]=vmax
                                    bounds = [-30, -15, -5, 0, 5, 15, 30]    # Caution! bounds[0] = vmin, bounds[-1]=vmax
                                elif srcType == 'abschange':
                                    bounds = [-20, -10, -5, 0, 5, 10, 20]
                                colors1 = divide([[213, 45, 37],[250,140, 87],[252,225,143],[221,241,246],[143,190,218],[ 68,116,180]],255.)
                                colors2 = divide([[210,210,210],[210,210,210],[210,210,210],[210,210,210],[210,210,210],[210,210,210]],255.)  # gray
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
                                    ###im3 = ax.contourf(lons, lats, signal3[:293], levels=bounds, norm=norm3, cmap=cmap3, transform=ccrs.PlateCarree())
                                    ###im2 = ax.contourf(lons, lats, signal2[:293], levels=bounds, norm=norm2, cmap=cmap2, transform=ccrs.PlateCarree())
                                    ###im1 = ax.contourf(lons, lats, signal1[:293], levels=bounds, norm=norm1, cmap=cmap1, transform=ccrs.PlateCarree())
                                    im2 = ax.contourf(lons, lats, signal2[:293], levels=bounds, cmap=cmap2, transform=ccrs.PlateCarree())
                                    im1 = ax.contourf(lons, lats, signal1[:293], levels=bounds, cmap=cmap1, transform=ccrs.PlateCarree())
                                    ax.add_feature(cfeature.COASTLINE, linewidth=coastal_line_width, edgecolor=coastal_line_color)

                                if tag == dict_variable[variable]+'3':
                                    #dx, dx2, dx3, width = 0.028, 0.015, 0.001, 0.1
                                    dx, dx2, dx3, width = 0.005, 0.015, 0.001, 0.113
                                    #####ax12 = fig.add_axes([ax_pos.x0+dx, ax_pos.y0+0.130, width, 0.03])
                                    ax11 = fig.add_axes([ax_pos.x0+dx, ax_pos.y0+0.100, width, 0.03])
                                    #####ax14 = fig.add_axes([ax_pos.x0+0.028, ax_pos.y0+0.190, 0.10, 0.03])  # dummy...
                                    cmap = [cmap1, cmap2]  #####, cmap3]
                                    #####for iax, axs in enumerate([ax12]): #####, ax13]):
                                    #####    norm = mpl.colors.BoundaryNorm(bounds, cmap[iax+1].N)
                                    #####    cb = mpl.colorbar.ColorbarBase(axs, cmap=cmap[iax+1], norm=norm, boundaries=bounds, spacing='uniform', orientation='horizontal')
                                    #####    cb.set_ticks(bounds)
                                    #####    cb.set_ticklabels([])
                                    #####    cb.ax.tick_params(labelsize=5, width=0.25, direction='in')
                                    #####    cb.outline.set_visible(False)
                                    norm = mpl.colors.BoundaryNorm(bounds, cmap[0].N)
                                    cb1 = mpl.colorbar.ColorbarBase(ax11, cmap=cmap[0], norm=norm, boundaries=bounds, spacing='uniform', orientation='horizontal')
                                    cb1.set_ticks(bounds)
                                    cb1.ax.set_xticklabels([str(int(ibound)) for ibound in bounds], rotation=-45, fontsize=7.5)
                                    cb1.ax.tick_params(labelsize=7.55, width=0.25, direction='in')
                                    if srcType == 'percentchange': cb1.set_label('Relative change [%]',                      fontsize=6.5, labelpad=-0.9)
                                    elif srcType == 'abschange':   cb1.set_label('{} [{}]'.format(srcType, dictUnit[index]), fontsize=6.5, labelpad=-0.6)
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

                        else:   # historical absolute or change

                            #if srcType == 'abschange':
                            #    aSrc = ma.masked_array(ma.masked_equal(mmeSrc[i],0), mask=dictKSmask[prd][ssn]['all'])
                            #    bounds = [-60, -30, 0, 30, 60]
                            #    cmap = cm.coolwarm
                            if srcType == 'absolute':
                                aSrc = mmeSrc[ifuture, iscn]
                                bounds = [0, 10, 20, 30, 40, 50]
                                cmap = cm.rainbow
                                norm1 = mpl.colors.BoundaryNorm(bounds, cmap.N)

                            ims1 = bm1.imshow(flipud(aSrc[11:293]), norm=norm1, cmap=cmap, vmin=bounds[0], vmax=bounds[-1], interpolation="nearest")
                            bm1.drawcountries(linewidth=0.2)
                            bm1.drawcoastlines(linewidth=0.1)

                            if tag == dict_variable[variable]+'3':
                                ax11 = fig.add_axes([ax_pos.x0+0.02, ax_pos.y0+0.080,  0.10, 0.03])
                                cb1 = mpl.colorbar.ColorbarBase(ax11, cmap=cmap, norm=norm1, orientation='horizontal')
                                cb1.set_ticks(bounds)
                                cb1.set_ticklabels([str(int(ibound)) for ibound in bounds])
                                cb1.ax.tick_params(labelsize=5, width=0.25, direction='in')
                                cb1.outline.set_visible(False)

                                if srcType == 'absolute':
                                    cb1.set_label('[{}]'.format(dictUnit[variable]), fontsize=5, labelpad=-0.6)

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

                        #ax.text(0., 1., '({}) {} {}'.format(tag,scn,future), va="top", ha="left", fontsize=8, transform=ax.transAxes)
                        ax.text(0.,  1., '({})'.format(tag),         va="bottom", ha="left",   fontsize=8, transform=ax.transAxes)
                        ax.text(0.5, 1., '{} {}'.format(scn,future), va="bottom", ha="center", fontsize=8, transform=ax.transAxes)

                        ## information of change
                        #if 'change' in srcType:
                        #    low_agreement_mask = ma.make_mask(agreementall[ifuture, iscn]<agreeThrsh[2])
                        #    _aSrc = ma.masked_array(aSrc, mask=low_agreement_mask)
                        #    decrease_rate = change_arearate('decrease', _aSrc)
                        #    increase_rate = change_arearate('increase', _aSrc)
                        #    df.loc['decrease', '{}-{}'.format(scn, future)] = decrease_rate
                        #    df.loc['increase', '{}-{}'.format(scn, future)] = increase_rate
                        #    if srcType == 'percentchange':
                        #        significant_level = 100
                        #        significant_increase_rate = change_arearate('significant_increase', _aSrc, significant_level)
                        #        df.loc['significant_increase', '{}-{}'.format(scn, future)] = significant_increase_rate

  
                for suffix in suffixes:
                    outFile = 'globMap.{}.Ensemble{}.{}.{}.{}'.format(variable, mmeType, srcType, prd, suffix)
                    outDir = os.path.join(fig_directory_main, variable, mmeType)
                    if not os.path.isdir(outDir): os.makedirs(outDir)
                    outPath = os.path.join(outDir, outFile)
                    if os.path.exists(outPath):
                        print("File exists, will be overwritten. ")
                        #raw_input("Print figure? Press key to continue...")
                    plt.savefig(outPath, dpi=dpi)
                    print('savefig: {}\n'.format(outPath))
                plt.close(1)

                change_excel_name = 'globMap.{}.Ensemble{}.{}.{}.xlsx'.format(variable, mmeType, srcType, prd)
                change_excel_path = os.path.join(outDir, change_excel_name)
                writer = pd.ExcelWriter(change_excel_path)
                df.to_excel(writer, 'change_area_rate')
                writer.close()
                print('ExcelWriter: {}\n'.format(change_excel_path))


            if S2NMap and s2nSrc is not None:
                #========#
                # 2. S2N #
                #========#
                
                print('Just a moment. preparing a figure object... (num=2)')
                fig = plt.figure(num=2, figsize=(8, 3.8))
                #plt.hold(True)
                print('GridSpec...')
                gs = gridspec.GridSpec(2, 2)  # (rows,cols)
                gs.update(left=0.007, right=0.995, bottom=0.01, top=0.975, hspace=0.01, wspace=0.01)

                ax1 = plt.subplot(gs[0,0])
                ax2 = plt.subplot(gs[0,1])
                ax3 = plt.subplot(gs[1,0])
                ax4 = plt.subplot(gs[1,1])
                axes = [ax1, ax2, ax3, ax4]

                istep = 0
                for ifuture, future in enumerate(futures):  # [mid, far]
                    for iscn, (scn, soc) in enumerate(SIMs):  # [rcp26, rcp85]
                        if not (ifuture == 0 and iscn == 0): istep = istep + 1
                        ax = axes[istep]
                        tag = dict_variable[variable]+['1', '2', '3', '4'][istep]

                        print('check... {}-{}'.format(s2nSrc[ifuture,iscn,11:293,:].min(), s2nSrc[ifuture,iscn,11:293,:].max()))

                        print('ax{}..'.format(istep+1))
                        plt.sca(ax)
                        ax.axis('off')

                        if   mmeType == 'median':
                            bounds = [0,0.5,1,1.5,2,2.5,3,3.5]
                            #colors1 = divide([[28,125,199], [0,0,180], [238,172,172], [228,121,121], [239,0,0], [198,0,0], [158,0,0]], 255.)
                            colors1 = divide([[153,204,255], [0,0,180], [238,172,172], [228,121,121], [239,0,0], [198,0,0], [158,0,0]], 255.)
                            labelName = 'signal to noize ratio [-]'
                        elif mmeType == 'mean':
                            bounds = [0,0.05,0.1,0.5,1,1.5,2,2.5]
                            colors1 = divide([[28,125,199], [0,0,180], [238,172,172], [228,121,121], [239,0,0], [198,0,0], [158,0,0]], 255.)
                            labelName = 'coefficient of variation [-]'
                        cmap = colors.ListedColormap(colors1)
                        norm = colors.BoundaryNorm(bounds, cmap.N)

                        ims = bm1.imshow(flipud(s2nSrc[ifuture,iscn,11:293,:]), norm=norm, cmap=cmap, vmin=bounds[0], vmax=bounds[-1], interpolation="nearest")
                        bm1.drawcountries(linewidth=0.1)
                        bm1.drawcoastlines(linewidth=0.1)

                        #ax.text(0., 1., '({}) {} {}'.format(tag,scn,future), va="top", ha="left", fontsize=8, transform=ax.transAxes)
                        ax.text(0., 1.,  '({})'.format(tag),         va="bottom", ha="left",   fontsize=8, transform=ax.transAxes)
                        ax.text(0.5, 1., '{} {}'.format(scn,future), va="bottom", ha="center", fontsize=8, transform=ax.transAxes)

                        if tag == dict_variable[variable]+'3':
                            ax_pos = ax.get_position()
                            ax11 = fig.add_axes([ax_pos.x0+0.02, ax_pos.y0+0.1,  0.10, 0.03])
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
                            cb1.set_label(labelName, fontsize=7)#, labelpad=-0.6)
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
                    outFile = 'globMap.{}.S2N.Ensemble{}.{}.{}.{}'.format(variable, mmeType, srcType, prd, suffix)
                    outDir = os.path.join(fig_directory_main, variable, mmeType)
                    if not os.path.isdir(outDir): os.makedirs(outDir)
                    outPath = os.path.join(outDir, outFile)
                    if os.path.exists(outPath): 
                        print("File exists, will be overwritten. ")
                        #raw_input("Print figure? Press key to continue...")
                    plt.savefig(outPath, dpi=dpi)
                    print('savefig: {}\n'.format(outPath))
                plt.close(2)


            if GHMMap:
                #==============#
                # 3. GHM table #
                #==============#
                print('\nJust a moment. preparing a figure object... (num=3)')
                fig = plt.figure(num=3, figsize=(len(ghms)+1, 4))
                #hold(True)
        
                print('GridSpec...')
                gs = gridspec.GridSpec(len(ghms)+1, 4)  # (rows,cols)
                gs.update(left=0.03, right=0.99, bottom=0.02, top=0.96, hspace=0.005, wspace=0.0025)
                """
                axes = [plt.subplot(gs[0,0]), plt.subplot(gs[0,1]), plt.subplot(gs[0,2]), plt.subplot(gs[0,3]),
                        plt.subplot(gs[1,0]), plt.subplot(gs[1,1]), plt.subplot(gs[1,2]), plt.subplot(gs[1,3]),
                        plt.subplot(gs[2,0]), plt.subplot(gs[2,1]), plt.subplot(gs[2,2]), plt.subplot(gs[2,3]),
                        plt.subplot(gs[3,0]), plt.subplot(gs[3,1]), plt.subplot(gs[3,2]), plt.subplot(gs[3,3]),
                        plt.subplot(gs[4,0]), plt.subplot(gs[4,1]), plt.subplot(gs[4,2]), plt.subplot(gs[4,3]),
                        plt.subplot(gs[5,0]), plt.subplot(gs[5,1]), plt.subplot(gs[5,2]), plt.subplot(gs[5,3])
                        ]
                """
                axes = []
                for ighm in range(len(ghms)):
                    axes.append(plt.subplot(gs[ighm,0]))
                    axes.append(plt.subplot(gs[ighm,1]))
                    axes.append(plt.subplot(gs[ighm,2]))
                    axes.append(plt.subplot(gs[ighm,3]))

                char_list = [chr(ichr) for ichr in range(97, 97+4*len(ghms))]
                norm1 = mpl.colors.Normalize(vmin=vmin2, vmax=vmax2)

                for g, ghm in enumerate(ghms):
                    istep = 0
                    for ifuture, future in enumerate(futures):  # [mid, far]
                        for iscn, (scn, soc) in enumerate(SIMs):  # [rcp26, rcp85]
                            if not (ifuture == 0 and iscn == 0): istep = istep + 1
                            gs = g*4+istep
                            print('{} {} {} (ax{})...'.format(ghm, future, scn, istep))

                            aSrc = ghmSRC[ifuture,iscn,g][11:293]

                            ax = axes[gs]
                            plt.sca(ax)
                            ax.axis('off')
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
                                if srcType == 'percentchange': 
                                    #bounds = [-10, -5, -1, 0, 1, 5, 10]  # Caution! bounds[0] = vmin, bounds[-1]=vmax
                                    bounds = [-20, -10, -5, 0, 5, 10, 20]  # Caution! bounds[0] = vmin, bounds[-1]=vmax
                                elif srcType == 'abschange': bounds = [-20, -10, -5, 0, 5, 10, 20]
                                ### blue >> green >> yellow >> orange >> red (5)
                                colors1 = divide([[170,  0,  0],[230,120,  0],[240,210,  0],[127,  0,255],[ 39,161,242],[ 23, 23,171]],255.)
                                colors2 = divide([[198, 84, 84],[238,164, 84],[244,224, 84],[175, 94,255],[125,202,253],[102,109,206]],255.)
                                colors3 = divide([[226,168,168],[246,209,168],[249,239,168],[215,175,255],[164,219,255],[195,195,240]],255.)
                                colors4 = divide([[210,210,210],[210,210,210],[210,210,210],[210,210,210],[210,210,210],[210,210,210]],255.)  # gray
                                cmap1 = mpl.colors.ListedColormap(colors1)
                                cmap2 = mpl.colors.ListedColormap(colors2)
                                cmap3 = mpl.colors.ListedColormap(colors3)
                                cmap4 = mpl.colors.ListedColormap(colors4)

                                ims3 = bm2.imshow(flipud(signal3), norm=norm1, cmap=cmap3, vmin=bounds[0], vmax=bounds[-1], interpolation="nearest");
                                ims2 = bm2.imshow(flipud(signal2), norm=norm1, cmap=cmap2, vmin=bounds[0], vmax=bounds[-1], interpolation="nearest");
                                ims1 = bm2.imshow(flipud(signal1), norm=norm1, cmap=cmap1, vmin=bounds[0], vmax=bounds[-1], interpolation="nearest");
                                ims4 = bm2.imshow(flipud(signal4), norm=norm1, cmap=cmap4, vmin=bounds[0], vmax=bounds[-1], interpolation="nearest");
                                #bm2.drawcountries(linewidth=0.05)
                                bm2.drawcoastlines(linewidth=0.05)

                            else:   # historical absolute or change
                                #if srcType == 'abschange':
                                #    ks_mask = dictKSmask[prd][ssn][ghm][11:293]
                                #    aSrc = ma.masked_array(ma.masked_equal(aSrc,0), mask=ks_mask)
                                #    bounds = [-60,-30,0,30,60]
                                #    cmap = cm.bwr
                                if srcType == 'absolute':
                                    bounds = [0,10,20,30,40,50,]
                                    cmap = cm.jet

                                ims1 = bm2.imshow(flipud(aSrc), norm=norm1, cmap=cmap, vmin=bounds[0], vmax=bounds[-1], interpolation='nearest')
                                #bm2.drawcountries(linewidth=0.05)
                                bm2.drawcoastlines(linewidth=0.05)

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
                        cb1.set_label('{} [{}]'.format(srcType, dictUnit[variable]), fontsize=4.5)  # labelpad=-0.6)
                    fig.text(ax_pos.x0+0.09,  ax_pos.y0-0.032, str(int(agreeThrsh[2]*1e2)), va='center', ha='center', fontsize=4.5)
                    fig.text(ax_pos.x0+0.09,  ax_pos.y0-0.052, str(int(agreeThrsh[1]*1e2)), va='center', ha='center', fontsize=4.5)
                    fig.text(ax_pos.x0+0.09,  ax_pos.y0-0.072, str(int(agreeThrsh[0]*1e2)), va='center', ha='center', fontsize=4.5)
                    fig.text(ax_pos.x0+0.075, ax_pos.y0-0.065, 'agreement [%]',             va='center', ha='center',  rotation='vertical', fontsize=4.5)

                else:
                    ax_pos = axes[-1].get_position()
                    ax11 = fig.add_axes([ax_pos.x0+0.06, ax_pos.y0-0.080,  0.15, 0.03])
                    cb1 = mpl.colorbar.ColorbarBase(ax11, cmap=cmap, norm=norm1, orientation='horizontal')
                    cb1.set_ticks(bounds)
                    cb1.set_ticklabels([str(int(ibound)) for ibound in bounds])
                    cb1.ax.tick_params(labelsize=5, width=0.25, direction='in')
                    cb1.outline.set_visible(False)
                    if srcType == 'absolute':
                        cb1.set_label('[{}]'.format(dictUnit[variable]), fontsize=5, labelpad=-0.6)
                
                for suffix in suffixes:
                    outFile = 'globMap.{}.GHM.Ensemble{}.{}.{}.{}'.format(variable, mmeType, srcType, prd, suffix)
                    outDir = os.path.join(fig_directory_main, variable, mmeType)
                    if not os.path.isdir(outDir): os.makedirs(outDir)
                    outPath = os.path.join(outDir,outFile)
                    if os.path.exists(outPath): 
                        print("File exists, will be overwritten.")
                        #raw_input("Print figure? Press key to continue...")
                    plt.savefig(outPath, dpi=dpi)
                    print('savefig: {}\n'.format(outPath))
                plt.close(3)

            if UncertaintySource and uncertainty_source_rate is not None:
                #=======================#
                # 4. Uncertainty source #
                #=======================#

                aSrc = uncertainty_source_rate[:,:,11:293,:] * 100

                bounds = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
                norm = colors.Normalize()

                print('Just a moment. preparing a figure object... (num=4)')
                fig = plt.figure(num=4, figsize=(8, 3.8))
                #plt.hold(True)
                print('GridSpec...')
                gs = gridspec.GridSpec(2, 2)  # (rows,cols)
                gs.update(left=0.007, right=0.995, bottom=0.01, top=0.975, hspace=0.01, wspace=0.01)
                ax1 = plt.subplot(gs[0,0])
                ax2 = plt.subplot(gs[0,1])
                ax3 = plt.subplot(gs[1,0])
                ax4 = plt.subplot(gs[1,1])
                axes = [ax1, ax2, ax3, ax4]

                istep = 0
                for ifuture, future in enumerate(futures):  # [mid, far]
                    for iscn, (scn, soc) in enumerate(SIMs):  # [rcp26, rcp85]
                        if (ifuture == 0 and iscn == 0): istep = istep + 1
                        ax = axes[istep]
                        tag = dict_variable[variable]+['1', '2', '3', '4'][istep]

                        print('ax{}..'.format(istep+1))
                        plt.sca(ax)
                        ax.axis('off')

                        #ax.set_title('Ratio of GCM variation to total variation', fontsize=8)
                        ims1 = bm1.imshow(flipud(aSrc[ifuture,iscn]), norm=norm, cmap=cm.RdYlBu, vmin=bounds[0], vmax=bounds[-1], interpolation="nearest")
                        bm1.drawcountries(linewidth=0.1)
                        bm1.drawcoastlines(linewidth=0.1)
                        #ax.text(0., 1., '({}) {} {}'.format(tag,scn,future), va="top", ha="left", fontsize=8, transform=ax.transAxes)
                        ax.text(0., 1., '({})'.format(tag),         va="bottom", ha="left",   fontsize=8, transform=ax.transAxes)
                        ax.text(0.5, 1., '{} {}'.format(scn,future), va="bottom", ha="center", fontsize=8, transform=ax.transAxes)

                        if tag == dict_variable[variable]+'4':
                            ax_pos = ax.get_position()
                            ax_pos = ax.get_position()
                            #ax2 = fig.add_axes([ax_pos.x0+0.04, ax_pos.y0+0.16,  0.20, 0.04])
                            #ax2 = fig.add_axes([ax_pos.x0+0.015, ax_pos.y0+0.14,  0.225, 0.04])
                            ax2 = fig.add_axes([ax_pos.x0+0.2, ax_pos.y0+0.055,  0.225, 0.02])
                            cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=cm.RdYlBu, norm=norm, orientation='horizontal')
                            cb1.set_ticks(bounds)
                            cb1.set_ticklabels([str(int(ibound)) for ibound in bounds])
                            cb1.ax.tick_params(labelsize=5, width=0.25, direction='in')
                            cb1.outline.set_visible(False)
                            cb1.set_label('[%]', fontsize=5)#, labelpad=-0.6)

                for suffix in suffixes:
                    outFile = 'globMap.{}.UNCSRC.Ensemble{}.{}.{}.{}'.format(variable, mmeType, srcType, prd, suffix)
                    outDir = os.path.join(fig_directory_main, variable, mmeType)
                    if not os.path.isdir(outDir): os.makedirs(outDir)
                    outPath = os.path.join(outDir,outFile)
                    if os.path.exists(outPath): 
                        print("File exists, will be overwritten.")
                        #raw_input("Print figure? Press key to continue...")
                    plt.savefig(outPath, dpi=dpi)
                    print('savefig: {}\n'.format(outPath))
                plt.close(4)

        endTime = datetime.datetime.now()
        diffTime = endTime - strTime
        print('end @{}'.format(endTime.strftime("%Y-%m-%d %H:%M:%S")))
        print('took {} min in total.'.format(int(diffTime.seconds/60)))

    return

if __name__=='__main__':
    main(*sys.argv)


