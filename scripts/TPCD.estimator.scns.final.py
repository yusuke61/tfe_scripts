#!/usr/local/bin/python

import os, sys, re
import itertools
import datetime
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy  as np
import pandas as pd
import matplotlib as mpl; mpl.use('Agg')
import decimal
import gc
import types
import psutil
from matplotlib import colors, ticker
from matplotlib.collections import LineCollection
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
from scipy import stats
from numpy import array, ma, fromfile, mean, median, average, where
from netCDF4 import Dataset
from multiprocessing import Pool
from weighted_kde import weighted_gaussian_kde
from copy import copy
from memory_profiler import profile
from utiltools import get_region_info
mpl.rcParams['axes.linewidth'] = 0.5
hostname = os.uname()[1]


# OPTIONS ==================================================================================
LOUD = True

TEST = False; DUMMY = False  # default
#TEST = True; DUMMY = True  # for test

REUSE = False
#REUSE = True


FigsMininum_for_NCommn = True

how_to_gen_member_pdf = 'newly create'
#how_to_gen_member_pdf = 'reuse'

how_to_gen_ensemble_pdf = 'median of member pdfs'
#how_to_gen_ensemble_pdf = 'reuse'
#how_to_gen_ensemble_pdf = 'create from all members'

which_member_timeseries = 'newly create'  # default
#which_member_timeseries = 'reuse'

how_to_plot_ensemble_timeseries = 'median of member timeseries'  # default
#how_to_plot_ensemble_timeseries = 'reuse'
#how_to_plot_ensemble_timeseries = 'extract from a pdf'

how_to_estimate_ensemble_TPCD = 'median of member TPCDs'
#how_to_estimate_ensemble_TPCD = 'directory from a representative timeseries'

READ_VMAX_OVERALL = False  # default
#READ_VMAX_OVERALL = True

which_hist_percentile = 'new'  # default
#which_hist_percentile = 'reuse'

if REUSE == True:
    how_to_gen_member_pdf = 'reuse'
    which_member_timeseries = 'reuse'
    which_hist_percentile = 'reuse'
    how_to_plot_ensemble_timeseries = 'reuse'
    how_to_gen_ensemble_pdf = 'reuse'
    which_date_pdf = '20210326'
    which_date_excel = '20210326'
    which_tChunk_reuse = int(sys.argv[8])
    which_threshold_reuse = sys.argv[7]

TPCD2EXCEL = True

VRANGE_DETAIL_TO_EXCEL = True
#VRANGE_DETAIL_TO_EXCEL = False

TPCD2BIN = False

if how_to_gen_member_pdf == 'newly create':
    PDFs2NC_OVERALL = True
    PDFs2NC_MEMBER = True
else:
    PDFs2NC_OVERALL = False
    PDFs2NC_MEMBER = False

# --- two options for the memory issue...
# Manual Gabage collection
#GC = False
GC = True

MEMORY_USAGE_NOISY = False
#MEMORY_USAGE_NOISY = True

experiment = sys.argv[1]

# select stats value(s) to pick up from a kernel density function -------------------------
smplTYPs = [sys.argv[2]]
#smplTYPs = ['Average', 'Extreme']
#smplTYPs = ['Average']  # Average from KED
#smplTYPs = ['Median']
#smplTYPs = ['Extreme']
#smplTYPs = ['Mean', 'Average', 'Median', 'Mode', 'Extreme']
#if TEST: smplTYPs = ['Average']

# gcms -------------------------------------------------------------------------------------
if sys.argv[3] == 'all': gcms = ['hadgem2-es',  'ipsl-cm5a-lr', 'gfdl-esm2m', 'miroc5']  # full members
else: gcms = [sys.argv[3]]
if TEST: gcms = ['ipsl-cm5a-lr']

# ghms -------------------------------------------------------------------------------------
if sys.argv[4] == 'all': ghms = ['cwatm', 'h08', 'lpjml', 'matsiro', 'watergap2']  # full members
else: ghms = [sys.argv[4]]
if TEST: ghms = ['matsiro']

# select scenarios ------------------------------------------------------------------------ 
if sys.argv[5] == 'all': scns = ['rcp85', 'rcp26']  # basic
else: scns = [sys.argv[5]]
if TEST: scns = ['rcp85']

# socs & co2s -----------------------------------------------------------------------------
socs = ['2005soc']
co2s = ['co2']

# --- update experimental setting 
if experiment == 'basic':
    pass
elif experiment == 'picontrol':
    scns = ['picontrol', 'rcp26', 'rcp85']
elif experiment == 'co2exp':  # CO2 experiment  (LPJmL & MATSIRO)
    ghms = ['lpjml', 'matsiro']
    co2s = ['co2', '2005co2']
elif experiment == 'rcp26soc':
    ghms = ['cwatm', 'h08', 'lpjml', 'matsiro']
    scns = ['rcp26']
    socs = ['2005soc', 'rcp26soc']
else:
    raise ValueError

# select a spatial scale you would like to estimate TPCD for ------------------------------
region_type = sys.argv[6]
#region_type = 'AR6_regions'
#region_type = 'SREX_regions'
#region_type = 'HydroBASINS_lev1'
#region_type = 'hydroregions'
#region_type = 'BasinCountryUnit'
#region_type = 'Basin'
#region_type = 'Nation'

# select if you'd like to process absolute value or change
value_type = sys.argv[7]
#value_type = 'abs'
#value_type = 'anm'

# select how to deal with yealy sampling to create a PDF  ----------------------------------
kdfSampling = 'window'; window = int(sys.argv[8])  # ?yr-sampling in total
#kdfSampling = 'window'; window = 5  # 5yr-sampling in total  default
#kdfSampling = 'window'; window = 7  # 7yr-sampling in total
#kdfSampling = 'yearly'

# which period do you want to refer in the TPCD analysis?  ----------------------------------
which_reference_period_for_TPCD = sys.argv[9]
#which_reference_period_for_TPCD = '2b_historical_full'  # 1861-2005
#which_reference_period_for_TPCD = '2b_bc_period'        # 1979-2013

# which threshold type do you want to apply?  -----------------------------------------------
threshold_type = sys.argv[10]
#threshold_type = 'max'
#threshold_type = '99percentile'
#threshold_type = '98percentile'
#threshold_type = '1sigma'
#threshold_type = '2sigma'

# select consecutive chunk size for the TPCD analysis (Exceedance for consecutively ??years
tChunkType = [int(sys.argv[11]) if sys.argv[11] != 'full' else sys.argv[11]]
#tChunkType = [1, 2, 5, 10, 20, 30, 'full']
#tChunkType = ['full']

# select which season you would like to analyze about --------------------------------------
#seasons = ['ALL','DJF','MAM','JJA','SON']
#seasons = ['ALL']
seasons = [sys.argv[12]]

#process_info = sys.argv[12]


print('='*10 + ' settings ' + '='*30)
for i, arg in enumerate(['experiment', 'gcms', 'ghms', 'scns', 'socs', 'co2s', 
                         'region_type', 'value_type', 'window', 
                         'which_reference_period_for_TPCD', 
                         'threshold_type', 'tChunkType']):
    print('arg{:>2} {}: {}'.format(i, arg, eval(arg)))
print('='*30)


# select figures to output ----------------------------------------------------------------
Figure1  = False  # PDF timeseries on basemap
Figure2  = False  # PDF timeseries in the list style                                  (SupFig7)
Figure3  = False  # plot of timeseries of aSrc (nYear)    ensemble
Figure4  = False  # single Kernel timeseries              ensemble                                  is this better than Figure2?
Figure5  = False  # plot of timeseries of aSrc (nYear)    member                      (Fig1)&(ExtDataFig3)
Figure6  = False  # single Kernel timeseries              member                               is this better than Figure2?
Figure7  = False  # TPCD global map (TPCD)                                            (Fig2)
Figure8  = False  # TPCD global map (STD in TPCD)
Figure9  = False  # TPCD global map (Earliest&latest TPCD among members)
Figure10 = False  # The member heat map of TPCD
Figure11 = False  # plot of timeseries of ensemble aSrcs (Global&Regional) on a basemap     x(Fig2) no longer used
Figure12 = False  # plot of timeseries of ensemble aSrcs for each region
Figure13 = False  # plot of timeseries of anomaly in ensemble aSrcs (Global&Regional) on a basemap
Figure14 = False  # population plot
Figure15 = False  # colorfulplot_shiftingPDF
 
if FigsMininum_for_NCommn:
    #Figure2  = True  # PDF timeseries in the list style
    #Figure3  = True  # plot of timeseries of aSrc (nYear)    ensemble 
    #Figure4  = True  # single Kernel timeseries              ensemble 
    #Figure5  = True  # plot of timeseries of aSrc (nYear)    member
    #Figure6  = True  # single Kernel timeseries              member
    Figure7  = True  # TPCD global map (TPCD)
    #Figure10 = True  # The member heat map of TPCD
    #Figure12 = True  # plot of timeseries of ensemble aSrcs for each region
    #Figure14 = True  # population plot

elif region_type in ['AR6_regions', 'SREX_regions', 'HydroBASINS_lev1', 'hydroregions'] and process_info == 'Fig2Fig15':
    Figure2  = True
    Figure15 = True
    #PDFs2NC_OVERALL = False
    #PDFs2NC_MEMBER  = False
    PDFs2NC_OVERALL = True
    PDFs2NC_MEMBER  = True
elif region_type in ['AR6_regions', 'SREX_regions', 'HydroBASINS_lev1', 'hydroregions']:
    #Figure7  = True
    #Figure10 = True
    #Figure11 = True
    #Figure12 = True
    #Figure14 = True
    pass
elif region_type == 'BasinCountryUnit':
    Figure7  = True
    Figure14 = True
    
if region_type == 'BasinCountryUnit':
    Figure10 = False  # The member heat map of TPCD
    Figure12 = False  # plot of timeseries of ensemble aSrcs for each region

if TEST:
    Figure1  = False # PDF timeseries on basemap
    Figure2  = False # PDF timeseries in the list style                                  (SupFig7)
    Figure3  = True  # plot of timeseries of aSrc (nYear)    ensemble
    Figure4  = False # single Kernel timeseries              ensemble                                  is this better than Figure2?
    Figure5  = True  # plot of timeseries of aSrc (nYear)    member                      (Fig1)&(ExtDataFig3)
    Figure6  = False # single Kernel timeseries              member                               is this better than Figure2?
    Figure7  = True  # TPCD global map (TPCD)                                            (Fig2)
    Figure8  = True  # TPCD global map (STD in TPCD)
    Figure9  = True  # TPCD global map (Earliest&latest TPCD among members)
    Figure10 = True  # The member heat map of TPCD
    Figure11 = False # plot of timeseries of ensemble aSrcs (Global&Regional) on a basemap     x(Fig2) no longer used
    Figure12 = True  # plot of timeseries of ensemble aSrcs for each region
    Figure13 = False # plot of timeseries of anomaly in ensemble aSrcs (Global&Regional) on a basemap
    Figure14 = True  # population plot
    Figure15 = False # colorfulplot_shiftingPDF


# select how you express uncertainty in a timeseries plot ---------------------------------
rangeType = 'IQR'
#rangeType = 'minmax'  # if there is no outlier
#rangeType = '2sigma'  # if it looks nomal distribution type...

# select how to sample from ensemble member for each grid. default is 'all' ---------------
#mmeTYPs = ['median','mean','all']
mmeTYPs = ['all']
#mmeTYPs = ['median']
#mmeTYPs = ['mean']

# select threshold type -------------------------------------------------------------------
#autQvalue = True    # automatically specify, depending on soc type: 'Qvale_hist_%s'%soc
autQvalue = False  # set qvaltype, below
#qvalTYPE = 'Qvale_hist_nosoc'
#qvalTYPE = 'Qvale_hist_pressoc'

# drought indices --------------------------------------------------------------------------
IDXs = ['nDayTot']
#IDXs = ['nDayTot', 'dfcTot', 'nEvent', 'nOnset', 'avrDSL', 'maxDSL']

# drought analysis parameters ---------------------------------------------------------------
# --- Qvalue
#Qs = [90, 85, 80, 75, 70, 65]
Qs = [80]
# --- window size for Qvalue 
#wins = [15, 10, 7, 5]
wins = [15]
# --- minimum drought days
#LENs = [180,90,60,30,14,7] 
#LENs = [30,14,60,90,180,7] 
LENs = [30] 
# --- Pooled duration
TAUs = [4]

# --- Parameter for Extreme rate
percentile_hist = 90  # calc XX percentile of historical period


# ------------------------------------------------- Temporal settings -------------------
# reference period of drought threshold
ref_syear, ref_eyear = 1861, 2005      # isimip2b full historical period
if which_reference_period_for_TPCD == '2b_bc_period':
    ref_syear, ref_eyear = 1979, 2013  # ISIMIP2b bias-correction reference period

# (historical) reference period for TPCD analysis
if which_reference_period_for_TPCD == '2b_historical_full':
    syear_hist_org, eyear_hist_org = ref_syear, ref_eyear  # default 1861-2005
elif which_reference_period_for_TPCD == '2b_long_length':
    syear_hist_org, eyear_hist_org = 1901, 2005
elif which_reference_period_for_TPCD == '2b_mid_length':
    syear_hist_org, eyear_hist_org = 1950, 2005
elif which_reference_period_for_TPCD == '2b_bc_period':  # Caution!! values for rcp85 are referred during 2006-2013 even for rcp26
    syear_hist_org, eyear_hist_org = 1979, 2013

# full period of drought data to analyze
syear0, eyear0 = ref_syear, 2099
years0 = range(syear0, eyear0+1)  # 1861-2099   full

if kdfSampling == 'window':  # moving average at the end of the window
    syear = syear0 + (window-1)         # 1865 = 1861 + (5-1)   # start year of the period for which yearly data is available. Start of the full period available.
    eyear = eyear0                      # 2099, as it is.     # end   year of the period for which yearly data is available. End of the full period available.
    syear_hist_org = syear_hist_org + (window-1)
    eyear_hist_org = eyear_hist_org
else:
    syear = syear0
    eyear = eyear0
years = range(syear, eyear+1)     # 1865-2099   with window
nyear = len(years)                # the number of total years analyzed  (235 years)

# the historical period of the TPCD analysis
syear_analysis_historical = max(syear_hist_org, syear)   # 1865 with the window sampling
eyear_analysis_historical = eyear_hist_org               # 2005 or 2013
years_analysis_historical = range(syear_analysis_historical, eyear_analysis_historical+1)  # 1865-2005
years_analysis_future = range(eyear_analysis_historical+1, eyear+1)  # 141 years
years_analysis = range(syear_analysis_historical, eyear+1)  # 1865-2099

# --- (just for TPCD maps)
# define ticklabels for TPCD analysis
if which_reference_period_for_TPCD == '2b_bc_period': syear_tpcdmap = 2013
else:                                                 syear_tpcdmap = eyear_hist_org + 1
eyear_tpcdmap = 2099
#ticklabels_year = range(2020, eyear_tpcdmap, 20)
ticklabels_year = range(2020, 2100, 20)

# the historical period for the time series plot
syear_historical_for_plot = syear  # 1865, 1983
eyear_historical_for_plot = 2005 # 2005
#syear_historical_for_plot = syear_analysis_historical  # 1865
#eyear_historical_for_plot = eyear_analysis_historical  # 2005
years_historical_for_plot = range(syear_historical_for_plot, eyear_historical_for_plot+1)  # ????-2005
years_future_for_plot = range(eyear_historical_for_plot, eyear+1)      # 2005-2099. To satart from the same point, it starts from 2005 here.

print(' ------ Drought detection ref ----------------------')
print('                              {}-{} ({})'.format(ref_syear, ref_eyear, len(range(ref_syear,ref_eyear+1))))
print(' ------ Original -----------------------------------')
print('  the full period:            {}-{} ({})'.format(syear0,                      eyear0,                    len(years0)))
print('  the TPCD historical period: {}-{} ({})'.format(syear_hist_org,              eyear_hist_org,            len(range(syear_hist_org, eyear_hist_org+1))))
print('  the TPCD future period:     {}-{} ({})'.format(eyear_hist_org+1,            eyear0,                    len(range(eyear_hist_org+1, eyear0+1))))
print(' ------ In the analysis with window_size={}years ---'.format(window))
print('  the full period:            {}-{} ({})'.format(syear,                       eyear,                     len(years)))
print('  the TPCD historical period: {}-{} ({})'.format(syear_analysis_historical,   eyear_analysis_historical, len(years_analysis_historical)))
print('  the TPCD future period:     {}-{} ({})'.format(eyear_analysis_historical+1, eyear,                     len(years_analysis_future)))
print(' ------ For a time series plot ---------------------')
print('  the historical period: {}-{} ({})'.format(syear_historical_for_plot,   eyear_historical_for_plot, len(years_historical_for_plot)))
print('  the future period:     {}-{} ({})'.format(eyear_historical_for_plot+1, eyear,                     len(years_future_for_plot)))
print(' ---------------------------------------------------')

# start and end year (index) of the (historical) reference period for the TPCD analysis.
si = years.index(syear_analysis_historical)       # index of the start year of the reference period in the TPCD analysis
ei = years.index(eyear_analysis_historical)       # the same but for the end of the reference period

# start and end year (index) of the historical period for the TPCD analysis.
index_hist_s_for_plot = years.index(syear_historical_for_plot)
index_hist_e_for_plot = years.index(eyear_historical_for_plot)

# for sample PDFs in colorfulplot_shiftingPDF 
sample_years = [2005, 2006, 2013, 2030, 2050, 2060, 2080]
if TEST: sample_years = [2030]

syear_pop = 2006
#hogehoge if kdfSampling == 'window': syear_pop = syear_pop + 2*window
if kdfSampling == 'window': syear_pop = syear_pop + (window-1)
eyear_pop = eyear
years_pop = range(syear_pop, eyear_pop+1)

days = range(366)  # 0,1,2,...,364,365

# --------------------------------------------------------------------------------------------
if region_type in ['AR6_regions', 'SREX_regions', 'HydroBASINS_lev1', 'hydroregions']:
    GLOBAL, DESERT = True, True   # default True
elif region_type == 'BasinCountryUnit': 
    GLOBAL, DESERT = False, False
else: raise ValueError('What do you want for GLOBAL and DESERT option??'); sys.exit()
if TEST: GLOBAL = False

TPCDmap_cmap = 'hot'
#TPCDmap_cmap = 'jet_r'
#TPCDmap_cmap = 'PuRd_r'
#TPCDmap_cmap = 'RdYlBu'

#pdf_cmap = 'gnuplot'  # This is!!
#pdf_cmap = 'stern_r'  
#pdf_cmap = 'cubehelix_r'  # NG
#pdf_cmap = 'rainbow'
pdf_cmap = 'terrain_r'  # in the first version
#pdf_cmap = 'gist_stern_r'
#pdf_cmap = 'viridis'
#pdf_cmap = 'plasma'
#pdf_cmap = 'YlOrRd'
#pdf_cmap = 'PuRd'
#pdf_cmap = 'pink_r'
#pdf_cmap = 'ocean_r'
#pdf_cmap = 'cubehelix_r'

#kernel_randomvariable_unit = 'day/year'
kernel_randomvariable_unit = '%'

# Basemap setting ------------------------------------------------------------------------
#-- Figure1: Background of the Kernel time series scattered over the basemap
fig1_proj = 'cyl'
#fig1_proj = 'robin'
fig1_backcolor = 'fillcontinents'
#fig1_backcolor = 'shadedrelief'
#fig1_backcolor = 'drawlsmask'

#-- Figure1: The region map for the Kernel time series list
fig2_proj = 'robin'
#fig2_proj = 'cyl'
fig2_backcolor = 'fillcontinents'
#fig2_backcolor = 'shadedrelief'
#fig2_backcolor = 'drawlsmask'

#-- Figure11, 11: The backgound for the drought stats time series
fig11_proj = 'cyl'
#fig11_proj = 'robin'
fig11_backcolor = 'fillcontinents'
#fig11_backcolor = 'shadedrelief'
#fig11_backcolor = 'drawlsmask'

#suffixes = ['.png']
suffixes = ['.png', '.eps', '.pdf']
dpi_figure  =  80
dpi_savefig = 300


# ----------------------------------------------------------------------------------------
if 'scs' in hostname or hostname == 'sun220s.nies.go.jp':
    if 'scs' in hostname:                  data_directory = '/data/rg001/sgec0017/data'
    elif hostname == 'sun220s.nies.go.jp': data_directory = '/sraid02/satoh/data'
    drougnt_src_directory = os.path.join(data_directory, 'isimip2b.drought')  # main droughtStats directory
    fig_directory_main    = os.path.join(data_directory, 'figure_box', 'drought_tpcd_isimip2b', 'TPCD.estimator.scns')
    population_directory  = os.path.join(data_directory, 'population')
    mapmask_directory     = os.path.join(data_directory, 'mapmask')
    landseamask_path      = os.path.join(mapmask_directory, 'ISIMIP2b_landseamask_generic.nc4')
    gaul_path             = os.path.join(mapmask_directory, 'GAUL', 'flt', 'gaul2014_05deg.flt')
    grid_index_path       = os.path.join(mapmask_directory, 'grid_index.nc4')
    gridarea_path         = os.path.join(mapmask_directory, 'grd_ara.hlf')
    drygrid_mask_path     = os.path.join(data_directory, 'figure_box', 'drought_tpcd_isimip2b', 
                                         'find_grids_always_zero_discharge',
                                         'drygridmask___combined_historical_histsoc_co2.nc4')
    if region_type == 'AR6_regions':
        region_shp_path = os.path.join(mapmask_directory, 'IPCC-AR6-WGI-reference-regions-v4_shapefile', 
                                       'shapefile_edited', 'land', 'reference-regions-AR6'
                                       )
    elif region_type == 'SREX_regions':
        region_shp_path = os.path.join(mapmask_directory, 'SREX.referenceRegions.editted', 
                                       'referenceRegions'
                                       )
    elif region_type == 'HydroBASINS_lev1':
        region_shp_path = os.path.join(mapmask_directory, 'HydroSHEDS', 'HydroBASINS', 
                                       'withLakes', 'Global_merge', 
                                       'hybas_lake____lev02_v1c_merge'
                                       )
    elif region_type == 'hydroregions':
        region_shp_path = os.path.join(mapmask_directory, 'Hydrobelt', 
                                       'hydroregions_shp', 
                                       'meybeck_et_al_2013_hydroregions'
                                       )
else:
    print('Ah... hostname is {}'.format(hostname))
    sys.exit()

# other global variables -----------------------------------------------------------------
landsea = Dataset(landseamask_path)['LSM'][:][0].filled(0)
greenland_mask = ma.make_mask(fromfile(gaul_path,'float32').reshape(360,720) == 98)
desert_mask = ma.make_mask(Dataset(grid_index_path).variables['grid_index'][:] == 11)
landsea_mask_original = ma.make_mask(Dataset(landseamask_path).variables['LSM'][:][0] != 1.0)
drygrid_mask = ma.make_mask(Dataset(drygrid_mask_path)['drygrid_mask'][:] == 1)
landsea_or_greenland_mask = ma.mask_or(landsea_mask_original, greenland_mask)    # landsea + greenland
landsea_or_greenland_mask = ma.mask_or(landsea_or_greenland_mask, drygrid_mask)  # landsea + greenland + drygrid
landsea_mask = copy(landsea_or_greenland_mask)                                   # landsea + greenland + drygrid
if DESERT: landsea_mask = ma.mask_or(landsea_mask, desert_mask)                  # landsea + greenland + drygrid + desert

gridArea = fromfile(gridarea_path,'float32').byteswap().reshape(360,720)

dictIndex = {'nDayTot': ['Drought Days', 'day/yr'], }

if region_type == 'SREX_regions':
              #  0          1              2                3         
    #regionID  Name    KernelTS(x,y)    for robin     Timeseries(x,y)
    dictIPCCregion = {                                             
           0: ['GLB', [0.035,0.130],  [0.035,0.330],  [0.030,0.200]],
           1: ['ALA', [0.050,0.670],  [0.200,0.870],  [0.060,0.785]],
           2: ['CGI', [0.285,0.680],  [0.370,0.870],  [0.220,0.785]],
           3: ['WNA', [0.115,0.650],  [0.200,0.740],  [0.025,0.640]],
           4: ['CNA', [0.175,0.650],  [0.246,0.705],  [0.150,0.640]],
           5: ['ENA', [0.230,0.630],  [0.295,0.700],  [0.275,0.640]],
           6: ['CAM', [0.155,0.330],  [0.240,0.550],  [0.115,0.495]],
           7: ['AMZ', [0.270,0.310],  [0.312,0.449],  [0.245,0.405]],
           8: ['NEB', [0.335,0.400],  [0.370,0.430],  [0.370,0.330]],
           9: ['WSA', [0.215,0.155],  [0.270,0.350],  [0.245,0.255]],
          10: ['SSA', [0.330,0.100],  [0.328,0.210],  [0.370,0.185]],
          11: ['NEU', [0.405,0.670],  [0.505,0.860],  [0.400,0.785]],
          12: ['CEU', [0.515,0.640],  [0.534,0.783],  [0.540,0.735]],
          13: ['MED', [0.460,0.585],  [0.510,0.695],  [0.415,0.640]],
          14: ['SAH', [0.390,0.360],  [0.500,0.605],  [0.402,0.490]],
          15: ['WAF', [0.445,0.270],  [0.466,0.486],  [0.495,0.330]],
          16: ['EAF', [0.555,0.270],  [0.570,0.486],  [0.618,0.330]],
          17: ['SAF', [0.500,0.208],  [0.515,0.320],  [0.568,0.185]],
          18: ['NAS', [0.805,0.670],  [0.720,0.825],  [0.820,0.785]],
          19: ['WAS', [0.570,0.585],  [0.596,0.680],  [0.550,0.538]],
          20: ['CAS', [0.625,0.660],  [0.650,0.707],  [0.680,0.640]],
          21: ['TIB', [0.680,0.660],  [0.700,0.707],  [0.675,0.785]],
          22: ['EAS', [0.750,0.590],  [0.804,0.665],  [0.805,0.640]],
          23: ['SAS', [0.640,0.320],  [0.690,0.570],  [0.680,0.495]],
          24: ['SEA', [0.730,0.290],  [0.826,0.513],  [0.805,0.495]],
          25: ['NAU', [0.790,0.230],  [0.840,0.355],  [0.805,0.350]],
          26: ['SAU', [0.850,0.150],  [0.840,0.240],  [0.820,0.205]],
          }
    n_regions = 26
elif region_type == 'AR6_regions':  # TODO
    dictIPCCregion = {
           0: [ 'GLB', [None,None], [None,None], [None,None]],  # original ID=0 is GIC in the IPCC regional definition but is replaced with GLB in the process.
           1: [ 'NWN', [None,None], [None,None], [None,None]], 
           2: [ 'NEN', [None,None], [None,None], [None,None]], 
           3: [ 'WNA', [None,None], [None,None], [None,None]], 
           4: [ 'CNA', [None,None], [None,None], [None,None]], 
           5: [ 'ENA', [None,None], [None,None], [None,None]], 
           6: [ 'NCA', [None,None], [None,None], [None,None]], 
           7: [ 'SCA', [None,None], [None,None], [None,None]], 
           8: [ 'CAR', [None,None], [None,None], [None,None]], 
           9: [ 'NWS', [None,None], [None,None], [None,None]], 
          10: [ 'NSA', [None,None], [None,None], [None,None]], 
          11: [ 'NES', [None,None], [None,None], [None,None]], 
          12: [ 'SAM', [None,None], [None,None], [None,None]], 
          13: [ 'SWS', [None,None], [None,None], [None,None]], 
          14: [ 'SES', [None,None], [None,None], [None,None]], 
          15: [ 'SSA', [None,None], [None,None], [None,None]], 
          16: [ 'NEU', [None,None], [None,None], [None,None]], 
          17: [ 'WCE', [None,None], [None,None], [None,None]], 
          18: [ 'EEU', [None,None], [None,None], [None,None]], 
          19: [ 'MED', [None,None], [None,None], [None,None]], 
          20: [ 'SAH', [None,None], [None,None], [None,None]], 
          21: [ 'WAF', [None,None], [None,None], [None,None]], 
          22: [ 'CAF', [None,None], [None,None], [None,None]], 
          23: ['NEAF', [None,None], [None,None], [None,None]], 
          24: ['SEAF', [None,None], [None,None], [None,None]], 
          25: ['WSAF', [None,None], [None,None], [None,None]], 
          26: ['ESAF', [None,None], [None,None], [None,None]], 
          27: [ 'MDG', [None,None], [None,None], [None,None]], 
          28: [ 'RAR', [None,None], [None,None], [None,None]], 
          29: [ 'WSB', [None,None], [None,None], [None,None]], 
          30: [ 'ESB', [None,None], [None,None], [None,None]], 
          31: [ 'RFE', [None,None], [None,None], [None,None]], 
          32: [ 'WCA', [None,None], [None,None], [None,None]], 
          33: [ 'ECA', [None,None], [None,None], [None,None]], 
          34: [ 'TIB', [None,None], [None,None], [None,None]], 
          35: [ 'EAS', [None,None], [None,None], [None,None]], 
          36: [ 'ARP', [None,None], [None,None], [None,None]], 
          37: [ 'SAS', [None,None], [None,None], [None,None]], 
          38: [ 'SEA', [None,None], [None,None], [None,None]], 
          39: [ 'NAU', [None,None], [None,None], [None,None]], 
          40: [ 'CAU', [None,None], [None,None], [None,None]], 
          41: [ 'EAU', [None,None], [None,None], [None,None]], 
          42: [ 'SAU', [None,None], [None,None], [None,None]], 
          43: [  'NZ', [None,None], [None,None], [None,None]], 
          }
    n_regions = 43
elif region_type == 'hydroregions': # TODO
    dictIPCCregion = {
            0: [         'GLB', [None,None], [None,None], [None,None]],  # original ID=0 is GIC in the IPCC regional definition but is replaced with GLB in the process.
           11: [     'BOR_Nam', [None,None], [None,None], [None,None]], 
           12: [     'NML_Nam', [None,None], [None,None], [None,None]], 
           13: [     'NDR_Nam', [None,None], [None,None], [None,None]], 
           14: [     'NST_Nam', [None,None], [None,None], [None,None]], 
           25: [     'EQT_Sam', [None,None], [None,None], [None,None]], 
           26: [     'SST_Sam', [None,None], [None,None], [None,None]], 
           27: [     'SDR_Sam', [None,None], [None,None], [None,None]], 
           28: [     'SML_Sam', [None,None], [None,None], [None,None]], 
           31: [     'BOR_Eur', [None,None], [None,None], [None,None]], 
           32: [     'NML_Eur', [None,None], [None,None], [None,None]], 
           43: [     'NDR_Afr', [None,None], [None,None], [None,None]], 
           44: [     'NST_Afr', [None,None], [None,None], [None,None]], 
           45: [     'EQT_Afr', [None,None], [None,None], [None,None]], 
           46: [     'SST_Afr', [None,None], [None,None], [None,None]], 
           47: [     'SDR_Afr', [None,None], [None,None], [None,None]], 
           48: [     'SML_Afr', [None,None], [None,None], [None,None]], 
          511: ['BOR_Asi(WSb)', [None,None], [None,None], [None,None]], 
          512: ['BOR_Asi(ESb)', [None,None], [None,None], [None,None]], 
           52: [     'NML_Asi', [None,None], [None,None], [None,None]], 
          531: ['NDR_Asi(MdE)', [None,None], [None,None], [None,None]], 
          532: ['NDR_Asi(CAs)', [None,None], [None,None], [None,None]], 
           54: [     'NST_Asi', [None,None], [None,None], [None,None]], 
           55: [     'EQT_Asi', [None,None], [None,None], [None,None]], 
           66: [     'SST_Aus', [None,None], [None,None], [None,None]], 
           67: [     'SDR_Aus', [None,None], [None,None], [None,None]], 
           68: [     'SML_Aus', [None,None], [None,None], [None,None]], 
          }
elif region_type == 'HydroBASINS_lev1':
    df_region = pd.read_excel(region_shp_path+'.xlsx', header=0, usecols=[0,1])
    df_region = df_region.drop_duplicates()
    dictIPCCregion = {}
    n_regions = df_region.shape[0]
    for iregion in range(n_regions):
        region_ID, long_name = df_region.iloc[iregion]
        dictIPCCregion[region_ID] = [long_name]


dictSeason = {
    'ALL': 365,
    'DJF':  90,
    'MAM':  92,
    'JJA':  92,
    'SON':  91,
    'DRY':  91,
    'WET':  91,
    }

dictSSP = {
    'rcp85': 'ssp5',
    'rcp26': 'ssp1',
    'picontrol': 'ssp1',
    }

dictSim = {
    'rcp26': '#0080ff',  # right blue      '#0000c1',
    'rcp85': '#ff6600',  #'#ff5500',  # orange  '#dc0d0d',
    #'rcp85': '#ffc400',  # dark yellow
    'picontrol': '#808080'  #'#b0b0b0'  # 
    }

dictMember = {  # this is for the check process at the member level.
    'hadgem2-es':  '#ff0000',
    'ipsl-cm5a-lr':'#800080',
    'gfdl-esm2m':  '#008000',
    'miroc5':      '#0000ff',
    # -------------------
    'cwatm':       '#ff0000',
    'h08':         '#800080',
    'lpjml':       '#008000',
    'matsiro':     '#0000ff',
    'watergap2':   '#0099ff',
    }

pT = re.compile('T+')
pF = re.compile('F+')

if region_type in ['AR6_regions', 'SREX_regions', 'HydroBASINS_lev1', 'hydroregions']:
    #multiprocessing = True; nMP = 12  # this will fail for AR6 regions
    multiprocessing = True
    if not experiment == 'picontrol': nMP = 10
    else:                             nMP = 6 
    if how_to_gen_member_pdf == 'reuse': nMP = 4
elif region_type == 'BasinCountryUnit':
    multiprocessing = True; nMP = 4
#multiprocessing = False

today = datetime.date.today().strftime('%Y%m%d')
if TEST: today = today + '_test'



# ---------------------------------------------------------------------------------------------------------------
def check_nan(src, check_id):
    try:
        if True in np.isnan(src):
            print('\n!!!!!!!!!!!!!!!!!!!!!!!!!')
            print('Caution #{} !!! There seems to be NaN in your input data !!!!'.format(check_id))
            print('Tentatively, NaN is replaced with missing_value and masked')
            print('!!!!!!!!!!!!!!!!!!!!!!!!!\n')
            src[where(np.isnan(src))] = 1e+20
            src = np.ma.masked_equal(src, 1e+20)
    except:
        print('skip check_nan.')
    return src

# ---------------------------------------------------------------------------------------------------------------
def load_nc(co2, soc, scn, gcm, ghm, dict_parameters):
    Q        = dict_parameters['Q']
    win      = dict_parameters['win']
    Len      = dict_parameters['Len']
    tau      = dict_parameters['tau']
    season   = dict_parameters['season']
    index    = dict_parameters['index']
    qvalType = dict_parameters['qvalType']
    # ---
    input_file_name = 'Q{:02}win{:02}_Len{:03}tau{}_{}_{}.nc4'.format(Q, win, Len, tau, season, index)
    input_directory = os.path.join(
        drougnt_src_directory,
        ghm, gcm, '{}.{}_{}'.format(qvalType, ref_syear, ref_eyear),
        'droughtStats.{}_{}_{}.{}_{}'.format(scn, soc, co2, syear0, eyear0)
        )
    srcpath = os.path.join(input_directory, input_file_name)
    if not os.path.isfile(srcpath): 
        if TEST and DUMMY:
            return np.random.rand(len(range(syear0, eyear0+1)), 360, 720)
        else:
            print('Error!! {} is not exist... Check!!'.format(srcpath)); sys.exit()  # TODO
    else:
        src = Dataset(srcpath).variables[index][:]
        print('Loaded : {} {}'.format(srcpath, src.shape))
        src = check_nan(src, 1)
        return src


# -----------------------------------------------------------------------------------------------------------
def absolute_to_anomaly(aSrc):  # aSrc: (nSCN, nGCM, nGHM, nYear, nY, nX)
    aSrc = np.transpose(aSrc, axes=(3,0,1,2,4,5))  # (nYear, nSCN, nGCM, nGHM, nY, nX)
    anomaly = aSrc - aSrc[si:ei].mean(axis=0)      # (nYear, nSCN, nGCM, nGHM, nY, nX)
    return np.transpose(anomaly, axes=(1,2,3,0,4,5))


# -----------------------------------------------------------------------------------------------------------
def get_sample(aSrc, gridArea):  # aSrc (nyear, ...)
    aSrc = ma.masked_equal(aSrc, 1e+20)  # (nyear, nwindow, ny, nx) or (nyear, nwindow, ngcm, nghm, ny, nx)
    #grid_area = array([area.compressed() for area in ma.masked_array(np.resize(gridArea, aSrc.shape), mask=aSrc.mask)])  # (nYear, nlandgrids)
    if len(aSrc.shape) == 6:
        mask = aSrc.mask[0,0,0,0]
        print(' @get_sample (in)\n aSrc.shape: {}\n mask.shape: {}'.format(aSrc.shape,mask.shape))
        grid_area = array([ma.masked_array(gridArea, mask=mask).compressed() for _win in range(window) 
                                                                                for _gcm in gcms 
                                                                                    for _ghm in ghms]).flatten()  # (nsample) nsample = window * ngcm * nghm * nlandgrids
    elif len(aSrc.shape) == 4:
        mask = aSrc.mask[0,0]
        print(' @get_sample (in)\n aSrc.shape: {}\n mask.shape: {}'.format(aSrc.shape,mask.shape))
        grid_area = array([ma.masked_array(gridArea, mask=mask).compressed() for _win in range(window)]).flatten()  # (nsample) nsample = window * nlandgrids
    elif len(aSrc.shape) == 3:
        grid_area = ma.masked_array(gridArea, mask=mask).compressed()
    aSrc = array([src.compressed() for src in aSrc])                                        # (nYear, nsample)  nsample = window * nlandgrids
    print(' @get_sample (out)\n aSrc.shape: {}\n grid_area.shape: {}'.format(aSrc.shape,grid_area.shape))
    return aSrc, grid_area


# -----------------------------------------------------------------------------------------------------------
def window_sampling(src):
    return array([src[i:i+window] for i, iyear in enumerate(range(syear,eyear+1))])  # (nYear, window, ...)


# ----------------------------------------------------------------------------------------------------------------------
def calc_probability(kdf):
    return kdf(days)


def calc_yearly_spatialstats(aSrc, smplType, gridarea, dict_info):   # aSrc: (nyear,???), gridarea (ngrid)
    """
    calculates regional/global stats.
    :param aSrc: array. stats samples.
    :param smplType: string. regional stats type >> Mean, Median, Average, Mode, Extreme
    :return aSrc: array. yearly spatial stats. shp=(nyear)
    :return PDFs: array of probability. shp=(nday, nyear)
    :return hist_min: single float value. historical range minimum
    :return hist_max: single float value. historical range maximum
    :return historicalPercentileDay: single integer value. This is for Extreme.
    """
    print(' --- calc_yearly_spatialstats.\n  smplType: {}...'.format(smplType))
    for item in ['regionID', 'region_name', 'ghm', 'gcm']: print('  {}: {}'.format(item, dict_info[item]))
    print('aSrc: {}\n gridarea: {}'.format(aSrc.shape, gridarea.shape))
    print(' -----------------------------')
    stime = datetime.datetime.now()

    regionID        = dict_info['regionID']
    region_name     = dict_info['region_name']
    ghm             = dict_info['ghm']
    gcm             = dict_info['gcm']
    Len             = dict_info['dict_parameters']['Len']
    index           = dict_info['dict_parameters']['index']
    tChunk          = dict_info['dict_parameters']['tChunk']
    drought_paras   = dict_info['dict_parameters']['drought_paras']
    mmeType         = dict_info['dict_parameters']['mmeType']
    smplType        = dict_info['dict_parameters']['smplType']
    scn             = dict_info['dict_parameters']['scn']
    soc             = dict_info['dict_parameters']['soc']
    co2             = dict_info['dict_parameters']['co2']
    season          = dict_info['dict_parameters']['season']
    fig_directory   = dict_info['dict_parameters']['fig_directory']
    ncout_directory = dict_info['dict_parameters']['ncout_directory']
    if how_to_gen_member_pdf == 'reuse':
        ex_nc_directory = dict_info['dict_parameters']['ex_nc_directory']
    if which_hist_percentile == 'reuse' or which_member_timeseries == 'reuse':
        ex_excel_directory = dict_info['dict_parameters']['ex_excel_directory']
    histricalPercentileDay = None

    def reuse_PDFs(_ghm=None, _gcm=None):
        """
        use this function, if you would like to use extisting PDF data which was generated befere.
        :return: pdfs
        """
        if _ghm is None: _ghm = ghm
        if _gcm is None: _gcm = gcm
        if PDFs2NC_MEMBER:                                       src_directory = ncout_directory
        elif how_to_gen_ensemble_pdf == 'median of member pdfs': src_directory = ncout_directory  # for Fig2Fig15
        else:                                                    src_directory = ex_nc_directory
        netcdf_file = 'drought_pdfs_{}.{}.{}.{:02}.{}.nc4'.format(drought_paras, _ghm, _gcm, regionID, region_name)
        netcdf_path = os.path.join(src_directory, netcdf_file)
        print('reuse_PDFs: {}'.format(netcdf_path))
        return Dataset(netcdf_path)['pdfs'][:].T  # (nday, nyear)

    def reuse_timeseries(_ghm=None, _gcm=None):
        if _ghm is None: _ghm = dict_info['ghm']
        if _gcm is None: _gcm = dict_info['gcm']
        # read aSrc and hist_pctl_day from an Excel file
        excel_name = 'TPCDs.{}.{:03}_{}.{}.{}.{}.{}_{}_{}.xlsx'.format(index, which_tChunk_reuse, which_threshold_reuse, 
                                                                    drought_paras, mmeType, smplType, scn, soc, co2)
        excel_path = os.path.join(ex_excel_directory, excel_name)
        row_name = '{}_{}_{:02}.{}'.format(ghm, gcm, regionID, region_name)
        return pd.read_excel(excel_path, sheet_name='data', index_col=0, usecols=[0]+list(range(8,8+nyear))).loc[row_name,:].values


    def reuse_hist_pctl():
        """
        read hist_pctl_day from an Excel file
        :return: hist_pctl_day
        """
        excel_name = 'TPCDs.{}.{:03}_{}.{}.{}.{}.{}_{}_{}.xlsx'.format(index, which_tChunk_reuse, which_threshold_reuse, 
                                                                 drought_paras, mmeType, smplType, scn, soc, co2)
        excel_path = os.path.join(ex_excel_directory, excel_name)
        row_name = '{}_{}_{:02}.{}'.format(ghm, gcm, regionID, region_name)
        return pd.read_excel(excel_path, sheet_name='data', index_col=0).loc[row_name, 'hist_90percentile_day']


    def calc_pdf(KDFs):
        """
        calculates probability of each value.
        :param KDFs: list of Kernel Density Functions.
        :return PDFs: array.
        """
        if multiprocessing:
            print('calculating PDFs of each NDD for each year  [MP={}]... This process will take time...'.format(multiprocessing))
            p = Pool(nMP)
            PDFs = np.array(list(p.map(calc_probability, KDFs))).T
            p.close()
            print('multiprocessing calc_probability process finished :)')
        else:
            PDFs = array([kdf(days) for kdf in KDFs]).T
        return PDFs  # (nday,nyear)


    def nan_exception(src):
        if np.nan in src:
            print('nan_exception:', src)
            return np.nan
        else:
            return src.tolist().index(np.max(src))


    def get_kdf(src, gridarea):
        ### without this process, we get "numpy.linalg.linalg.LinAlgError: singular matrix"...
        if (src[0] == 0) and len(list(set(src))) == 1:   src[-1] = 1e-20
        if (src[0] == 365) and len(list(set(src))) == 1: src[-1] = 365+1e-4
        return weighted_gaussian_kde(src, weights=np.resize(gridarea,src.shape))  # weighted by grid area in Karnel density estimate

    def getNearestValueIndex(List, num):
        return np.abs(np.array(List)-num).argmin()


    def get_percentileDay(kdf):
        _percentile_hist = percentile_hist * 1e-2
        cdmax = kdf.integrate_box_1d(0,365)            # total cumelative density
        cds = [kdf.integrate_box_1d(0,day) for day in days]  # list of cumelative density from 0 to x.
        cd_rates = np.array(cds)/cdmax                  # list of rate
        return getNearestValueIndex(cd_rates, _percentile_hist)  # int


    def get_percentileDay2(pdf):
        _percentile_hist = percentile_hist * 1e-2
        cdmax = pdf.sum()
        cds = [pdf[:days.index(day)].sum() for day in days]
        cd_rates = np.array(cds)/cdmax
        return getNearestValueIndex(cd_rates, _percentile_hist)  # int


    # main ---------------------------------------------------------------

    # KDE and PDF
    if smplType in ['Average', 'Mode', 'Extreme']:
        if ghm == 'all' and gcm == 'all':
            if how_to_gen_ensemble_pdf == 'create from all members':
                KDFs = [get_kdf(src, gridarea) for src in aSrc]
                PDFs = calc_pdf(KDFs)  # (nday,nyear2)
                del KDFs  # TODO: is this necessary in a function??
            elif how_to_gen_ensemble_pdf == 'median of member pdfs':
                PDFs = np.median([[reuse_PDFs(_ghm, _gcm) for _ghm in ghms] for _gcm in gcms], axis=(0,1))
            elif how_to_gen_ensemble_pdf == 'reuse':
                PDFs = reuse_PDFs()  # (nday, nyear2)
        else:  # for members
            if how_to_gen_member_pdf == 'newly create':  # newly create
                KDFs = [get_kdf(src, gridarea) for src in aSrc]
                PDFs = calc_pdf(KDFs)  # (nday,nyear2)
                del KDFs  # TODO: is this necessary in a function??
            elif how_to_gen_member_pdf == 'reuse':
                PDFs = reuse_PDFs()  # (nday, nyear2)
    else:
        print('PDFs...')
        N = aSrc.shape[1]
        _days = range(dictSeason[season]+1)
        PDFs = np.array([[np.where(_src==_day,1,0).sum() for _day in _days] for _src in aSrc]) / N  # (nyear, nday)
        PDFs = PDFs.T  # (nday, nyear)
        #PDFs = None


    # calcurate stats and gen a timeseries
    print(' --- calcurate stats and gen a timeseries')
    if ghm == 'all' and gcm == 'all' and how_to_plot_ensemble_timeseries == 'median of member timeseries':
        if which_member_timeseries == 'newly create':  # use the value the current df_stats. reuse time series in df which is made in this process.
            timeseries = [[dict_info['df_stats'][years].loc['{}_{}_{:02}.{}'.format(_ghm, _gcm, regionID, region_name)].values for _ghm in ghms] for _gcm in gcms]
        elif which_member_timeseries == 'reuse':  # load from an exsting excel file. reuse time series made in previous data processing...
            timeseries = [[reuse_timeseries(_ghm=_ghm, _gcm=_gcm) for _ghm in ghms] for _gcm in gcms]
        timeseries = np.median(timeseries, axis=(0,1))  # (nyear2)
    elif ghm == 'all' and gcm == 'all' and how_to_plot_ensemble_timeseries == 'reuse':
        timeseries = reuse_timeseries()
    elif ghm != 'all' and gcm != 'all' and which_member_timeseries == 'reuse':
        timeseries = reuse_timeseries()
    else:  # extract from a pdf... calculate one by one here...
        # generate a statistic time-series...
        if smplType == 'Mean':  # area-weighted spatial mean
            area_weight = gridarea / np.sum(gridarea)  # (ngrid)
            timeseries = np.sum(area_weight*aSrc, axis=1)  # (nyear)
            #timeseries = mean(aSrc, axis=1)  # (_nYear)
        elif smplType == 'Median':  # spatial median
            timeseries = median(aSrc, axis=1)  # (_nYear)
        #elif smplType == 'Mode': # spatial mode  TODO
        #    bottomT = Len
        #    topT = 40
        #    SS = range(bottomT, 365+1-topT)          # search scope for mode (cut out ~bottomT and 325~)
        #    print("Let's find yearly mode.\nobtaining mode (MP={})... This process takes time..".format(multiprocessing))
        #    if multiprocessing:
        #        p = Pool(nMP)
        #        aSrc = np.array(list(p.map(wrapper_calc_probability, itertools.product(KDFs,SS)))).reshape(-1, len(SS)) 
        #        p.close() 
        #        timeseries = array([bottomT + nan_exception(src) for src in aSrc])
        #    else:
        #        timeseries = array([bottomT + kdf(SS).tolist().index(np.max(kdf(SS))) for kdf in KDFs])  # (_nYear)
        elif smplType == 'Average':  # spatial average with KDE weight
            _days = days[Len:]
            Weights = PDFs.T  # (nyear, nday)
            timeseries = array([np.divide(np.sum(_days*weight[Len:]),(weight[0]+np.sum(weight[Len:]))) for weight in Weights])  # _nYear
            # Caution: Because NDD between 1-Len never occurs, this analysis excludes the range.
        elif smplType == 'Extreme':
            # TODO: Just a tentative approach!!!!!
            if which_hist_percentile == 'reuse':
                histricalPercentileDay = reuse_hist_pctl()
                print('histricalPercentileDay: {}'.format(histricalPercentileDay))
            elif which_hist_percentile == 'new':
                histricalPercentileDay = max([get_percentileDay2(pdf) for pdf in PDFs.T[si:ei]])
            i_histricalPercentileDay = range(366).index(histricalPercentileDay )
            timeseries = PDFs[i_histricalPercentileDay:,:].sum(axis=0)  # (nyear2)
            """    Note: Tentatively, this process is replaced with a simple way...
            if multiprocessing:
                def wrapper_calc_ExtremeArea(args):
                    def calc_ExtremeArea(kdf, histricalPercentileDay):
                        return kdf.integrate_box_1d(histricalPercentileDay, 365)
                    return calc_ExtremeArea(*args)
                print('get_percentileDay... (MP={})'.format(multiprocessing))
                p = Pool(nMP)
                histricalPercentileDay = max(list(p.map(get_percentileDay, KDFs[si:ei])))
                p.close()
                print('histricalPercentileDay: {}'.format(histricalPercentileDay))
                p = Pool(nMP)
                timeseries = array(list(p.map(wrapper_calc_ExtremeArea, itertools.product(KDFs, [histricalPercentileDay]))))
                p.close()
            else:
                print('get_percentileDay...')
                histricalPercentileDay = max([get_percentileDay(kdf) for kdf in KDFs[si:ei]])
                print('histricalPercentileDay: {}'.format(histricalPercentileDay))
                timeseries = array([kdf.integrate_box_1d(histricalPercentileDay,365) for kdf in KDFs])
            """
    if MEMORY_USAGE_NOISY:
        print("\n\n====== in calc_yearly_spatialstats\n{}{: >45}{}{: >10}{}".format('|','Variable Name','|','  Size','|'))
        print(" -------------------------- ")
        for k, v in locals().items():
            if hasattr(v, 'size') and not k.startswith('_') and not isinstance(v,types.ModuleType):
                print("{}{: >45}{}{: >10}{}".format('|',k,'|',str(v.size),'|'))
            elif hasattr(v, '__len__') and not k.startswith('_') and not isinstance(v,types.ModuleType):
                print("{}{: >45}{}{: >10}{}".format('|',k,'|',str(len(v)),'|'))
    print('calc_yearly_spatialstats took {} minutes in total...'.format((datetime.datetime.now() - stime).seconds / 60))

    if LOUD: print('@ calc_yearly_spatialstats\n{}\n{}'.format(timeseries, PDFs))
    return timeseries, PDFs, histricalPercentileDay  # (_nYear), (nday,nyear), int


# ----------------------------------------------------------------------------------------------------------------------
def find_TPCD(aSrc, tChunk, dict_info, dict_threshold_rcp85):
    """
    To find TPCD from a time series of regional stats.  for each region, each member and ensemble.
    :param aSrc: ndarray (nyear2)
    :param tChunk: int, single value
    :return tpcd:
    :return unprecedented_mask: True  = unprecedented!! = abnormal = more than the hiscorical max
                                False = normal within historical range.
    :return unexceptional_mask: True  = 
                                False = 
    :return hist_min:
    :return hist_max:
    """
    regionID      = dict_info['regionID']
    region_name   = dict_info['region_name']
    ghm           = dict_info['ghm']
    gcm           = dict_info['gcm']
    scn           = dict_info['scn']
    df_TPCD       = dict_info['df_TPCD']
    

    if ghm == 'all' and gcm == 'all' and how_to_estimate_ensemble_TPCD == 'median of member TPCDs':
        # use "already estimated member TPCD" from df_TPCD
        INDX = df_TPCD.index.tolist()[1:]  # exclude all_all
        print(df_TPCD.loc[INDX]['{}.{}'.format(regionID, region_name)])
        TPCD = int(df_TPCD.loc[INDX]['{}.{}'.format(regionID, region_name)].median())
        if TPCD > eyear_tpcdmap: TPCD = 1e+20
        unprecedented_mask, unexceptional_mask, threshold, hist_min, hist_max, n_event, n_year = None, None, None, None, None, None, None

    else:
        # --- define a histrical variation range
        hist_min, hist_max = aSrc[si:ei+1].min(), aSrc[si:ei+1].max()
        hist_mean = aSrc[si:ei+1].mean()
         
        if eyear_analysis_historical > 2005 and scn == 'rcp26': 
            threshold = dict_threshold_rcp85['{}_{}_{:03}.{}'.format(ghm, gcm, regionID, region_name)]
        else:
            if threshold_type == 'max':
                threshold = hist_max
            elif 'percentile' in threshold_type:
                threshold = np.percentile(aSrc[si:ei+1], float(threshold_type.replace('percentile','')))
            elif 'sigma' in threshold_type:
                factor = int(threshold_type.replace('sigma',''))
                threshold = hist_mean + factor*np.std(aSrc[si:ei+1])


        # --- detect unprecedented (find worsen cases)
        #   True:  unprecedented!! = abnormal = more than the hiscorical max
        #   False: normal, within historical range.
        unprecedented_mask = np.ma.make_mask(aSrc > threshold)
        unprecedented = ''.join(['F' if i == b'F' else 'T' for i in unprecedented_mask.astype('S1')])  # 'FFFFFFTTTTTTFFTTTFFFFFFFF...'   'T' means 'unprecedented'
        if len(unprecedented) == 1: 
            unprecedented = 'F'*aSrc.shape[0]

        # --- Caution!! At this point, the first item in "unprecedented" is NOT always F in case you do not select the entire historical preiod as the reference period.
        if unprecedented[0] == 'F': 
            initial_flug = 'F'
        elif unprecedented[0] == 'T': 
            initial_flug = 'T'

        # --- split the seaquence
        Fs = pF.findall(unprecedented)  # ['FFFFFF', 'FF', 'FFFFFFFF', ...]
        Ts = pT.findall(unprecedented)  # ['TTTTTT', 'TTT, ... ]

        # --- remove shorter unprecedented period than tChunk
        Ts = ['F'*len(ts) if len(ts) < tChunk else ts for ts in Ts]  # ['TTTTTT', 'FFF', ...]
        
        # merge
        if initial_flug == 'F':
            item1, item2 = Fs, Ts
        elif initial_flug == 'T': 
            item1, item2 = Ts, Fs
        if len(item1) > len(item2): item2.append('')
        unprecedented = []
        for _item1, _item2 in zip(item1, item2):
            unprecedented.append(_item1)
            unprecedented.append(_item2)
        unprecedented_mask = [False if i == 'F' else True for i in list(''.join(unprecedented))]  # [False, ..., True, True, True, ...]

        # --- find TPCD
        # --- get all non-future years False.
        _unprecedented_mask = [False]*(years.index(years_analysis_historical[-1])+1) + unprecedented_mask[years.index(years_analysis_future[0]):]
        if True in _unprecedented_mask:
            # --- find an index of the first True in the future period
            TPCD = int(years[_unprecedented_mask.index(True)])
        else:
            TPCD = 1e+20
        # --- limit with tChunk size
        if TPCD > (2100-tChunk):
            TPCD = 1e+20

        # --- make a mask to calculate the potential population affected under the unprecedented drought condition.
        #     the mask maskes out NON unprecedented grids. -->  If not unprecedented, then this mask is True.
        #     shorter unprecedented condition than tChunk are removed here.
        unexceptional_mask = np.full(aSrc.shape, False)
        unexceptional_mask[np.where(np.array(unprecedented_mask)==False)] = True  # [True ... False False False ...]

        _Ts = [len(ts) for ts in Ts if 'T' in ts]
        n_event = len(_Ts)
        n_year = sum(_Ts)

    if LOUD:
        print('@find_TPCD')
        print('TPCD     : {}'.format(TPCD))
        print('threshold: {}'.format(threshold))
    return TPCD, unprecedented_mask, unexceptional_mask, threshold, hist_min, hist_max, n_event, n_year


# -----------------------------------------------------------------------------------------------------------
def create_df_TPCD(dictRegion, RegionIDs):  # to create an Excel sheet for TPCD
    CLMNs = ['GHM', 'GCM'] + ['{}.{}'.format(regionID,dictRegion[regionID]) for regionID in RegionIDs]#; print(CLMNs)
    INDXs = ['all_all'] + ['{}_{}'.format(ghm,gcm) for ghm, gcm in itertools.product(ghms,gcms)]#; print(INDXs)
    df_TPCD = pd.DataFrame(index=INDXs, columns=CLMNs)
    for clmn in ['GHM', 'GCM']: df_TPCD[clmn] = df_TPCD[clmn].astype(str)
    df_TPCD['GHM']['all_all'] = 'all'
    df_TPCD['GCM']['all_all'] = 'all'
    for ghm, gcm in itertools.product(ghms,gcms): 
        df_TPCD['GHM']['{}_{}'.format(ghm,gcm)] = ghm
        df_TPCD['GCM']['{}_{}'.format(ghm,gcm)] = gcm
    # substitute missing value, as an initial value just for now
    for CLMN in CLMNs[2:]:
        df_TPCD[CLMN]['all_all'] = 1e+20
    for CLMN, ghm, gcm in itertools.product(CLMNs[2:], ghms, gcms):
        df_TPCD[CLMN]['{}_{}'.format(ghm,gcm)] = 1e+20
    return df_TPCD


# ----------------------------------------------------------------------------------------------------------
def create_df_stats(dictRegion, RegionIDs):  # to create an Excel sheet for timeseries information with TPCD and historical min/max value
    info_CLMNs = ['GHM', 'GCM', 'Region', 'TPCD', 'threshold', 'min', 'max', 'hist_{}percentile_day'.format(percentile_hist)]
    CLMNs = info_CLMNs + list(years)  # 1973-2097
    INDXs = ['{}_{}_{:02}.{}'.format(ghm, gcm, regionID,  dictRegion[regionID]) for regionID in RegionIDs for ghm in ghms for gcm in gcms]
    df_stats = pd.DataFrame(index=INDXs, columns=CLMNs)
    for clmn in ['GHM', 'GCM', 'Region']: 
        df_stats[clmn] = df_stats[clmn].astype(str) 
    return df_stats

def add2df_stats(df_stats, ghm, gcm, rgn, tpcd, threshold, hist_min, hist_max, hist_percentile_day, aSrc):

    idx = '{}_{}_{}'.format(ghm, gcm, rgn)

    df_stats.loc[idx,'GHM'] = ghm
    df_stats.loc[idx,'GCM'] = gcm
    df_stats.loc[idx,'Region'] = rgn
    df_stats.loc[idx,'TPCD'] = tpcd
    df_stats.loc[idx,'threshold'] = threshold
    df_stats.loc[idx,'min'] = hist_min
    df_stats.loc[idx,'max'] = hist_max
    df_stats.loc[idx,'hist_{}percentile_day'.format(percentile_hist)] = hist_percentile_day
    for y, year in enumerate(years): 
        df_stats.loc[idx, year] = aSrc[y]


# -----------------------------------------------------------------------------------------------------------
def add2df_unprecedented_flag_member(df, rgn, ghm, gcm, flag):
    idx = '{}_{}_{}'.format(rgn, ghm, gcm)
    df.loc[idx, 'region'] = rgn
    df.loc[idx, 'ghm'] = ghm
    df.loc[idx, 'gcm'] = gcm
    if not len(flag) == 0:
        for y, year in enumerate(years):
            df.loc[idx, year] = 1 if flag[y] else 0
    else:
        for y, year in enumerate(years):
            df.loc[idx, year] = 0


def add2df_unprecedented_info_member(df, rgn, ghm, gcm, n_event, n_year):
    idx = '{}_{}_{}'.format(rgn, ghm, gcm)
    df.loc[idx, 'region'] = rgn
    df.loc[idx, 'ghm'] = ghm
    df.loc[idx, 'gcm'] = gcm
    df.loc[idx, 'n_event'] = n_event
    df.loc[idx, 'n_year'] = n_year


def add2df_unprecedented_info_median(df, rgn, n_event, n_year, df_member):
    if n_event is not None and n_year is not None:
        df.loc[rgn, 'n_event'] = n_event
        df.loc[rgn, 'n_year'] = n_year
    else:
        df.loc[rgn, 'n_event'] = df_member[df_member.region==rgn]['n_event'].median(axis=0)
        df.loc[rgn, 'n_year'] = df_member[df_member.region==rgn]['n_year'].median(axis=0)


# -----------------------------------------------------------------------------------------------------------
def gen_vmax(pdf):
    """
    estimate vmax of a pdf within a nday range.
    :param pdf: array. (nday)
    :return vmax: single float value.
    """
    #print('gen_vmax...')
    dp = 0.0001
    head = 30
    tail = 366-5
    #tail_max = int(365*0.8)
    for_simple_decrease = int(365*0.15)
    x_threshold = int(365*0.25)

    start_value = dp * (pdf.max() // dp)  # integer
    criteria_series = np.arange(start_value, 0, -dp)                                                     # (ncriterion) from the largest value
    #pdf = pdf[:tail_max]
    #pdf = pdf[head:]
    pdf = pdf[head:tail]
    n_intersection = np.array([sum([1 if pdf[i:i+2].min() < criteria < pdf[i:i+2].max() else 0 for i in range(pdf.shape[0])]) for criteria in criteria_series]) # (ncriterion) 

    if not np.where(n_intersection>=2)[0].shape[0] == 0:   
        index_first_multi_intersection = np.where(n_intersection>=2)[0][0]
        vmax_tmp = criteria_series[index_first_multi_intersection]
        vmax = max(vmax_tmp, pdf[x_threshold])
        vmax_type = 'type1'
    else:  # a convex downword curve
        #vmax = pdf[for_simple_decrease]
        #vmax = np.percentile(pdf, 70)  # TODO. Just for now. a test...
        #vmax = np.percentile(pdf, 50)  # TODO. Just for now. a test...
        #vmax = np.percentile(pdf, 30)  # TODO. Just for now. a test...
        #vmax = np.percentile(pdf, 20)  # TODO. Just for now. a test...
        vmax = np.percentile(pdf, 10)  # TODO. Just for now. a test...
        #vmax = np.percentile(pdf, 5)  # TODO. Just for now. a test...
        #vmax = np.percentile(pdf, 2)  # TODO. Just for now. a test...
        #vmax = np.percentile(pdf, 1)  # TODO. Just for now. a test...
        #vmax = np.percentile(pdf, 0.5)  # TODO. Just for now. a test...
        vmax_type = 'type2'

    return vmax, vmax_type


# -----------------------------------------------------------------------------------------------------------
def update_dict_stats(dict_stats, scn, soc, co2, region_name, ghm, gcm):  # TODO   modified
    if not scn in dict_stats:                                  dict_stats[scn] = {}
    if not soc in dict_stats[scn]:                             dict_stats[scn][soc] = {}
    if not co2 in dict_stats[scn][soc]:                        dict_stats[scn][soc][co2] = {}
    if not region_name in dict_stats[scn][soc][co2]:           dict_stats[scn][soc][co2][region_name] = {}
    if not ghm in dict_stats[scn][soc][co2][region_name]:      dict_stats[scn][soc][co2][region_name][ghm] = {}
    if not gcm in dict_stats[scn][soc][co2][region_name][ghm]: dict_stats[scn][soc][co2][region_name][ghm][gcm] = {}
    return dict_stats

def allocate_item_2_dict_stats(dict_stats, key, item, dict_info):  # TODO
    scn         = dict_info['scn']
    soc         = dict_info['soc']
    co2         = dict_info['co2']
    region_name = dict_info['region_name']
    ghm         = dict_info['ghm']
    gcm         = dict_info['gcm']
    dict_stats[scn][soc][co2][region_name][ghm][gcm][key] = item
    return dict_stats

def draw_colorfulplot_and_shiftingPDF(dict_stats, dict_info):
    """
    creates the colorfulplot_shiftingPDF figure.
    :param dict_stats: dictionary.
                key)
                    pdfs (ndarray, (nday,nyear))
                    spatial_average (array)
                    historical_percentile (single integer)
    :param year_index: integer. a target year for sample colorplot
    """
    scn           = dict_info['scn']
    soc           = dict_info['soc']
    co2           = dict_info['co2']
    regionID      = dict_info['regionID']
    region_name   = dict_info['region_name']
    ghm           = dict_info['ghm']
    gcm           = dict_info['gcm']
    sample_year   = dict_info['sample_year']
    qvalType      = dict_info['dict_parameters']['qvalType']
    index         = dict_info['dict_parameters']['index']
    drought_paras = dict_info['dict_parameters']['drought_paras']
    fig_directory_common = dict_info['fig_directory_common']
    print('\ndraw_colorfulplot_and_shiftingPDF: {}, {}, {}, {}, {}, {}, {}'.format(scn, soc, co2, region_name, ghm, gcm, sample_year))

    fig_directory = os.path.join(fig_directory_common, today)

    fig_id_loc_x, fig_id_loc_y = -0.25, 1.03
    #ndd_title = 'Frequency of drought days (FDD) in a year'
    ndd_title = 'FDD'
    labelsize = 7.5
    unitfontsize = 7.5

    colorbar_position = 'lower_horizontal'
    # colorbar_position = 'center_vertical'
    # colorbar_position = 'right_vertical'

    exremecolor = '#808080'
    # exremecolor = '#dd0000'

    x_scale = 3  # to get the map wider    for 1971-2099 (129)
    #x_scale = 2  # to get the map wider    for 1865-2099 (235)
    y_scale = 1.1
    #ticklabels_days = [365, 300, 250, 200, 150, 100, 50, 0]  # some days
    ticklabels_days = [100, 80, 60, 40, 20, 0]  # [%]
    #ticklabels_years = [1975, 2000, 2030, 2060, 2090]
    #ticklabels_years = [1900, 1950, 2000, 2050]
    ticklabels_years = [syear_historical_for_plot, 2013, 2030, 2050, 2080]

    _days = range(365,-1,-1)  # [365, 364, ..., 2, 1, 0]
    year_index = years.index(sample_year)

    pdfs = dict_stats['PDFs']  # (nday,nyear)
    pdf = pdfs[:, year_index]  # pdf of the year (nday)
    print('pdfs.shape: {}'.format(pdfs.shape))
    print('pdf.shape:  {}  for {}'.format(pdf.shape, sample_year))
    spatial_average = round(dict_stats['spatial_average'][year_index])  # single integer value
    #historical_percentile = round(dict_stats['historical_percentile'])  # single integer value

    nday = pdf.shape[0]
    xs = range(nday)  # [0, 1, 2, ..., 365]

    vmin = pdf.min()
    #vmax, dummy = gen_vmax(pdf)  # generate vmax based on the target pdf...
    vmax, dummy = gen_vmax_overall(pdfs)
    # create a continuous "shared" norm to map from data points to colors
    norm = plt.Normalize(vmin=vmin, vmax=vmax)

    if   colorbar_position == 'center_vertical':  figsize = (10, 4);    left = 0.08; bottom = 0.1;  right = 0.96; top = 0.9; wspace = 0.3
    elif colorbar_position == 'right_vertical':   figsize = (10, 4);    left = 0.1;  bottom = 0.1;  right = 0.9;  top = 0.9; wspace = 0.25
    #elif colorbar_position == 'lower_horizontal': figsize = (8.5, 4.2); left = 0.02; bottom = 0.27; right = 0.98; top = 0.95; wspace = 0.3
    elif colorbar_position == 'lower_horizontal': figsize = (8.5, 5); left = 0.10; bottom = 0.27; right = 0.98; top = 0.90; wspace = 0.1


    # draw a figure ---------------------------------------------------------------------------------------------------------------
    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=figsize)
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace)
    # --- (left) colorful plot
    print('drawing ax1 figure...')
    ax1.text(fig_id_loc_x, fig_id_loc_y, '(a)', fontsize=10, transform=ax1.transAxes)
    # create a set of line segments
    points = np.array([xs, pdf]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    # color gradation plot
    lc = LineCollection(segments, cmap=pdf_cmap, norm=norm)
    lc.set_array(pdf)
    lc.set_linewidth(2)
    line = ax1.add_collection(lc)
    ax1.text(0.75, 0.90, 'Sample:\n@Year={}'.format(sample_year), fontsize=labelsize, transform=ax1.transAxes)
    # add an average line
    average_index = xs.index(spatial_average)
    ax1.plot([xs[average_index], xs[average_index]], [0, pdf[average_index]], linewidth=0.2, color='#000000')
    ## add extreme shade
    #percentile_index = xs.index(historical_percentile)
    #xs_percentile = xs[percentile_index:]
    #ax1.fill_between(xs_percentile, pdf[percentile_index:], 0, facecolor=exremecolor)  # , alpha=0.5)
    # figure parameters
    ax1.set_xlabel(ndd_title)
    ax1.text(0.95, -0.11, '[{}]'.format(kernel_randomvariable_unit), ha='left', va='center', fontsize=unitfontsize, transform=ax1.transAxes)
    ax1.set_xticks([xs.index(round(365*i*1e-2)) for i in ticklabels_days])
    ax1.set_xticklabels(ticklabels_days)
    ax1.set_xlim(xmin=0)  # TODO: this is just for the test
    #ax1.set_ylabel('Probability')
    ax1.set_ylabel('Probability of FDD')
    ax1.set_ylim([0, vmax * y_scale])
    ax1.tick_params(axis='both', which='major', labelsize=labelsize)
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    # --- (right) shifting PDFs
    print('drawing ax2 figure...')
    ax2.text(fig_id_loc_x-0.015, fig_id_loc_y, '(b)', fontsize=10, transform=ax2.transAxes)
    # create src (get the pdf wider)
    src = np.array([pdfs[:, year_index] for year_index in range(nyear) for repitation in range(x_scale)]).T  # (nday, x_scale*nyear)
    src = src[::-1]
    # --- draw a shifting PDFs
    im = ax2.imshow(src, vmin=vmin, vmax=vmax, cmap=pdf_cmap)
    ax2.axvline(x=x_scale*year_index, linewidth=0.4, linestyle='--', color='#ffffff')
    
    #### plot regional proability-weighted timeseries  TODO: plots needs to be region-member specific. currently, ensemble time series is used... 
    ###xs4timeseries = x_scale * np.arange(0,len(dict_stats['spatial_average'])) + x_scale/2
    ###ax2.plot(xs4timeseries, dict_stats['spatial_average'], linewidth=0.4, linestyle='-', color='#d0d0d0', label='The regional average FDD')
    ###ax2.legend(bbox_to_anchor=(1, 1.05), loc='lower right', borderaxespad=0, fontsize=labelsize, frameon=False)
    # figure parameters
    xticklabels = ticklabels_years
    xticks = [x_scale * years.index(year) for year in xticklabels]
    xticklabels.append(sample_year)
    yticklabels = ticklabels_days
    yticks = [_days.index(round(365*percent*1e-2)) for percent in yticklabels]
    print('figure parametes...')
    ax2.set_xlabel('Year')
    ax2.set_xticks(xticks)
    ax2.set_xticklabels(xticklabels)
    ax2.set_ylabel(ndd_title)
    ax2.text(-0.10, 1.03, '[{}]'.format(kernel_randomvariable_unit), ha='left', va='bottom', fontsize=unitfontsize, transform=ax2.transAxes)
    ax2.set_yticks(yticks)
    ax2.set_yticklabels(yticklabels)
    ax2.tick_params(axis='both', which='major', labelsize=labelsize)
    # --- annotation
    to_position_x = x_scale*years.index(sample_year)
    if sample_year <= 2040: dx = 15
    else:                   dx = -15
    ax2.annotate('Sample:\n@Year={}'.format(sample_year), xy=(to_position_x,0), xytext=(to_position_x+dx, -20),
                 xycoords='data', textcoords='data', size=labelsize, va='center', ha='left',
                 arrowprops = dict(arrowstyle='->', connectionstyle='arc3,rad=+0.2', fc='w'),
                 )
    # --- colorbar
    if 'vertical' in colorbar_position:
        if colorbar_position == 'center_vertical':
            divider = make_axes_locatable(ax1); pad = 0.2
        elif colorbar_position == 'right_vertical':
            divider = make_axes_locatable(ax2); pad = 0.05
        ax_cb = divider.new_horizontal(size="2%", pad=pad)
        fig.add_axes(ax_cb)
        clb = plt.colorbar(im, cax=ax_cb)
        clb.ax.tick_params(labelsize=labelsize)
        if colorbar_position == 'center_vertical':
            #clb.ax.set_title('Probability', fontsize=10)
            clb.ax.set_title('Probability of FDD', fontsize=10)
        elif colorbar_position == 'right_vertical':
            #clb.set_label('Probability', rotation=270, labelpad=10)
            clb.set_label('Probability of FDD', rotation=270, labelpad=10)
    elif colorbar_position == 'lower_horizontal':
        ax_cb = fig.add_axes([0.3, 0.115, 0.4, 0.02])
        clb = plt.colorbar(im, cax=ax_cb, orientation='horizontal')
        # clb = fig.colorbar(im, cax=ax_cb, orientation='horizontal')
        clb.ax.tick_params(labelsize=labelsize)
        #clb.set_label('Probability', size=10, labelpad=4)
        clb.set_label('Probability of FDD', size=10, labelpad=4)

    if ghm == 'all' and gcm == 'all': fig_num = '01'
    else:                             fig_num = '02'
    figname = 'colorfulplot_shiftingPDS.{}_{}_{}.{:02}{}.{}_{}_{}_{}'.format(scn, soc, co2, regionID, region_name, fig_num, ghm, gcm, sample_year)
    for suffix in suffixes:
        _fig_directory = os.path.join(fig_directory, 'colorfulplot_shiftingPDF', suffix[1:])
        if not os.path.isdir(_fig_directory): os.makedirs(_fig_directory)
        figpath = os.path.join(_fig_directory, figname+suffix)
        plt.savefig(figpath, dpi=dpi_savefig)
        if suffix == suffixes[0]: print('savefig: {}\n'.format(figpath))
    plt.close(3)


# -----------------------------------------------------------------------------------------------------------
def create_figureobject1(num):  # Kernel timeseries on basemap

    print('Just a moment. Preparing figure object_1 (num={})...'.format(num))
    fig = plt.figure(num=num, figsize=(10,5), dpi=dpi_figure)
    print('adding ax...')
    ax = fig.add_subplot(1,1,1)
    plt.subplots_adjust(left=0.02, bottom=0.02, right=0.98, top=0.98)
    ax.axis('off')
    ax_pos = ax.get_position()

    if fig1_proj == 'robin': bm = Basemap(projection='robin',lon_0=0, ax=ax); bm.drawmapboundary(linewidth=0)
    elif fig1_proj == 'cyl': bm = Basemap(projection='cyl',llcrnrlat=-90,llcrnrlon=-180,urcrnrlat=90,urcrnrlon=180,resolution='l',ax=ax)
    #elif fig1_proj == 'cyl': bm = Basemap(projection='cyl',llcrnrlat=-60,llcrnrlon=-180,urcrnrlat=90,urcrnrlon=180,resolution='l',ax=ax)
    else: print('check fig1_proj...'); sys.exit()

    if fig1_backcolor == 'fillcontinents': bm.fillcontinents(color='sandybrown')
    elif fig1_backcolor == 'shadedrelief': bm.shadedrelief()
    else: print('check fig1_backcolor'); sys.exit()

    if region_type == 'AR6_regions':
        bm.readshapefile(region_shp_path, 'reference-regions-AR6', linewidth=0.3, color='k')
    elif region_type == 'SREX_regions':
        bm.readshapefile(region_shp_path, 'referenceRegions', linewidth=0.3, color='k')
    elif region_type == 'HydroBASINS_lev1':
        bm.readshapefile(region_shp_path, 'hybas_lake____lev02_v1c_merge', linewidth=0.3, color='k')
    elif region_type == 'hydroregions':
        bm.readshapefile(region_shp_path, 'meybeck_et_al_2013_hydroregions', linewidth=0.3, color='k')

    return fig, ax_pos 


def updatefig1_addaxes_and_imshow(PDFs, fig, ax_pos, kdfSampling, dict_info):

    print('preparing KDF matrix (nDay x years)...')
    rgnid = dict_info['regionID']
    rgn = dict_info['region_name']
    Len = dict_info['dict_parameters']['Len']

    ticks_ndd = [0, (365-Len)/4*1, (365-Len)/4*2, (365-Len)/4*3, 365-Len]
    if kernel_randomvariable_unit == 'day/year':
        ticklabels_ndd = [Len, 
                          Len+int((365-Len)/4*1), 
                          Len+int((365-Len)/4*2), 
                          Len+int((365-Len)/4*3), 365][::-1]
    elif kernel_randomvariable_unit == '%':
        ticklabels_ndd = [int(float(Len)/365*1e2), 
                          int((Len+float(365-Len)/4*1)/365*1e2), 
                          int((Len+float(365-Len)/4*2)/365*1e2), 
                          int((Len+float(365-Len)/4*3)/365*1e2), 100][::-1]
        
    #vmax           = dictIPCCregion[rgnid][2]     Caution!! No longer used!!
    dist_x, dist_y = dictIPCCregion[rgnid][1]
    if rgnid == 0: hight = 0.45
    else:          hight = 0.25
    norm = colors.Normalize(vmin=0, vmax=vmax)
    bounds = [0,vmax]
    aspect = 0.34

    if kdfSampling == 'window':
        xticklabels = [2006,2050]
        xticks = [range(syear, eyear+1).index(year) for year in xticklabels]
        xticklabels = ["'{}".format(str(year)[-2:]) for year in xticklabels]
    else:
        xticks = [0,35,79,129]
        xticklabels = [1971,2006,2050,2099]

    ax02 = fig.add_axes([ax_pos.x0+dist_x, ax_pos.y0+dist_y, aspect*hight, hight])
    ax02_pos = ax02.get_position()
    im = ax02.imshow(PDFs[::-1], norm=norm, vmin=0, vmax=vmax, cmap=pdf_cmap)
    if Len == 1: Len = 0
    ax02.set_yticks(ticks_ndd)
    ax02.set_xticks(xticks)

    if rgn == 'GLB': 
        ax02.text(0.01, 1.08, 'Global', va='top', ha='left', fontsize=8, transform=ax02.transAxes)
        ax02.set_xlabel('Year', fontsize=8)
        ax02.set_xticklabels(xticklabels, fontsize=7)
        ax02.set_ylabel(kernel_randomvariable_unit, fontsize=8)
        ax02.set_yticklabels(ticklabels_ndd, fontsize=8)
        ax02.xaxis.set_label_coords(0.5,-0.14)
        cax = fig.add_axes([ax02_pos.x0+0.038, ax02_pos.y0-0.05, ax02_pos.x1-ax02_pos.x0-0.076, 0.02])
    else:               
        ax02.text(0.01, 1.085, '{} {}'.format(rgnid,rgn), va='top', ha='left', fontsize=8, transform=ax02.transAxes)
        ax02.set_xticklabels(xticklabels, fontsize=0)
        ax02.set_yticklabels([])
        cax = fig.add_axes([ax02_pos.x0+0.021, ax02_pos.y0-0.014, ax02_pos.x1-ax02_pos.x0-0.047, 0.01])

    ax02.tick_params(axis='y', direction='in', length=3)
    ax02.tick_params(axis='x', direction='in', length=3)

    cb = fig.colorbar(im, cax=cax, orientation='horizontal')
    cb.set_ticks(bounds)
    cb.set_ticklabels(bounds)
    cb.ax.tick_params(labelsize=5, direction='in')

    #return fig


# -----------------------------------------------------------------------------------------------------------
def create_figureobject2(num):  # PDF time series in the list style

    # currently, this figure is not used...
    print('Caution!! You need to modify this function for AR6 regions!!')
    print('26 -> 44 regions')
    sys.exit()

    print('preparing a figure object_2 (num={}). Just a moment...'.format(num))

    # --- for a longer period
    #fig = plt.figure(num=num, figsize=(8,4.5), dpi=dpi_figure)
    fig = plt.figure(num=num, figsize=(11,4.5), dpi=dpi_figure)
    #plt.subplots_adjust(left=0.1, bottom=0.08, right=0.98, top=0.98)  #, wspace=0.06, hspace=0.44)
    plt.subplots_adjust(left=0.1, bottom=0.08, right=0.98, top=0.97)  #, wspace=0.06, hspace=0.44)

    # add the region category map, at first
    ax = plt.subplot2grid((3,10), (0,7), colspan=3)
    ax.axis('off')
    plt.text(0.5, -0.1, '(a) Region category', transform=ax.transAxes, fontsize=7, horizontalalignment='center')
    if fig2_proj == 'robin': bm = Basemap(projection='robin',lon_0=0, ax=ax); bm.drawmapboundary(linewidth=0)
    elif fig2_proj == 'cyl': bm = Basemap(projection='cyl',llcrnrlat=-90,llcrnrlon=-180,urcrnrlat=90,urcrnrlon=180,resolution='l',ax=ax)
    else: print('check fig2_proj...'); sys.exit()
    if fig2_backcolor == 'fillcontinents': bm.fillcontinents(color='sandybrown')
    elif fig2_backcolor == 'shadedrelief': bm.shadedrelief()
    else: print('check fig2_backcolor'); sys.exit()
    if region_type == 'AR6_regions':
        bm.readshapefile(region_shp_path, 'reference-regions-AR6', linewidth=0.3, color='k')
    elif region_type == 'SREX_regions':
        bm.readshapefile(region_shp_path, 'referenceRegions', linewidth=0.3, color='k')
    elif region_type == 'HydroBASINS_lev1':
        bm.readshapefile(region_shp_path, 'hybas_lake____lev02_v1c_merge', linewidth=0.3, color='k')
    elif region_type == 'hydroregions':
        bm.readshapefile(region_shp_path, 'meybeck_et_al_2013_hydroregions', linewidth=0.3, color='k')
    for iregion in range(1, n_regions+1):
        ix, iy = dictIPCCregion[iregion][2]   # position of regions name 
        plt.text(ix, iy, iregion, transform=ax.transAxes, fontsize=6)

    return fig


def read_vmax_overall(dict_info):

    scn                = dict_info['scn']
    soc                = dict_info['soc']
    co2                = dict_info['co2']
    regionID           = dict_info['regionID']
    region_name        = dict_info['region_name']
    index              = dict_info['dict_parameters']['index']
    drought_paras      = dict_info['dict_parameters']['drought_paras']
    mmeType            = dict_info['dict_parameters']['mmeType']
    smplType           = dict_info['dict_parameters']['smplType']
    ex_excel_directory = dict_info['dict_parameters']['ex_excel_directory']

    excel_name = 'vrange_of_PDFs.{}.{}.{}.{}.{}_{}_{}.xlsx'.format(index, drought_paras, mmeType, smplType, scn, soc, co2)
    excel_path = os.path.join(ex_excel_directory, excel_name)
    column_name = '{}.{}'.format(regionID, region_name)
    return pd.read_excel(excel_path, index_col=0)[column_name]['vmax']


def gen_vmax_overall(pdfs, region=None, df_vrange=None):
    """
    calculate appropriate vmax for a PDFs
    :param pdfs: ndarray (float)   shape=(nDayTot, nyear)
    :return _vmax: single float value
    """
    print('gen_vmax_overall...')
    nday, nyear = pdfs.shape
    yearly_vmaxs = []
    n_type1_counter,  n_type2_counter = 0, 0

    for iyear in range(nyear):
        vmax_iyear, vmax_type = gen_vmax(pdfs[:, iyear])
        yearly_vmaxs.append(vmax_iyear)
        if vmax_type == 'type1': n_type1_counter += 1
        elif vmax_type == 'type2': n_type2_counter += 1
        else: print('check vmax_type...'); sys.exit()
        if VRANGE_DETAIL_TO_EXCEL and region is not None and df_vrange is not None:
            df_vrange.loc['n_type1', region] = n_type1_counter
            df_vrange.loc['n_type2', region] = n_type2_counter
            df_vrange.loc[years[iyear], region] = vmax_iyear

    #vmax = np.median(yearly_vmaxs)
    #vmax = np.max(yearly_vmaxs)
    #vmax = np.percentile(yearly_vmaxs, 50)
    #vmax = np.percentile(yearly_vmaxs, 60)
    #vmax = np.percentile(yearly_vmaxs, 75)
    #vmax = np.percentile(yearly_vmaxs, 80)
    #vmax = np.percentile(yearly_vmaxs, 90)
    #vmax = np.percentile(yearly_vmaxs, 92)
    #vmax = np.percentile(yearly_vmaxs, 93)
    vmax = np.percentile(yearly_vmaxs, 94)
    #vmax = np.percentile(yearly_vmaxs, 95)
    #vmax = np.percentile(yearly_vmaxs, 99)
    #vmax = min(yearly_vmaxs)
    if vmax == 0: vmax = 0.004
    #vmax = round(vmax, 5)
    vmax = float(decimal.Decimal(str(vmax)).quantize(decimal.Decimal('0.0001'), rounding=decimal.ROUND_UP))

    return vmax, df_vrange


def updatefig2_addaxes_and_imshow(PDFs, fig, iregion, kdfSampling, dict_info, df_vrange):
    print('\nupdatefig2_addaxes_and_imshow...')
    print('preparing KDF matrix (nDay x years)...')

    regionID = dict_info['regionID']
    region_name = dict_info['region_name']
    Len = dict_info['dict_parameters']['Len']

    if kernel_randomvariable_unit == 'day/year':
        ticks_ndd = [0, (365-Len)/4*1, (365-Len)/4*2, (365-Len)/4*3, 365-Len]
        ticklabels_ndd = [Len, Len+int((365-Len)/4*1), Len+int((365-Len)/4*2), Len+int((365-Len)/4*3), 365][::-1]
    elif kernel_randomvariable_unit == '%':
        #ticklabels_ndd = [int(float(Len)/365*1e2), 
        #                  int((Len+float(365-Len)/4*1)/365*1e2), 
        #                  int((Len+float(365-Len)/4*2)/365*1e2), 
        #                  int((Len+float(365-Len)/4*3)/365*1e2), 
        #                  100][::-1]
        ticklabels_ndd = [25, 50, 75, 100]
        #ticks_ndd = [(365 - 365*i*1e-2)/365*100 for i in ticklabels_ndd]  # [%]  NG?
        ticks_ndd = [365*(1-i*1e-2) for i in ticklabels_ndd]  # [%]

    print('updatefig2 addaxes and imshow...')
    if iregion >= 7: ax = fig.add_subplot(3,10,iregion+1+3)
    else:            ax = fig.add_subplot(3,10,iregion+1)
    ax_pos = ax.get_position()
    if region_name == 'GLB': title = ax.set_title('Global', fontsize=7)
    else:                    title = ax.set_title('{} {}'.format(regionID,region_name), fontsize=7)
    title.set_position([0.5, 0.95])

    vmin = round(PDFs.min(), 4)
    if READ_VMAX_OVERALL:
        vmax = read_vmax_overall(dict_info)
    else:
        vmax, df_vrange = gen_vmax_overall(PDFs, '{}.{}'.format(regionID, region_name), df_vrange)
        df_vrange.loc['vmax', '{}.{}'.format(regionID, region_name)] = vmax
        df_vrange.loc['vmin', '{}.{}'.format(regionID, region_name)] = vmin
    norm = colors.Normalize(vmin=vmin, vmax=vmax)
    bounds = [vmin,vmax]

    im = ax.imshow(PDFs[::-1], norm=norm, vmin=vmin, vmax=vmax, cmap=pdf_cmap)
    
    # y-tick
    ax.tick_params(axis='y', direction='out', length=2)
    if Len == 1: Len = 0
    ax.set_yticks(ticks_ndd)
    if iregion == 0 or iregion == 7 or iregion == 17:
        ax.set_yticklabels(ticklabels_ndd, fontsize=6)
    else: ax.set_yticklabels(['']*len(ticks_ndd))
    ax.set_ylim([365-Len,0])
    if iregion == 0:  # y-axis label and unit
        #ax.set_ylabel('FDD')
        ax.text(0, 1, '[%] \n', va='bottom', ha='right', fontsize=7)

    # x-tick
    ax.tick_params(axis='x', direction='out', length=2)
    if kdfSampling == 'window':
        xticklabels = [2000,2050]
        xticks = [range(syear, eyear+1).index(year) for year in xticklabels]   # TODO
        if 17 <= regionID <=26:
            #xticklabels = ["'{}".format(str(year)[-2:]) for year in xticklabels]
            xticklabels = [year for year in xticklabels]
        else:
            xticklabels = ['' for year in xticklabels]
    else:
        xticks      = [0,35,79,129]
        xticklabels = [1971,2006,2050,2099]
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels, rotation=-45, fontsize=6)

    # present vmax
    if regionID == 0: this_is_its_vrange = ' vmax={}\n vmin={}'.format(vmax,vmin)
    else:             this_is_its_vrange = '  {}\n  {}'.format(vmax,vmin)
    #fig.text(ax_pos.x0, ax_pos.y0-0.055, this_is_its_vrange,
    if 17 <= regionID <=26: y_pad = 0.025  #0.05
    else:                   y_pad = 0.010  #0.03
    fig.text(ax_pos.x0, ax_pos.y0-y_pad, this_is_its_vrange,
             fontsize=6, ha='left', va='top')

    # colorbar
    if regionID == 7:  # AMZ (the most left column in the middle row)
        cax = fig.add_axes([ax_pos.x0-0.050, ax_pos.y0-0.2, 0.012, ax_pos.y1-ax_pos.y0+0.4])
        cb = fig.colorbar(im, cax=cax, orientation='vertical')
        cb.set_label('Probability', fontsize=8, labelpad=-15, rotation=90)  # !!!!
        cax.yaxis.set_label_position('left')
        cax.yaxis.set_ticks_position('left')
        cax.spines['right'].set_visible(False)
        cax.spines['left'].set_visible(False)
        cb.set_ticks(bounds)
        cb.set_ticklabels(['vmin', 'vmax'])
        cb.ax.tick_params(labelsize=6,direction='in')


# -----------------------------------------------------------------------------------------------------------
def figure3(ghm, gcm, aSrc, TPCD, threshold, hist_min, hist_max, dict_info):
    # plot of timeseries of aSrc (nYear)    ensemble
    print('\nFiguring #3 or #5...')
    
    scn           = dict_info['scn']
    soc           = dict_info['soc']
    co2           = dict_info['co2']
    regionID      = dict_info['regionID']
    region_name   = dict_info['region_name']
    index         = dict_info['dict_parameters']['index']
    drought_paras = dict_info['dict_parameters']['drought_paras']
    tChunk        = dict_info['dict_parameters']['tChunk']
    mmeType       = dict_info['dict_parameters']['mmeType']
    fig_directory = dict_info['dict_parameters']['fig_directory']

    #fig = plt.figure(num=3, figsize=(10,4), dpi=dpi_figure)
    fig = plt.figure(num=3, figsize=(10,8), dpi=dpi_figure)
    ax = fig.add_subplot(111)
    ax.set_title('{} {}: {}  ({}_{}_{})'.format(regionID,region_name,TPCD,ghm,gcm,mmeType))
    ax.plot(aSrc)
    if hist_min is not None: ax.axhline(y=hist_min, lw=0.3, color='r')
    if hist_max is not None: ax.axhline(y=hist_max, lw=0.3, color='r')
    if threshold is not None: ax.axhline(y=hist_max, lw=0.45, color='b')
    ax.axvline(x=si, lw=0.3, color='g')
    ax.axvline(x=ei, lw=0.3, color='g')
    ax.axvline(x=TPCD-syear, lw=0.3, color='m')
    ax.set_xticks(range(-1,2099-syear+1,10))  # TODO
    ax.set_xticklabels(range(syear-1,2099+1,10))  # TODO
    ax.set_xlim(0,2100-syear+1)
    #ax.set_ylim(0,350)

    figname = '{}.{}.chunk{:02}.plot.{:02}{}.{}_{}.{}_{}_{}.{}'.format(index, drought_paras, tChunk, regionID, region_name, ghm, gcm, scn, soc, co2, mmeType)
    for suffix in suffixes:
        _fig_directory = os.path.join(fig_directory, suffix[1:], 'timeseries_spatialstats_splitted')
        if not os.path.isdir(_fig_directory): os.makedirs(_fig_directory)
        figpath = os.path.join(_fig_directory, figname+suffix)
        plt.savefig(figpath, dpi=dpi_savefig)
        if suffix == suffixes[0]: print('savefig: {}'.format(figpath))
    plt.close(3)


# -----------------------------------------------------------------------------------------------------------
def figure4(ghm, gcm, PDFs, TPCD, dict_info):
    # single Kernel timeseries  (Global & regional)
    print('\nFiguring #4 or #6...')
    scn           = dict_info['scn']
    soc           = dict_info['soc']
    co2           = dict_info['co2']
    regionID      = dict_info['regionID']
    region_name   = dict_info['region_name']
    index         = dict_info['dict_parameters']['index']
    drought_paras = dict_info['dict_parameters']['drought_paras']
    tChunk        = dict_info['dict_parameters']['tChunk']
    mmeType       = dict_info['dict_parameters']['mmeType']
    fig_directory = dict_info['dict_parameters']['fig_directory']

    if READ_VMAX_OVERALL:
        vmax = read_vmax_overall(dict_info)
    else:
        vmax, dummy = gen_vmax_overall(PDFs)
    vmin = PDFs.min()

    fig = plt.figure(num=4, figsize=(4,10), dpi=dpi_figure)
    ax = fig.add_subplot(111)
    ax.set_title('{:02} {}  {}'.format(regionID, region_name, TPCD))
    im = ax.imshow(PDFs[::-1], vmin=vmin, vmax=vmax, cmap=pdf_cmap)  # TODO: is vmax OK with this?  How about vmax??
    plt.colorbar(im, pad=0.01, aspect=35)
    ax.set_xticks(range(-1,2099-syear+1,10))                # TODO
    ax.set_xticklabels(range(1970,2099+1,10), rotation=90)  # TODO
    ax.set_xlim(0,2100-1970)                                # TODO

    filename ='{}.{}.chunk{:03}.gridChart.{:02}{}.{}_{}_{}.{}_{}_{}.{}'.format(index, drought_paras, tChunk, regionID, region_name, ghm, gcm, scn, soc, co2, mmeType)
    for suffix in suffixes:
        _fig_directory = os.path.join(fig_directory, suffix[1:], 'timeseries_shiftingPDFs_splitted')
        if not os.path.isdir(_fig_directory): os.makedirs(_fig_directory)
        figpath = os.path.join(_fig_directory, filename+suffix)
        plt.savefig(figpath, dpi=dpi_savefig)
        if suffix == suffixes[0]: print('savefig: {}'.format(figpath))
    plt.close(4)


# -----------------------------------------------------------------------------------------------------------
def get_df_pop(rcp, dict_unexceptional_global_mask):
    """
    :param rcp:
    :param unexceptional_global_mask: to maskout normal grids. =(True at normal grids.)
    :return df_pop (DataFrame): about population information
    """

    def population(ssp, population_type, stat=None):
        """
        :retur POPs (list)
        """
        filename = 'population_{}soc_0p5deg_annual_2006-2100.nc4'.format(ssp)
        population_path = os.path.join(population_directory, '{}soc'.format(ssp), filename)
        population = Dataset(population_path)['number_of_people'][:]  # The full period for pop is 2006-2099 (nyear3, ny, nx)
        years_pop_data = range(2006, 2099+1)

        if population_type == 'total':
            POPs = []
            for year in years_pop:  # A yearly loop for 2006-2097  
                POPs.append(population[years_pop_data.index(year)].sum())

        elif population_type == 'affected':
            if how_to_estimate_ensemble_TPCD == 'median of member TPCDs':
                POPs = np.zeros((len(years_pop), len(gcms), len(ghms)), 'float32')
                for (ighm, ghm), (igcm, gcm), (iyear, year) in itertools.product(enumerate(ghms), enumerate(gcms), enumerate(years_pop)):  # 2006-2097
                    POPs[iyear,igcm,ighm] = ma.masked_array(population[years_pop_data.index(year)],
                                                            mask=dict_unexceptional_global_mask['{}_{}'.format(ghm, gcm)][years.index(year)]
                                                            ).sum()
                ##print('-------\nmember info on affected population\n-------')
                ##for (ighm, ghm), (igcm, gcm) in itertools.product(enumerate(ghms), enumerate(gcms)):  # 2006-2097
                ##    print(ghm, gcm, POPs[:,igcm,ighm])

                if stat is None:    POPs = np.median(POPs, axis=(1,2)).tolist()
                elif stat == 'Q25': POPs = np.percentile(POPs, 25, axis=(1,2)).tolist()
                elif stat == 'Q75': POPs = np.percentile(POPs, 75, axis=(1,2)).tolist()

            else:
                POPs = []
                for year in years_pop:  # 2006-2097
                    POPs.append(ma.masked_array(population[years_pop.index(year)], mask=dict_unexceptional_global_mask['all_all'][years.index(year)]).sum())

        return POPs

    # pre-process
    ssp = dictSSP[rcp]
    df_pop = pd.DataFrame(columns=years_pop)

    # global total population
    df_pop['{}_total'.format(ssp)] = population(ssp, 'total')

    # the potential population affected (overall)
    # exclude the population in normal grids
    df_pop['{}_affected'.format(ssp)] = population(ssp, 'affected')
    if how_to_estimate_ensemble_TPCD == 'median of member TPCDs':
        for stat in ['Q25', 'Q75']:
            df_pop['{}_affected_{}'.format(ssp, stat)] = population(ssp, 'affected', stat)

    return df_pop


# -----------------------------------------------------------------------------------------------------------
def figure7(src, info, df_pop, dict_parameters):
    print('\nFiguring # 7...')
    # 7: Global map of TPCD

    if TEST:
        src = 2050 + np.random.rand(360, 720)[:293]
        src = src.astype('int')

    rcp           = dict_parameters['scn']
    soc           = dict_parameters['soc']
    co2           = dict_parameters['co2']
    index         = dict_parameters['index']
    smplType      = dict_parameters['smplType']
    mmeType       = dict_parameters['mmeType']
    tChunk        = dict_parameters['tChunk']
    drought_paras = dict_parameters['drought_paras']
    fig_directory = dict_parameters['fig_directory']

    if   rcp == 'rcp26' and smplType == 'Average': fig_id = '(a)'
    elif rcp == 'rcp85' and smplType == 'Average': fig_id = '(b)'
    elif rcp == 'rcp26' and smplType == 'Extreme': fig_id = '(c)'
    elif rcp == 'rcp85' and smplType == 'Extreme': fig_id = '(d)'
    else: fig_id = ''

    fig7 = plt.figure(num=7, figsize=(4,1.75), dpi=dpi_figure)
    plt.subplots_adjust(left=0.02, bottom=0.02, right=0.98, top=0.98)


    ax = fig7.add_subplot(1,1,1, projection=ccrs.Robinson())
    ax.set_extent([-180, 180, -60, 90], ccrs.PlateCarree())
    ax.outline_patch.set_linewidth(0.1)
    ax.add_feature(cfeature.LAND, facecolor='#d4d4d4')
    ax_pos = ax.get_position()

    # TPCD map
    cmap = cm.get_cmap(TPCDmap_cmap)
    vmin, vmax = syear_tpcdmap, eyear_tpcdmap
    #boundaries = [vmin, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, vmax]
    boundaries = [vmin, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100, 2105]  # added extra to avoid white in the colorbar
    norm = colors.BoundaryNorm(boundaries, cmap.N)


    #img_extent = [-180, 180, -60, 90]
    #im = ax.imshow(src, vmin=vmin, vmax=vmax, norm=norm, cmap=cmap, transform=ccrs.PlateCarree())
    #im = ax.imshow(src, origin='upper', extent=img_extent, vmin=vmin, vmax=vmax, norm=norm, cmap=cmap)  #, transform=ccrs.PlateCarree())
    #im = ax.imshow(src)
    #im = ax.contourf(np.arange(-179.75, 179.75+0.5, 0.5), np.arange(-56.25, 89.75+0.5, 0.5), src)
    lons, lats = np.arange(-179.75, 179.75+0.5, 0.5), np.arange(89.75, -56.25-0.5, -0.5)
    im = plt.contourf(lons, lats, src, levels=boundaries, norm=norm, cmap=cmap, transform=ccrs.PlateCarree())
    ax.add_feature(cfeature.COASTLINE, linewidth=0.3, edgecolor='#000000')

    # regional shape file
    shape_feature = ShapelyFeature(Reader(region_shp_path).geometries(), ccrs.PlateCarree(), edgecolor='#808080', linewidth=0.2, facecolor='none')
    ax.add_feature(shape_feature)
    ###if 'rcp' in rcp:  # TODO. to avoid picontrol
    ###    if region_type == 'AR6_regions':
    ###        region_shp_path2 = os.path.join(mapmask_directory,
    ###                                      'IPCC-AR6-WGI-reference-regions-v4_shapefile',
    ###                                      'shapefile_edited',
    ###                                      'earliestTFE_{}'.format(rcp), 
    ###                                      'reference-regions-AR6')
    ###    elif region_type == 'SREX_regions':
    ###        region_shp_path2 = os.path.join(mapmask_directory, 
    ###                                      'SREX.referenceRegions.editted_TFE08regions_{}'.format(rcp), 
    ###                                      'referenceRegions')
    #shape_feature2 = ShapelyFeature(Reader(region_shp_path2).geometries(), ccrs.PlateCarree(), edgecolor='#000000', linewidth=0.5, facecolor='none')
    #ax.add_feature(shape_feature2)

        
    # add a color bar
    cax = fig7.add_axes([ax_pos.x0+0.4, ax_pos.y0+0.08, 0.37, 0.03])
    cb = fig7.colorbar(im, cax=cax, orientation='horizontal')
    cb.ax.set_title('$TFE_{}  $'.format(tChunk), fontsize=5, pad=0.05)
    cb.set_ticks(ticklabels_year)
    cb.set_ticklabels(ticklabels_year)
    cb.ax.tick_params(labelsize=7, direction='in')
    cb.outline.set_visible(False)

    figname = 'TPCD.globmap.{}.{}.{}.{:03}_{}.{}.{}.{}_{}_{}'.format(info, index, drought_paras, tChunk, threshold_type, mmeType, smplType, rcp, soc, co2)
    for suffix in suffixes:
        _fig_directory = os.path.join(fig_directory, suffix[1:], TPCDmap_cmap)
        if not os.path.isdir(_fig_directory): os.makedirs(_fig_directory)
        figpath = os.path.join(_fig_directory, figname+suffix)
        plt.savefig(figpath , dpi=dpi_savefig)
        if suffix == suffixes[0]: print('savefig: {}'.format(figpath))
    plt.close(7)


# -----------------------------------------------------------------------------------------------------------
def figure89(src, info, df_pop, dict_parameters):
    print('\nFiguring # 8,9...')
    # 8: Global map of STD in TPCD
    # 9: Global map of earliest&latest TPCD among members

    if TEST:
        if info == 'Ensmbl' or info == 'Earliest' or info == 'Latest':
            src = 2050 + np.random.rand(360, 720)[:293]
        else:
            src = 20 + np.random.rand(360, 720)[:293]
        src = src.astype('int')

    rcp           = dict_parameters['scn']
    soc           = dict_parameters['soc']
    co2           = dict_parameters['co2']
    index         = dict_parameters['index']
    smplType      = dict_parameters['smplType']
    mmeType       = dict_parameters['mmeType']
    tChunk        = dict_parameters['tChunk']
    drought_paras = dict_parameters['drought_paras']
    fig_directory = dict_parameters['fig_directory']

    xticklabels_pop = [2010, 2030, 2050, 2070, 2090]

    if   rcp == 'rcp26' and smplType == 'Average': fig_id = '(a)'
    elif rcp == 'rcp85' and smplType == 'Average': fig_id = '(b)'
    elif rcp == 'rcp26' and smplType == 'Extreme': fig_id = '(c)'
    elif rcp == 'rcp85' and smplType == 'Extreme': fig_id = '(d)'
    else: fig_id = ''


    fig7 = plt.figure(num=7, figsize=(4,1.75), dpi=dpi_figure)
    plt.subplots_adjust(left=0.02, bottom=0.02, right=0.98, top=0.98)
    ax = fig7.add_subplot(1,1,1, projection=ccrs.Robinson())
    ax.set_extent([-180, 180, -60, 90], ccrs.PlateCarree())
    ax.outline_patch.set_linewidth(0.1)
    ax.add_feature(cfeature.LAND, facecolor='#d4d4d4')
    ax_pos = ax.get_position()

    if info == 'Ensmbl' or info == 'Earliest' or info == 'Latest':  # TPCD map
        cmap = cm.get_cmap(TPCDmap_cmap)
        vmin, vmax = syear_tpcdmap, eyear_tpcdmap
        #boundaries = [vmin, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, vmax]
        boundaries = [vmin, 2020, 2030, 2040, 2050, 2060, 2070, 2080, 2090, 2100, 2105]  # added extra to avoid white in the colorbar
        norm = colors.BoundaryNorm(boundaries, cmap.N)
        #img_extent = [-180, 180, -60, 90]
        #im = ax.imshow(src, vmin=vmin, vmax=vmax, norm=norm, cmap=cmap, transform=ccrs.PlateCarree())
        #im = ax.imshow(src, origin='upper', extent=img_extent, vmin=vmin, vmax=vmax, norm=norm, cmap=cmap)  #, transform=ccrs.PlateCarree())
        #im = ax.imshow(src)
        #im = ax.contourf(np.arange(-179.75, 179.75+0.5, 0.5), np.arange(-56.25, 89.75+0.5, 0.5), src)
        lons, lats = np.arange(-179.75, 179.75+0.5, 0.5), np.arange(89.75, -56.25-0.5, -0.5)
        im = plt.contourf(lons, lats, src, levels=boundaries, norm=norm, cmap=cmap, transform=ccrs.PlateCarree())
        ax.add_feature(cfeature.COASTLINE, linewidth=0.3, edgecolor='#000000')
        shape_feature = ShapelyFeature(Reader(region_shp_path).geometries(), ccrs.PlateCarree(), edgecolor='#808080', linewidth=0.2, facecolor='none')
        ax.add_feature(shape_feature)
    else:  # STD
        cmap = cm.get_cmap('jet')
        vmin, vmax = 5, 40
        boundaries = range(5, 40+1, 5)  # added extra to avoid white in the colorbar
        norm = colors.BoundaryNorm(boundaries, cmap.N)
        #img_extent = [-180, 180, -60, 90]
        #im = ax.imshow(src, vmin=vmin, vmax=vmax, norm=norm, cmap=cmap, transform=ccrs.PlateCarree())
        #im = ax.imshow(src, origin='upper', extent=img_extent, vmin=vmin, vmax=vmax, norm=norm, cmap=cmap)  #, transform=ccrs.PlateCarree())
        #im = ax.imshow(src)
        #im = ax.contourf(np.arange(-179.75, 179.75+0.5, 0.5), np.arange(-56.25, 89.75+0.5, 0.5), src)
        lons, lats = np.arange(-179.75, 179.75+0.5, 0.5), np.arange(89.75, -56.25-0.5, -0.5)
        im = plt.contourf(lons, lats, src, levels=boundaries, norm=norm, cmap=cmap, transform=ccrs.PlateCarree())
        ax.add_feature(cfeature.COASTLINE, linewidth=0.3, edgecolor='#000000')
        shape_feature = ShapelyFeature(Reader(region_shp_path).geometries(), ccrs.PlateCarree(), edgecolor='#808080', linewidth=0.2, facecolor='none')
        ax.add_feature(shape_feature)

    # add a color bar
    cax = fig7.add_axes([ax_pos.x0+0.4, ax_pos.y0+0.085, 0.35, 0.03])
    cb = fig7.colorbar(im, cax=cax, orientation='horizontal')
    cb.ax.set_title('TPCD', fontsize=10, pad=0.05)

    if info == 'Ensmbl' or info == 'Earliest' or  info == 'Latest':
        cb.set_ticks(ticklabels_year)                   # TODO
        cb.set_ticklabels(ticklabels_year)
    else:
        cb.set_ticks([5,10,15,20,25,30,35,40])
        cb.set_ticklabels([5,10,15,20,25,30,35,40])
    cb.ax.tick_params(labelsize=7, direction='in')
    cb.outline.set_visible(False)

    figname = 'TPCD.globmap.{}.{}.{}.{:03}.{}.{}.{}_{}_{}'.format(info, index, drought_paras, tChunk, mmeType, smplType, rcp, soc, co2)
    for suffix in suffixes:
        _fig_directory = os.path.join(fig_directory, suffix[1:], TPCDmap_cmap)
        if not os.path.isdir(_fig_directory): os.makedirs(_fig_directory)
        figpath = os.path.join(_fig_directory, figname+suffix)
        plt.savefig(figpath , dpi=dpi_savefig)
        if suffix == suffixes[0]: print('savefig: {}'.format(figpath))
    plt.close(7)


# -----------------------------------------------------------------------------------------------------------
def figure10(aSrc, dictRegion, RegionIDs, dict_parameters):
    print('\nFiguring #10...')  # member heatmap of TPCD

    scn           = dict_parameters['scn']
    soc           = dict_parameters['soc']
    co2           = dict_parameters['co2']
    index         = dict_parameters['index']
    smplType      = dict_parameters['smplType']
    mmeType       = dict_parameters['mmeType']
    tChunk        = dict_parameters['tChunk']
    drought_paras = dict_parameters['drought_paras']
    fig_directory = dict_parameters['fig_directory']

    fig = plt.figure(num=10, figsize=(12,6), dpi=dpi_figure)
    #plt.suptitle('TPCD in regional {}'.format(smplType.lower()), fontsize=8)
    plt.subplots_adjust(left=0.15, bottom=0.08, right=0.99, top=0.95)
    ax = fig.add_subplot(111)
    im = ax.imshow(aSrc, vmin=ticklabels_year[0], vmax=eyear_tpcdmap, cmap=cm.jet_r)
    cbar = plt.colorbar(im, pad=0.01, aspect=35, ticks=ticklabels_year, orientation='horizontal')
    cbar.ax.tick_params(labelsize=8)

    ax.xaxis.set_major_formatter(ticker.NullFormatter())
    ax.xaxis.set_minor_locator(ticker.FixedLocator(np.arange(n_regions)))
    ax.xaxis.set_minor_formatter(ticker.FixedFormatter(['{} {}'.format(regionID, dictRegion[regionID]) for regionID in RegionIDs]))
    labels = ax.get_xticklabels(minor=True)
    plt.setp(labels, rotation=-90, fontsize=6)

    ax.yaxis.set_major_formatter(ticker.NullFormatter())
    ax.yaxis.set_minor_locator(ticker.FixedLocator(np.arange(len(ghms)*len(gcms)+1)))
    ax.yaxis.set_minor_formatter(ticker.FixedFormatter((['{:<22}'.format('Ensemble median')] + ['{:<9} {:<12}'.format(ghm, gcm) for ghm in ghms for gcm in gcms])))
    labels = ax.get_yticklabels(minor=True)
    plt.setp(labels, fontsize=6)

    ax.tick_params(axis='y',    which='major', direction='out', length=0.2, width=0.5, colors='black')
    ax.tick_params(axis='x',    which='major', length=0)
    ax.tick_params(axis='both', which='minor', length=0)

    ax.axhline(y=0.5, color='black', lw=0.5)
    ax.axvline(x=0.5, color='black', lw=0.5)
    for ighm, ghm in enumerate(ghms[::-1]):
        if ighm != 0: ax.axhline(y=0.5+len(gcms)*ighm, color='black', ls='--', lw=0.3)
        #ax.text(-0.175, 0.1+0.192*ighm, ghm, transform=ax.transAxes, fontsize=7, va='center', ha='center',  rotation='vertical')

    figname = 'memberresults.{}.{}.{:02}.Emsenble{}.{}.{}_{}_{}'.format(index, drought_paras, tChunk, mmeType, smplType, scn, soc, co2)
    for suffix in suffixes:
        _fig_directory = os.path.join(fig_directory, suffix[1:])
        if not os.path.isdir(_fig_directory): os.makedirs(_fig_directory)
        figpath = os.path.join(_fig_directory, figname+suffix)
        plt.savefig(figpath, dpi=dpi_savefig)
        if suffix == suffixes[0]: print('savefig: {}'.format(figpath))
    plt.close(10)


# -----------------------------------------------------------------------------------------------------------
def create_figureobject11(num, plotType, RegionIDs):
    # 11: plot of timeseries of ensemble aSrcs (Global&Regional) on a basemap
    # 12: plot of timeseries of anomaly in ensemble aSrcs (Global&Regional)

    print('Just a moment. Preparing figure object_11or12 (num={})...'.format(num))
    fig = plt.figure(num=num, figsize=(10,5), dpi=dpi_figure)
    ax = fig.add_subplot(1,1,1)
    plt.subplots_adjust(left=0.02, bottom=0.02, right=0.98, top=0.98)
    ax.axis('off')
    ax_pos = ax.get_position()
   
    # add a nice background map :)
    if fig11_proj == 'robin': bm = Basemap(projection='robin',lon_0=0, ax=ax); bm.drawmapboundary(linewidth=0)
    elif fig11_proj == 'cyl': bm = Basemap(projection='cyl',llcrnrlat=-90,llcrnrlon=-180,urcrnrlat=90,urcrnrlon=180,resolution='l',ax=ax)
    else: print('check fig11_proj...'); sys.exit()
    #if fig11_backcolor == 'fillcontinents': bm.fillcontinents(color='sandybrown')
    if fig11_backcolor == 'fillcontinents': bm.fillcontinents(color='gainsboro')
    elif fig11_backcolor == 'shadedrelief': bm.shadedrelief()
    else: print('check fig11_backcolor'); sys.exit()
    if region_type == 'AR6_regions':
        bm.readshapefile(region_shp_path, 'reference-regions-AR6', linewidth=0.2, color='#656565')
    elif region_type == 'SREX_regions':
        bm.readshapefile(region_shp_path, 'referenceRegions', linewidth=0.2, color='#656565')
    elif region_type == 'HydroBASINS_lev1':
        bm.readshapefile(region_shp_path, 'hybas_lake____lev02_v1c_merge', linewidth=0.3, color='k')
    elif region_type == 'hydroregions':
        bm.readshapefile(region_shp_path, 'meybeck_et_al_2013_hydroregions', linewidth=0.3, color='k')

    # add brank axes
    #aspect = 0.93
    aspect = 0.9
    #aspect = 0.85
    #aspect = 0.8
    #aspect = 0.7
    dictAxes = {}
    for regionID in RegionIDs:
        if regionID == 0:
            ix, iy = ax_pos.x0+dictIPCCregion[0][3]  # position of regional box
            dictAxes[0] = fig.add_axes([ix, iy, aspect*0.2, 0.24])
        else:
            ix, iy = ax_pos.x0+dictIPCCregion[regionID][3]  # position of regional box
            dictAxes[regionID] = fig.add_axes([ix, iy, aspect*0.13, 0.13])

    # add a horizontal zero line if this is an anomaly plot
    if plotType == 'anomaly': 
        for regionID in RegionIDs:
            if regionID == 0: dictAxes[0].axhline(y=0, lw=0.5, linestyle='--',color='gray')
            else: dictAxes[regionID].axhline(y=0, lw=0.4, linestyle='--',color='gray')

    return fig, dictAxes

def add_timeseries2axis(plotType, aSRC, fig, ax, dict_info):
    print('\nadd_timeseries2axis...')

    def mann_whitney_u_test(src_h, src_f):
        """
        To check if changes are significant.
        :param src_h: historical time series (full period)
        :param src_f: future time series
        :return: True or False
        """
        n_sample = 30
        weight = np.ones(n_sample)/n_sample
        index_max = np.argmax(np.convolve(src_f, weight, mode='valid'))
        p_value = stats.mannwhitneyu(src_h, src_f[index_max: index_max+n_sample], alternative='two-sided')[1]
        return True if p_value < 0.05 else False

    scn      = dict_info['scn']
    soc      = dict_info['soc']
    co2      = dict_info['co2']
    regionID = dict_info['regionID']
    smplType = dict_info['dict_parameters']['smplType']
    regionID = dict_info['regionID']

    #if kernel_randomvariable_unit == '%':
    #    aSRC = aSRC / 365
    aSrc  = np.array(aSRC[-1])  / 365. * 100  # ensemble timeseries (nYear)
    aSrcs = np.array(aSRC[:-1]) / 365. * 100  # each member (nGHM*nGCM,nYear)

    if plotType == 'anomaly': 
        if not smplType == 'Extreme': 
            aSrc  = aSrc - np.mean(aSrc[si:ei])
            aSrcs = (aSrcs.T - np.mean(aSrcs[:,si:ei], axis=1)).T
        elif smplType == 'Extreme':
            aSrc  = aSrc / np.mean(aSrc[si:ei])
            aSrcs = (aSrcs.T / np.mean(aSrcs[:,si:ei], axis=1)).T

    if   rangeType == 'minmax':
        aMinTS = np.min(aSrcs, axis=0)
        aMaxTS = np.max(aSrcs, axis=0)
    elif rangeType == 'IQR':
        aMinTS = np.percentile(aSrcs, 25., axis=0)
        aMaxTS = np.percentile(aSrcs, 75., axis=0)
    elif rangeType == '2sigma':
        aSigma = np.std(aSrcs, axis=0)
        aMinTS = aSrc - 2*aSigma
        aMaxTS = aSrc + 2*aSigma

    ##if (len(ghms) == 5 and len(gcms) == 4) or (len(ghms) == 1 and len(gcms) == 1):  # for full member. This is for the main process.
    if regionID == 0: lw=0.6  #0.5
    else            : lw=0.5  #0.4
    if scn == 'picontrol': lw = 0.25

    if scn == 'picontrol':
        color_hist = dictSim[scn]
        color_future = dictSim[scn]
        alpha_hist = 0.1
        alpha_future = 0.1
    elif 'rcp' in scn:
        color_hist = '#2b2b2b'
        color_future = dictSim[scn]
        alpha_hist = 0.05
        alpha_future = 0.25
    else: raise ValueError

    print('check')
    for variable in ['scn', 'lw', 'color_hist', 'alpha_hist', 'alpha_future']:
        print('{:<12}: {}')
    sys.exit()


    # uncertainty shade
    ax.fill_between(years_historical_for_plot, 
                    aMinTS[index_hist_s_for_plot:index_hist_e_for_plot+1], 
                    aMaxTS[index_hist_s_for_plot:index_hist_e_for_plot+1], 
                    facecolor=color_hist, alpha=alpha_hist)
    ax.fill_between(years_future_for_plot, 
                    aMinTS[index_hist_e_for_plot:], 
                    aMaxTS[index_hist_e_for_plot:], 
                    facecolor=color_future, alpha=alpha_future)
    src_hist =    aSrc[index_hist_s_for_plot:index_hist_e_for_plot+1]
    src_future =  aSrc[index_hist_e_for_plot:]
    ax.plot(years_historical_for_plot, src_hist,   lw=lw, color=color_hist)
    ax.plot(years_future_for_plot,     src_future, lw=lw, color=color_future, label=scn)
    ##else:
    ##    if len(gcms) == 1 or len(ghms) == 1:  # This is for the check process.
    ##        if len(ghms) == 5 and len(gcms) == 1: members = ghms
    ##        elif len(ghms) == 1 and len(gcms) == 4: members = gcms
    ##        for imember, member in enumerate(members):
    ##            src_hist =    aSrcs[imember, index_hist_s_for_plot:index_hist_e_for_plot+1]
    ##            src_future =  aSrcs[imember, index_hist_e_for_plot:]
    ##            ax.plot(years_historical_for_plot, src_hist,   lw=lw/3, color=dictMember[member])
    ##            if scn == 'rcp26':
    ##                ax.plot(years_future_for_plot,     src_future, lw=lw/3, color=dictMember[member], label=member)
    ##            else:
    ##                ax.plot(years_future_for_plot,     src_future, lw=lw/3, color=dictMember[member])

    #significance = mann_whitney_u_test(src_hist, src_future)  # True or False
    #if significance:
    #    if regionID == 0: 
    #        if scn == 'rcp26': loc_x = 0.465
    #        elif scn == 'rcp85': loc_x = 0.535
    #        ax.text(loc_x, 0.90, '*', color=dictSim[scn], va='top', ha='center', fontsize=10, transform=ax.transAxes)
    #    else: 
    #        if scn == 'rcp26': loc_x = 0.08
    #        elif scn == 'rcp85': loc_x = 0.15
    #        ax.text(loc_x, 0.90, '*', color=dictSim[scn], va='top', ha='center', fontsize=10, transform=ax.transAxes)
    #    print('added timeseries to axis for {} in {}. Change is significant!!!'.format(regionID,scn))
    #else:
    #    print('added timeseries to axis for {} in {}.'.format(regionID,scn))


def custum_timeseries_axis(fig, ax, plotType, regionID):
    print('custum_timeseries_axis... {}'.format(regionID))
    plt.sca(ax)

    years_on_xaxis = [syear, 2006, 2050, eyear]
    ax.set_xlim([syear, eyear])
    ax.set_xticks(years_on_xaxis)
    ax.set_xticklabels(years_on_xaxis)
    if regionID == 0:
        ax.tick_params(axis='x', direction='in', length=2, pad=0.8, labelsize=6)
        ax.tick_params(axis='y', direction='in', length=1, pad=0.4, labelsize=6)
    else:
        ax.tick_params(axis='x', direction='in', length=2, pad=0.1, labelsize=0)
        ax.tick_params(axis='y', direction='in', length=1, pad=0.4, labelsize=4.5)

    if plotType == 'absolute' and ax.get_ylim()[0] < 0: ax.set_ylim(bottom=0)

    spine_color = '#525252'
    lw = 0.45
    ax.patch.set_alpha(0.7)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_color(spine_color)
    ax.spines['bottom'].set_color(spine_color)
    ax.spines['left'].set_linewidth(lw)
    ax.spines['bottom'].set_linewidth(lw)

    if regionID == 0: 
        ax.text(0.50, 0.99, 'Global', va="top", ha="center", fontsize=6, transform=ax.transAxes)
        ax.text(0, 1.1, '[%]', va='bottom', ha='right', fontsize=6, transform=ax.transAxes)
        ax.set_xlabel('Year', fontsize=6)
        ax.set_ylabel('Regional average FDD', fontsize=6)
    else: 
        regionName = dictIPCCregion[regionID][0]
        ax.text(0.05, 0.99, '{} {}'.format(regionID, regionName), va="top", ha="left", fontsize=6, transform=ax.transAxes)

    if regionID == 0: ax.legend(prop={'size' : 6},labelspacing=0.5,loc='upper left', frameon=False)

    return fig, ax


# ----------------------------------------------------------------------------------------------------------------------
def create_figureobject12(_num, RegionIDs):
    # 12: plot of timeseries of ensemble aSrcs (Global&Regional) for each regions

    for regionID in RegionIDs:
        num = int('{}{:02}'.format(_num, regionID))
        print('Just a moment. Preparing figure object_12 (num={})...'.format(num))

        # --- Fig settings
        #figsize=(3, 1.2); left=0.065; bottom=0.03; right=0.98; top=0.95  # for Fig2
        figsize=(2.655, 1.031); left=0.09; bottom=0.15; right=0.97; top=0.85  # for ExtendedDataFig3

        fig = plt.figure(num=num, figsize=figsize, dpi=dpi_figure)
        plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top)
        ax = fig.add_subplot(111)
        #ax.axis('off')

def update_indevidual_timeseries(_num, aSRC, dict_info):

    scn      = dict_info['scn']
    soc      = dict_info['soc']
    co2      = dict_info['co2']
    regionID = dict_info['regionID']

    #if (len(ghms) == 5 and len(gcms) == 4) or (len(ghms) == 1 and len(gcms) == 1):  # for full member. This is for the main process.

    #if regionID == 0: lw = 0.5  #0.9  #
    #else            : lw = 0.4  #0.8  #
    ##if scn == 'picontrol': lw = 0.3
    #if scn == 'rcp85': lw = 0.4
    lw = 0.5

    linestyle = '-'
    if co2 == '2005soc': linestyle = '--'

    if scn == 'picontrol':
        color_hist = dictSim[scn]
        color_future = dictSim[scn]
        alpha_hist = 0.15
        alpha_future = 0.15
    elif 'rcp' in scn:
        color_hist = '#2b2b2b'
        color_future = dictSim[scn]
        alpha_hist = 0.05
        alpha_future = 0.25
    else: raise ValueError

    aSrc  = np.array(aSRC[-1])  / 365. * 100  # ensemble timeseries (nYear)
    aSrcs = np.array(aSRC[:-1]) / 365. * 100  # each member (nGHM*nGCM,nYear)

    if   rangeType == 'minmax':
        aMinTS = np.min(aSrcs, axis=0)
        aMaxTS = np.max(aSrcs, axis=0)
    elif rangeType == 'IQR':
        aMinTS = np.percentile(aSrcs, 25., axis=0)
        aMaxTS = np.percentile(aSrcs, 75., axis=0)
    elif rangeType == '2sigma':
        aSigma = np.std(aSrcs, axis=0)
        aMinTS = aSrc - 2*aSigma
        aMaxTS = aSrc + 2*aSigma

    num = int('{}{:02}'.format(_num, regionID))
    fig = plt.figure(num)
    ax = fig.axes[0]
    #plt.sca(ax)
    # uncertainty shade
    ax.fill_between(years_historical_for_plot,
                    aMinTS[index_hist_s_for_plot:index_hist_e_for_plot+1],
                    aMaxTS[index_hist_s_for_plot:index_hist_e_for_plot+1],
                    facecolor=color_hist, 
                    alpha=alpha_hist)
    ax.fill_between(years_future_for_plot,
                    aMinTS[index_hist_e_for_plot:],
                    aMaxTS[index_hist_e_for_plot:],
                    facecolor=color_future, 
                    alpha=alpha_future)
    src_hist   = aSrc[index_hist_s_for_plot:index_hist_e_for_plot+1]
    src_future = aSrc[index_hist_e_for_plot:]
    ax.plot(years_historical_for_plot, src_hist,   linestyle=linestyle, lw=lw, color=color_hist)
    ax.plot(years_future_for_plot,     src_future, linestyle=linestyle, lw=lw, color=color_future, label=scn.upper())

def custum_indevidual_timeseries_axis(_num, regionID, region_name):

        num = int('{}{:02}'.format(_num, regionID))
        years_on_xaxis = [syear, 2006, 2050, eyear]
        spine_color = '#525252'
        lw = 0.45

        fig = plt.figure(num)
        ax = fig.axes[0]
        ax.yaxis.get_major_locator().set_params(integer=True)
        ax.set_xlim([syear, eyear])
        ax.set_xticks(years_on_xaxis)
        #if regionID == 3:  # WNA
        #    ax.set_xticklabels(years_on_xaxis, ha='center')  #, rotation=20)
        #else:
        ax.set_xticklabels(['']*len(years_on_xaxis))
        ax.tick_params(axis='x', direction='in', length=3.5, pad=1, labelsize=7.5)
        #ax.tick_params(axis='x', direction='in', length=3.5, pad=1, labelsize=9)
        ax.tick_params(axis='y', direction='in', length=3,   pad=1, labelsize=9)
        ax.patch.set_alpha(0.7)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_visible(False)
        ax.spines['left'].set_color(spine_color)
        #ax.spines['bottom'].set_color(spine_color)
        ax.spines['left'].set_linewidth(lw)
        #ax.spines['bottom'].set_linewidth(lw)
        #if regionID in [9,10,15,26]: _y = -0.065; va = 'top'
        #else:                        _y =  1.015; va = 'bottom'
        _y = 1; va = 'top'
        _x = 0.05; ha = 'left' 
        #ax.text(0.05, 0.99, '{} {}'.format(regionID, region_name), va="top", ha="left", fontsize=6, transform=ax.transAxes)
        #ax.text(0.07, 1.03, '{} {}'.format(regionID, region_name), ha='left', va='bottom', fontsize=9, transform=ax.transAxes)
        #ax.text(_x, _y, '{} {}'.format(regionID, region_name), ha=ha, va=va, fontsize=9, transform=ax.transAxes)
        #ax.text(0, 1.05, '[%]', ha='right', va='bottom', fontsize=8, transform=ax.transAxes)
        ax.text(_x, _y, '{} {}'.format(regionID, region_name), ha=ha, va=va, fontsize=8, transform=ax.transAxes)
        #ax.set_xlabel('Year', fontsize=10)
        #ax.set_ylabel('Regional average FDD', fontsize=10)
        #if regionID == 3: 
        #    ax.legend(prop={'size' : 9}, 
        #              bbox_to_anchor=(0, 1.08), loc='upper left', 
        #              labelspacing=0.2, 
        #              frameon=False)


def close_figureobject12(_num, dict_parameters, regionID, region_name):

    index         = dict_parameters['index']
    drought_paras = dict_parameters['drought_paras']
    tChunk        = dict_parameters['tChunk']
    mmeType       = dict_parameters['mmeType']
    smplType      = dict_parameters['smplType']
    fig_directory = dict_parameters['fig_directory']
    # figure names...
    figname = 'ts.abs.{}.{}.{}.{:03}.Emsenble_{}.{}.{:02}.{}'.format(rangeType, index,
                                                                     drought_paras,
                                                                     tChunk, mmeType, smplType,
                                                                     regionID, region_name)
    num = int('{}{:02}'.format(_num, regionID))
    print('\n================\nclose_figureobject ({})'.format(num))
    plt.figure(num)
    for suffix in suffixes:
        _fig_directory = os.path.join(fig_directory, suffix[1:], 'timeseries_indevidual')
        if not os.path.isdir(_fig_directory): os.makedirs(_fig_directory)
        figpath = os.path.join(_fig_directory, figname+suffix)
        plt.savefig(figpath, dpi=dpi_savefig)
        if suffix == suffixes[0]: print('savefig: {}'.format(figpath))
    plt.close(num)
    print('\n')


# -----------------------------------------------------------------------------------------------------------
def figure13(dict_pops, dict_parameters):
    # population plot

    index           = dict_parameters['index']
    smplType        = dict_parameters['smplType']
    mmeType         = dict_parameters['mmeType']
    tChunk          = dict_parameters['tChunk']
    drought_paras   = dict_parameters['drought_paras']
    fig_directory   = dict_parameters['fig_directory']
    xticklabels_pop = [2010, 2030, 2050, 2070, 2090]

    fig7 = plt.figure(num=13, figsize=(3.5,4), dpi=dpi_figure)
    plt.subplots_adjust(left=0.15, bottom=0.1, right=0.97, top=0.9)
    ax = fig7.add_subplot(111)
    #ax.text(0, 1, '(b)', va='bottom', ha='right', fontsize=10, transform=ax.transAxes)

    for scn in scns:
        ssp = dictSSP[scn]
        ax.plot(years_pop, dict_pops[scn]['{}_total'.format(ssp)].values, 
                color=dictSim[scn], linestyle='--', lw=1)  #, label='global total')          # gray
        if how_to_estimate_ensemble_TPCD == 'median of member TPCDs':
            ax.fill_between(years_pop,
                            dict_pops[scn]['{}_affected_Q25'.format(ssp)].values, 
                            dict_pops[scn]['{}_affected_Q75'.format(ssp)].values, 
                            color=dictSim[scn], alpha=0.1)
        ax.plot(years_pop, dict_pops[scn]['{}_affected'.format(ssp)].values, 
                color=dictSim[scn], linestyle='-', lw=1, label='{}_{}'.format(scn, ssp))  # orange

    ax.set_xticks(xticklabels_pop)
    ax.set_xticklabels(xticklabels_pop)  #, rotation=20)
    #ax.set_xlim([years_pop[0], years_pop[-1]])
    ax.set_xlim([years_analysis_future[0], years_pop[-1]])
    ax.set_xlabel('Year', fontsize=10)
    #ax.set_ylim(ymin=0, ymax=9e+9)
    ax.set_ylim(ymin=0, ymax=10e+9)
    #ax.set_ylim(ymax=9e+9)
    ax.tick_params(axis='both', labelsize=10, direction='in')#, pad=0.20)
    plt.legend(prop={'size': 10}, labelspacing=0.2, loc='upper left', frameon=False)
    ax.set_ylabel('Population [cap]', fontsize=10)
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=True))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

    figname = 'PopUnderUD.{}.{}.{:03}.{}.{}'.format(index, drought_paras, tChunk, mmeType, smplType)
    for suffix in suffixes:
        _fig_directory = os.path.join(fig_directory, suffix[1:])
        if not os.path.isdir(_fig_directory): os.makedirs(_fig_directory)
        figpath = os.path.join(_fig_directory, figname+suffix)
        plt.savefig(figpath , dpi=dpi_savefig)
        if suffix == suffixes[0]: print('savefig: {}'.format(figpath))
    plt.close(13)


# -----------------------------------------------------------------------------------------------------------
def close_figureobject(num, dict_parameters):
    print('\n================\nclose_figureobject ({})'.format(num))

    scn           = dict_parameters['scn']
    soc           = dict_parameters['soc']
    co2           = dict_parameters['co2']
    index         = dict_parameters['index']
    drought_paras = dict_parameters['drought_paras']
    tChunk        = dict_parameters['tChunk']
    mmeType       = dict_parameters['mmeType']
    smplType      = dict_parameters['smplType']
    fig_directory = dict_parameters['fig_directory']

    # figure names...
    if num == 1:
        figname = 'RegionalKernelMap.{}.{}.{:03}.Emsenble_{}.{}.{}_{}_{}.{}'.format(index, drought_paras, tChunk,
                                                                                    mmeType, smplType, scn, soc, co2,
                                                                                    fig1_proj)
    elif num == 2:
        figname = 'RegionalKernelList.{}.{}.{:03}.Emsenble_{}.{}.{}_{}_{}.{}'.format(index, drought_paras, tChunk,
                                                                                     mmeType, smplType, scn, soc, co2,
                                                                                     fig2_proj)
    elif num == 11:
        figname = 'timeseriesBasemap.absolute.{}.{}.{}.{:03}.Emsenble_{}.{}.{}'.format(rangeType, index, drought_paras,
                                                                                       tChunk, mmeType, smplType,
                                                                                       fig11_proj)
    elif num == 12:
        figname = 'timeseriesBasemap.anomaly.{}.{}.{}.{:03}.Emsenble_{}.{}.{}'.format(rangeType, index, drought_paras,
                                                                                      tChunk, mmeType, smplType,
                                                                                      fig11_proj)
    plt.figure(num)
    for suffix in suffixes:
        _fig_directory = os.path.join(fig_directory, suffix[1:])
        if not os.path.isdir(_fig_directory): os.makedirs(_fig_directory)
        figpath = os.path.join(_fig_directory, figname+suffix)
        plt.savefig(figpath, dpi=dpi_savefig)
        if suffix == suffixes[0]: print('savefig: {}'.format(figpath))
    plt.close(num)
    print('\n')


# -----------------------------------------------------------------------------------------------------------
def pdf2nc(PDFs, dict_info):
    #print('pdf2nc... PDFs {}'.format(PDFs.shape))

    regionID        = dict_info['regionID']
    region_name     = dict_info['region_name']
    ghm             = dict_info['ghm']
    gcm             = dict_info['gcm']
    drought_paras   = dict_info['dict_parameters']['drought_paras']
    ncout_directory = dict_info['dict_parameters']['ncout_directory']

    _nday, _nyear = PDFs.shape
    PDFs = PDFs.T

    filename = 'drought_pdfs_{}.{}.{}.{:02}.{}.nc4'.format(drought_paras, ghm, gcm, regionID, region_name)
    outpath = os.path.join(ncout_directory, filename)
    
    rootgrp = Dataset(outpath, 'w', format='NETCDF4')
    import time
    rootgrp.title = 'Drought PDF at {}.{}'.format(regionID, region_name)
    rootgrp.description = 'Regional Drougnt PDF for {:02}{} ({}-{}): {} faced by {}, {}'.format(regionID, region_name, syear, eyear, ghm, gcm, drought_paras)
    rootgrp.history = 'Created ' + time.ctime(time.time())
    rootgrp.source = 'ISI-MIP FT: {}_{}'.format(ghm, gcm)
    rootgrp.memo = 'These druoght PDFs derive from the average number of drought days (NDD)'
    rootgrp.memo2 = ncout_directory
    rootgrp.institution = 'NIES'
    rootgrp.contact = 'satoh.yusuke@nies.go.jp'

    _year = rootgrp.createDimension('year', size=_nyear)
    _day  = rootgrp.createDimension('day', size=_nday)

    _years = rootgrp.createVariable('year', 'int', ('year',))
    _years.axis = 'Y'
    _years.units = 'year'
    _days = rootgrp.createVariable('day', 'int', ('day',))
    _days.axis = 'X'
    srcs = rootgrp.createVariable('pdfs', 'f8', ('year', 'day'), fill_value=1e+20, chunksizes=(_nyear, _nday), zlib=True)
    srcs.long_name = 'drought pdf of {}.{}'.format(regionID, region_name)
    srcs.unit = '(probability)'  # unitName

    #=== Allocate data ===
    _years[:] = years  # range(syear, eyear+1)
    _days[:] = range(_nday)
    srcs[:] = PDFs
    #=== close ===
    rootgrp.close()
    

# -----------------------------------------------------------------------------------------------------------
def writeoutinbinary(aSrc, name, fig_directory):
    fig_directory = os.path.join(fig_directory, 'bin')
    if not os.path.isdir(fig_directory): os.makedirs(fig_directory)
    outpath = os.path.join(fig_directory, 'TPCD.result.{}.bin'.format(name))
    if   isinstance(aSrc,np.ma.MaskedArray): aSrc.filled(1e+20).tofile(outpath)
    elif isinstance(aSrc,np.ndarray): aSrc.tofile(outpath)


# -----------------------------------------------------------------------------------------------------------
def main(*args):
    print('\n\nStart Job !!! d(^.^)\n')
    strTime = datetime.datetime.now()

    #=== preparations
    #RegionMap, dictRegion = get_regionInfo(region_type)
    RegionMap, dictRegion = get_region_info(region_type, TEST=TEST)
    RegionIDs = list(dictRegion.keys())
    RegionIDs.sort()             # region ID for Global, 0, is not included here.
    if GLOBAL:
        dictRegion[0] = 'GLB'
        RegionIDs = [0] + RegionIDs
    regions = ['{}.{}'.format(regionID, dictRegion[regionID]) for regionID in RegionIDs]

    # main processes...
    #for win, Len, tau, Q, soc, season, index, mmeType, tChunk in itertools.product(wins, LENs, TAUs, Qs, socs, seasons, IDXs, mmeTYPs, tChunkType):
    for win, Len, tau, Q, season, index, mmeType, tChunk in itertools.product(wins, LENs, TAUs, Qs, seasons, IDXs, mmeTYPs, tChunkType):
        dict_stats = {}
        for smplType in smplTYPs:

            if index != 'nDayTot': print('Ah, index is not nDayTot...'); sys.exit()  # tentative !!!!!!!!!!!!
            if tChunk == 'full': tChunk = nyear
            if season == 'ALL': _season = 'annual'
            else: _season = season
            drought_paras = 'Q{:02}win{:02}Len{:02}tau{}{}'.format(Q, win, Len, tau, _season)

            #if autQvalue: 
            #    qvalType = 'Qvalue_historical_{}'.format(soc)
            #else:
            #    if soc == 'nosoc': qvalType = 'Qvalue_historical_nosoc'
            #    else:              qvalType = 'Qvalue_historical_histsoc'
            qvalType = 'Qvalue_historical_histsoc_co2'  # use consisntent reference throughout this analysis.

            member_subdirectory = 'all_all'
            #if   len(ghms) == 5 and len(gcms) == 4: member_subdirectory = 'all_all' 
            #elif len(ghms) == 5 and len(gcms) == 1: member_subdirectory = 'all_{}'.format(gcms[0])
            #elif len(ghms) == 1 and len(gcms) == 4: member_subdirectory = '{}_all'.format(ghms[0])
            #elif len(ghms) == 1 and len(gcms) == 1: member_subdirectory = '{}_{}'.format(ghms[0], gcms[0])
            #else: print('hmm... How do you want to deal with member_subdirectory??'); sys.exit()

            # ---------- figure directory common ----------
            fig_directory_common = os.path.join(
                fig_directory_main, region_type, member_subdirectory,
                drought_paras,
                '{}.{}.{:02}.{}.{}-{}.{}.{:02}.{}_{}-{}.{}-{}.plt{}-'.format(
                    experiment,
                    value_type, window, smplType[:3],
                    syear_analysis_historical, eyear_analysis_historical,  # The TPCD (historical) reference period
                    threshold_type[:3], tChunk, #soc, co2,
                    qvalType[7:], ref_syear, ref_eyear,
                    syear, eyear,                                          # The full period of the analysis
                    syear_historical_for_plot   # The period for time series plot(s)
                    ),
                )
            # ---------------------------------------------

            # a directory for figures
            fig_directory = os.path.join(fig_directory_common, today)
            if not os.path.isdir(fig_directory): os.makedirs(fig_directory)

            if which_hist_percentile == 'reuse' or which_member_timeseries == 'reuse' or READ_VMAX_OVERALL:
                ex_excel_directory = os.path.join(fig_directory_common, which_date_excel)
                if not os.path.isdir(ex_excel_directory): print('ex_excel_directory: {} does not exist!!'.format(ex_excel_directory)); sys.exit()

            # put Logs of reuse options
            for reuse_option in ['how_to_gen_ensemble_pdf', 'how_to_gen_member_pdf',
                                 'how_to_plot_ensemble_timeseries', 'how_to_estimate_ensemble_TPCD',
                                 'which_hist_percentile', 'which_member_timeseries']:
                if   reuse_option is 'how_to_gen_member_pdf'   and how_to_gen_member_pdf   == 'reuse': which_date = which_date_pdf
                elif reuse_option is 'which_hist_percentile'   and which_hist_percentile   == 'reuse': which_date = which_date_excel
                elif reuse_option is 'which_member_timeseries' and which_member_timeseries == 'reuse': which_date = which_date_excel
                else: which_date = ''
                empty_file_name = '{}-----{}{}'.format(reuse_option, str(eval(reuse_option)).replace(' ', '_'), which_date)
                empty_file_path = os.path.join(fig_directory, empty_file_name)
                with open(empty_file_path, 'w') as f: pass

            # create figure object for figure11 and figure12  (time-series plots)
            if Figure11: fig11, dictAxes11 = create_figureobject11(11, 'absolute', RegionIDs)
            if Figure12: create_figureobject12(12, RegionIDs)
            if Figure13: fig12, dictAxes12 = create_figureobject11(13, 'anomaly',  RegionIDs)

            # reading data
            dict_parameters = {}
            for parameter in ['Q', 'win', 'Len', 'tau', 'season', 'index', 'qvalType']:  #, 'soc']:
                dict_parameters[parameter] = eval(parameter)
            if how_to_gen_member_pdf == 'reuse' and how_to_gen_ensemble_pdf == 'reuse':
                print('skip loading nc files...')
            else:
                print('\nLoading data...')
                SRC = array([[[[[load_nc(co2, soc, scn, gcm, ghm, dict_parameters) for ghm in ghms] for gcm in gcms] for scn in scns] for soc in socs] for co2 in co2s])
                SRC = ma.masked_equal(SRC, 1e+20)  # orignally, missing_value in each data is 1e+20
                SRC = ma.masked_array(SRC, mask=np.resize(landsea_mask, SRC.shape))  #(nCO2, nSOC, nSCN, nGCM, nGHM, nYear, nY, nX)
                SRC = SRC.filled(1e+20)
                SRC = check_nan(SRC, 2)
                print('SRC.shape: {}'.format(SRC.shape))  # nYear = 2099 - 1861 + 1 = 239   (all period of input data)

                # conv. absolute to anomaly
                if value_type == 'anm':
                    SRC = absolute_to_anomaly(SRC)

            dict_pops = {rcp: None for rcp in scns}
            if eyear_analysis_historical > 2005: dict_threshold_rcp85 = {}

            # ------------- #
            # scenario loop #
            # ------------- #
            for (iscn, scn), (isoc, soc), (ico2, co2) in itertools.product(enumerate(scns), enumerate(socs), enumerate(co2s)):  # loop for rcp85 & rcp26

                # a directory to store PDFs in netcdf
                ncout_directory = os.path.join(fig_directory_common, today, 'nc', scn)
                if not os.path.isdir(ncout_directory): os.makedirs(ncout_directory)
                if how_to_gen_member_pdf == 'reuse':
                    ex_nc_directory = os.path.join(fig_directory_common, which_date_pdf, 'nc', scn)
                    if not os.path.isdir(ex_nc_directory):
                        print('ex_nc_directory: {} does not exist!!'.format(ex_nc_directory)); sys.exit()

                parameters = ['index', 'smplType', 'tChunk', 'mmeType', 'win', 'Len', 'tau', 'Q',
                              'scn', 'soc', 'co2', 
                              'season', 'drought_paras', 'qvalType', 'fig_directory']
                if which_hist_percentile == 'reuse' or which_member_timeseries == 'reuse' or READ_VMAX_OVERALL:
                    parameters.append('ex_excel_directory')
                if how_to_gen_member_pdf == 'reuse':
                    parameters.append('ex_nc_directory')
                dict_parameters = {}
                for parameter in parameters+['ncout_directory']:
                    dict_parameters[parameter] = eval(parameter)
                print('\n\n--------------------------------------------')
                for k in parameters: print('{:<8}: {}'.format(k, dict_parameters[k]))
                print('Regions: {}\n--------------------------------------------\n'.format(RegionIDs))

                df_TPCD = create_df_TPCD(dictRegion, RegionIDs)  # TPCD table
                df_stats = create_df_stats(dictRegion, RegionIDs)  # data table
                EnsmOUT = np.full((360,720), 1e+20)                                            # To map TPCD (TPCD)
                mmbrOUT = np.full((len(gcms),len(ghms),360,720), 1e+20)                        # To map TPCD (TPCD)
                mmbrChart = np.full((len(gcms)*len(ghms)+1, len(dictRegion.keys())+1), 1e+20)  # compare each TPCD


                if Figure1:  # PDF timeseries on basemap
                    fig1, ax01_pos = create_figureobject1(1)
                if Figure2:  # PDF timeseries in the list style
                    fig2 = create_figureobject2(2)
                    if READ_VMAX_OVERALL:
                        df_vrange = None  # read here??
                    else:
                        if VRANGE_DETAIL_TO_EXCEL: vrange_index = ['vmax', 'vmin', 'n_type1', 'n_type2']+list(range(syear, eyear+1))
                        else:                      vrange_index = ['vmax', 'vmin']
                        df_vrange = pd.DataFrame(index=vrange_index, columns=regions)

                df_unprecedented_flag_member = pd.DataFrame(index=['{}_{}_{}'.format(dictRegion[regionID], ghm, gcm) for regionID in RegionIDs for ghm in ghms for gcm in gcms],
                                                            columns=['region', 'ghm', 'gcm']+list(years))
                df_unprecedented_info_member = pd.DataFrame(index=['{}_{}_{}'.format(dictRegion[regionID], ghm, gcm) for regionID in RegionIDs for ghm in ghms for gcm in gcms],
                                                            columns=['region', 'ghm', 'gcm', 'n_event', 'n_year'])
                df_unprecedented_info_median = pd.DataFrame(index=['{}'.format(dictRegion[regionID]) for regionID in RegionIDs], 
                                                            columns=['n_event', 'n_year'])

                # mask a mask to estimate the potential population affected.  start from all True (= masked out)
                dict_unexceptional_global_mask = {}
                for ghm, gcm in itertools.product(ghms, gcms):
                    dict_unexceptional_global_mask['{}_{}'.format(ghm, gcm)] = np.full((nyear, 360, 720), True)
                dict_unexceptional_global_mask['all_all'] = np.full((nyear, 360, 720), True)
                #unexceptional_global_mask = np.full((nyear, 360, 720), True)

                # ----------- #
                # region loop #
                # ----------- #
                if GLOBAL: start = 0
                else:      start = 1
                for iregion, regionID in enumerate(RegionIDs, start=start):
                    region_name = dictRegion[regionID]

                    # Make region mask
                    if region_name == 'GLB': region_mask = landsea_mask
                    else: region_mask = ma.make_mask(RegionMap != regionID)           #(nY,nX)
                    YY,XX = where(region_mask == 0)   # Grids with non-masked     # To fill in EnsmOUT and map TPCD (TPCD)
                    print('\n=================\n{:03} {} ({})\n================='.format(regionID,region_name, len(YY)))

                    n_regionallandgrids = where(np.ma.mask_or(landsea_mask, region_mask) == 0)[0].shape[0]
                    print('n_regionallandgrids = {}'.format(n_regionallandgrids))
                    if n_regionallandgrids == 0:  # count the number of land grid in the region
                        print('Caution!!!! This area has no land grids !?!?!? Hence, we skip it for now...\n\n')
                        pass
                    else:
                        timeserieses = []  # a list of time serieses for plots
                        # ------------------------------ #
                        #   processes for each  GCM/GHM  #
                        # ------------------------------ #
                        for (ighm, ghm), (igcm, gcm) in itertools.product(enumerate(ghms), enumerate(gcms)):
                            print('\n\nprocessing  {:03} {} with  {} - {}'.format(regionID, region_name, ghm, gcm))
                            dict_info = {'regionID': regionID, 'region_name': region_name, 
                                         'scn': scn, 'soc': soc, 'co2': co2,
                                         'ghm': ghm, 'gcm': gcm, 
                                         'dict_parameters': dict_parameters,
                                         'df_TPCD': df_TPCD, 'df_stats': df_stats,
                                         'fig_directory_common': fig_directory_common}
                            # pre-process
                            if LOUD: print('01: pre-process')
                            if how_to_gen_member_pdf == 'newly create' or index in ['Mean', 'Median', 'Mode']:
                                aSrc = copy(SRC[ico2, isoc, iscn, igcm, ighm])                                       # (nYear,nY,nX)
                                aSrc = ma.masked_array(aSrc, mask=np.resize(region_mask, aSrc.shape)).filled(1e+20)  # (nYear,nY,nX)
                                aSrc = check_nan(aSrc, 3)
                                print('SRC.shape {} >> aSrc.shape: {}'.format(SRC.shape, aSrc.shape))
                                # resample aSrc for the period of TPCD analysis with moving window sampling
                                if kdfSampling == 'window':
                                    if LOUD: print('01-2: pre-process')
                                    aSrc = window_sampling(aSrc)  # (_nYear,???)  _nYear = nYear - (window-1)
                                    aSrc = check_nan(aSrc, 4)
                                    print('window sampe {}  aSrc.shape: {}'.format(window, aSrc.shape))
                                aSrc, gridarea = get_sample(aSrc, gridArea)
                            else:
                                print('None & None of dummy')
                                aSrc, gridarea = None, None  # dummy
                            if LOUD: print('aSrc: {}-{} ({})'.format(aSrc.min(), aSrc.max(), aSrc.shape))
                            aSrc = check_nan(aSrc, 5)

                            # calulate stats from samples...
                            if LOUD: print('02: calc_yearly_spatialstats')
                            timeseries, PDFs, hist_pctl_day = calc_yearly_spatialstats(aSrc, smplType, gridarea, dict_info)  #(_nYear or nYear)
                            timeseries = check_nan(timeseries, 6)
                            del aSrc
                            # the TPCD analysis
                            if eyear_analysis_historical > 2005 and scn == 'rcp26': ref_thresholds = dict_threshold_rcp85
                            else:                                                   ref_thresholds = None
                            if LOUD: print('03: find_TPCD')
                            TPCD, unprecedented_mask, unexceptional_mask, threshold, hist_min, hist_max, n_event, n_year = find_TPCD(timeseries, tChunk, dict_info, ref_thresholds)
                            print('{} {} {:03} {:02} {} {}'.format(ghm, gcm, tChunk, regionID, region_name, TPCD))
                            if eyear_analysis_historical > 2005 and scn == 'rcp85': 
                                dict_threshold_rcp85['{}_{}_{:03}.{}'.format(ghm, gcm, regionID, region_name)] = threshold

                            # record ----------
                            timeserieses.append(timeseries)  # aSrc (nyear)
                            if PDFs2NC_MEMBER: pdf2nc(PDFs, dict_info)

                            update_dict_stats(dict_stats, scn, soc, co2, region_name, ghm, gcm)
                            if smplType == 'Average':
                                allocate_item_2_dict_stats(dict_stats, 'PDFs', PDFs, dict_info)
                                allocate_item_2_dict_stats(dict_stats, 'spatial_average', timeseries, dict_info)
                            elif smplType == 'Extreme':
                                allocate_item_2_dict_stats(dict_stats, 'historical_percentile', hist_pctl_day, dict_info)

                            df_TPCD['{}.{}'.format(regionID,region_name)]['{}_{}'.format(ghm,gcm)] = TPCD
                            add2df_stats(df_stats, ghm, gcm, '{:02}.{}'.format(regionID, region_name), TPCD, threshold, hist_min, hist_max, hist_pctl_day, timeseries)
                            mmbrChart[ighm*len(gcms)+igcm+1, iregion] = TPCD

                            add2df_unprecedented_flag_member(df_unprecedented_flag_member, region_name, ghm, gcm, unprecedented_mask)
                            add2df_unprecedented_info_member(df_unprecedented_info_member, region_name, ghm, gcm, n_event, n_year)

                            if not region_name == 'GLB':
                                mmbrOUT[igcm, ighm, YY, XX] = TPCD
                                # update dict_unexceptional_global_mask
                                for iy, ix in zip(YY,XX):
                                    dict_unexceptional_global_mask['{}_{}'.format(ghm, gcm)][:,iy,ix] = unexceptional_mask


                            if Figure5: figure3(ghm, gcm, timeseries, TPCD, threshold, hist_min, hist_max, dict_info)  # plot of timeseries of aSrc (nYear)  indevidual member
                            if Figure6: figure4(ghm, gcm, PDFs, TPCD, dict_info)                            # PDF time series  Regional&member  indevidual member
                            del timeseries, gridarea, PDFs
                            if GC: gc.collect()
                            if MEMORY_USAGE_NOISY: 
                                print("\n\n====== in main\n{}{: >45}{}{: >10}{}".format('|','Variable Name','|','  Size','|'))
                                print(" -------------------------- ")
                                for k, v in locals().items():
                                    if hasattr(v, 'size') and not k.startswith('_') and not isinstance(v,types.ModuleType):
                                        print("{}{: >45}{}{: >10}{}".format('|',k,'|',str(v.size),'|'))
                                    elif hasattr(v, '__len__') and not k.startswith('_') and not isinstance(v,types.ModuleType):
                                        print("{}{: >45}{}{: >10}{}".format('|',k,'|',str(len(v)),'|'))
                            print('\ncheck >> memory usage: {}'.format(psutil.virtual_memory().percent))

                        # ------- #
                        # Overall #
                        # ------- #
                        print('\n\nprocessing {:03} {} with all members...'.format(regionID, region_name))
                        dict_info = {'regionID': regionID, 'region_name': region_name, 
                                     'scn': scn, 'soc': soc, 'co2': co2,
                                     'ghm': 'all', 'gcm': 'all', 
                                     'dict_parameters': dict_parameters, 
                                     'df_TPCD': df_TPCD, 'df_stats': df_stats,
                                     'fig_directory_common': fig_directory_common}
                        # pre-process
                        if LOUD: print('01: pre-process')
                        #if how_to_gen_member_pdf == 'extract from a pdf':  bug??
                        if how_to_gen_member_pdf == 'newly create':
                            if LOUD: print('01-1: copy')
                            aSrc = copy(SRC[ico2, isoc, iscn])
                            if LOUD: print('01-2: mask')
                            aSrc = ma.masked_array(aSrc, mask=np.resize(region_mask,aSrc.shape)).filled(1e+20)  #(nGCM,nGHM,nYear,nY,nX)
                            if LOUD: print('01-3: transpose')
                            aSrc = aSrc.transpose(2,0,1,3,4)                                                    #(nYear,nGCM,nGHM,nY,nX)
                            # resample aSrc for the period of TPCD analysis with moving window sampling
                            if kdfSampling == 'window':
                                aSrc = window_sampling(aSrc)                                                    #(_nYear, ???)
                            aSrc, gridarea = get_sample(aSrc, gridArea)
                        else:
                            print('None & None of dummy. oh...')
                            aSrc, gridarea = None, None  # dummy

                        # calculate stats from samples...
                        if LOUD: print('02: calc_yearly_spatialstats')
                        timeseries, PDFs, hist_pctl_day = calc_yearly_spatialstats(aSrc, smplType, gridarea, dict_info)  #(nYear or _nYear)
                        del aSrc
                        # the TPCD analysis overall
                        # how_to_estimate_ensemble_TPCD: 'median of member TPCDs' or else
                        if eyear_analysis_historical > 2005 and scn == 'rcp26': ref_thresholds = dict_threshold_rcp85
                        else:                                                   ref_thresholds = None
                        if LOUD: print('03: find_TPCD')
                        TPCD, unprecedented_mask, unexceptional_mask, threshold, hist_min, hist_max, n_event, n_year = find_TPCD(timeseries, tChunk, dict_info, ref_thresholds)
                        print('Ensemble {:03} {:02} {} {} {}'.format(tChunk, regionID, region_name, TPCD, type(TPCD)))
                        if eyear_analysis_historical > 2005 and scn == 'rcp85': 
                            dict_threshold_rcp85['all_all_{:03}.{}'.format(regionID, region_name)] = threshold 

                        # record ----------
                        timeserieses.append(timeseries)  # aSrc (nyear)
                        if PDFs2NC_OVERALL: pdf2nc(PDFs, dict_info)

                        update_dict_stats(dict_stats, scn, soc, co2, region_name, 'all', 'all')
                        if smplType == 'Average':
                            allocate_item_2_dict_stats(dict_stats, 'PDFs', PDFs, dict_info)
                            allocate_item_2_dict_stats(dict_stats, 'spatial_average', timeseries, dict_info)
                        elif smplType == 'Extreme':
                            allocate_item_2_dict_stats(dict_stats, 'historical_percentile', hist_pctl_day, dict_info)

                        df_TPCD['{}.{}'.format(regionID,region_name)]['all_all'] = TPCD
                        add2df_stats(df_stats, 'all', 'all', '{:02}.{}'.format(regionID,region_name), TPCD, threshold, hist_min, hist_max, hist_pctl_day, timeseries)
                        mmbrChart[0, iregion] = TPCD

                        add2df_unprecedented_info_median(df_unprecedented_info_median, region_name, n_event, n_year, df_unprecedented_info_member)

                        if not region_name == 'GLB':
                            EnsmOUT[YY, XX] = TPCD
                            # update dict_unexceptional_global_mask
                            if unexceptional_mask is not None:
                                for iy, ix in zip(YY,XX):
                                    dict_unexceptional_global_mask['all_all'][:,iy,ix] = unexceptional_mask

                        if Figure1: updatefig1_addaxes_and_imshow(PDFs, fig1, ax01_pos, kdfSampling, dict_info)             # PDF timeseries on basemap
                        if Figure2: updatefig2_addaxes_and_imshow(PDFs, fig2, iregion, kdfSampling, dict_info, df_vrange)   # PDF timeseries in the list style
                        if Figure3: figure3('all', 'all', timeseries, TPCD, threshold, hist_min, hist_max, dict_info)       # plot of timeseries of aSrc (nYear)  ensemble
                        if Figure4: figure4('all', 'all', PDFs, TPCD, dict_info)                                            # PDF time series   ensemble
                        if Figure11: add_timeseries2axis('absolute', timeserieses, fig11, dictAxes11[regionID], dict_info)  # plot of timeseries of ensemble aSrcs on a basemap
                        if Figure12: update_indevidual_timeseries(12, timeserieses, dict_info)                              # plot of timeseries of ensemble aSrcs for each region
                        if Figure13: add_timeseries2axis('anomaly', timeserieses, fig12, dictAxes12[regionID], dict_info)   # plot of timeseries of anomaly in ensemble aSrcs on a basemap

                        del timeseries, gridarea, PDFs
                        if GC: gc.collect()
                        if MEMORY_USAGE_NOISY:  
                            print("\n\n====== in main\n{}{: >45}{}{: >10}{}".format('|','Variable Name','|','  Size','|'))
                            print(" -------------------------- ")
                            for k, v in locals().items():
                                if hasattr(v, 'size') and not k.startswith('_') and not isinstance(v,types.ModuleType):
                                    print("{}{: >45}{}{: >10}{}".format('|',k,'|',str(v.size),'|'))
                                elif hasattr(v, '__len__') and not k.startswith('_') and not isinstance(v,types.ModuleType):
                                    print("{}{: >45}{}{: >10}{}".format('|',k,'|',str(len(v)),'|'))
                        print('\ncheck >> memory usage: {}'.format(psutil.virtual_memory().percent))

                # ------------------------------------------------------------------------------------------------------
                print('\n'*3)
                EarliestOUT = mmbrOUT.min(axis=(0,1))
                LatestOUT = mmbrOUT.max(axis=(0,1))
                EnsmOUT = ma.masked_array(ma.masked_equal(EnsmOUT, 1e+20), mask=np.resize(landsea_or_greenland_mask, EnsmOUT.shape))  #(nY,nX)
                mmbrOUT = ma.masked_array(ma.masked_equal(mmbrOUT, 1e+20), mask=np.resize(landsea_or_greenland_mask, mmbrOUT.shape))  #(nGCM,nGHM,nY,nX)
                EarliestOUT = ma.masked_array(ma.masked_equal(EarliestOUT,1e+20), mask=np.resize(landsea_or_greenland_mask, EarliestOUT.shape))  #(nY,nX)
                LatestOUT =   ma.masked_array(ma.masked_equal(LatestOUT,  1e+20), mask=np.resize(landsea_or_greenland_mask, LatestOUT.shape))  #(nY,nX)
                StdDev = np.std(mmbrOUT, axis=(0,1))                                                                            #(nY,nX)

                # write out list in Excel format
                df_TPCD = df_TPCD.replace(1e+20, np.nan)
                df_stats = df_stats.replace(1e+20, np.nan)
                df_pop_ensmbl = get_df_pop(scn, dict_unexceptional_global_mask)  # TODO
                dict_pops[scn] = df_pop_ensmbl                                   # TODO

                if TPCD2EXCEL:
                    outExcel = 'TPCDs.{}.{:03}_{}.{}.{}.{}.{}_{}_{}.xlsx'.format(index,tChunk,threshold_type,drought_paras,mmeType,smplType,scn,soc,co2)
                    outpath = os.path.join(fig_directory, outExcel)
                    writer = pd.ExcelWriter(outpath)
                    # ---
                    df_TPCD.to_excel(writer,'TPCD table')
                    df_stats.to_excel(writer,'data')
                    df_pop_ensmbl.T.to_excel(writer,'totPop(ensmbl)')
                    df_unprecedented_info_member.to_excel(writer, 'unprecedented_info_member')
                    df_unprecedented_info_median.to_excel(writer, 'unprecedented_info_median')
                    df_unprecedented_flag_member.to_excel(writer, 'unprecedented_flag')
                    # ---
                    writer.save()
                    print('write out: {}'.format(outpath))

                if Figure2 and not READ_VMAX_OVERALL:
                    outExcel = 'vrange_of_PDFs.{}.{}.{}.{}.{}_{}_{}.xlsx'.format(index,drought_paras,mmeType,smplType,scn,soc,co2)
                    outpath = os.path.join(fig_directory, outExcel)
                    writer = pd.ExcelWriter(outpath)
                    df_vrange.to_excel(writer,'vrange')
                    writer.save()
                    print('write out: {}'.format(outpath))

                # write out results in binary
                if TPCD2BIN:
                    writeoutinbinary(EnsmOUT,  'ensemble{}.{}.{}.{:03}.{}.{}.{}'.format(mmeType,smplType,index,tChunk,drought_paras,soc,co2), fig_directory)
                    writeoutinbinary(mmbrOUT,      'member.{}.{}.{:03}.{}.{}.{}'.format(        smplType,index,tChunk,drought_paras,soc,co2), fig_directory)
                    writeoutinbinary(mmbrChart, 'mmbrChart.{}.{}.{:03}.{}.{}.{}'.format(        smplType,index,tChunk,drought_paras,soc,co2), fig_directory)

                # draw maps
                EnsmOUT = EnsmOUT[:293]
                EarliestOUT = EarliestOUT[:293]
                LatestOUT = LatestOUT[:293]
                mmbrOUT = mmbrOUT[:,:,:293,:]
                StdDev = StdDev[:293]
                mmbrChart = ma.masked_equal(mmbrChart, 1e+20)

                if Figure7: figure7(EnsmOUT,      'Ensmbl',   df_pop_ensmbl, dict_parameters)  # TPCD global map
                if Figure8: figure89(StdDev,      'STD',      None,          dict_parameters)  # TPCD global map (STD)
                if Figure9: figure89(EarliestOUT, 'Earliest', None,          dict_parameters)  # TPCD global map (Earliest&latest TPCD among members)
                if Figure9: figure89(LatestOUT,   'Latest',   None,          dict_parameters)  # TPCD global map (Earliest&latest TPCD among members)
                if Figure10: figure10(mmbrChart, dictRegion, RegionIDs, dict_parameters)

                if Figure1: close_figureobject(1, dict_parameters)  # PDF timeseries on basemap
                if Figure2: close_figureobject(2, dict_parameters)  # PDF timeseries in the list style

            if Figure11:  # plot of timeseries of ensemble aSrcs (Global&Regional) on a basemap
                if GLOBAL:                 custum_timeseries_axis(fig11, dictAxes11[0],        'absolute', 0)
                for regionID in RegionIDs: custum_timeseries_axis(fig11, dictAxes11[regionID], 'absolute', regionID)
                close_figureobject(11, dict_parameters)
            if Figure12:  # plot of timeseries of ensemble aSrcs for each region
                for regionID in RegionIDs:
                    region_name = dictRegion[regionID]
                    custum_indevidual_timeseries_axis(12, regionID, region_name)
                    close_figureobject12(12, dict_parameters, regionID, region_name)
            if Figure13:  # plot of timeseries of anomaly in ensemble aSrcs (Global&Regional) on a basemap
                if GLOBAL:                 custum_timeseries_axis(fig12, dictAxes12[0],      'anomaly', 0)
                for regionID in RegionIDs: custum_timeseries_axis(fig12,dictAxes12[regionID],'anomaly', regionID)
                close_figureobject(13, dict_parameters)

        if Figure14: figure13(dict_pops, dict_parameters)  # population plot

        if Figure15:  # colorfulplot_shiftingPDF
            for scn, regionID, ghm, gcm, sample_year in itertools.product(scns, RegionIDs, ghms, gcms, sample_years):
                #avoid_these = [0,1,2,11,21,18]  # TODO: just tentative!!
                avoid_these = []
                if regionID == 0: region_name = 'GLB'
                else: region_name = dictRegion[regionID]
                if regionID in avoid_these or scn == 'rcp26':
                    pass
                else:
                    dict_info = {'scn':scn, 'soc': soc, 'co2': co2,
                                 'regionID':regionID, 'region_name':region_name, 
                                 'ghm':ghm, 'gcm':gcm, 'sample_year':sample_year, 
                                 'dict_parameters':dict_parameters,
                                 'fig_directory_common': fig_directory_common}
                    draw_colorfulplot_and_shiftingPDF(dict_stats[scn][soc][co2][region_name][ghm][gcm], dict_info)

    print('\n\n=====================================\nThe process successfully finished!! (^.^)b')
    print('Multiprocessing: {}'.format(multiprocessing))
    if multiprocessing: print('(nMP: {})'.format(nMP))
    endTime  = datetime.datetime.now()
    diffTime = endTime - strTime
    print('took {} min in total.'.format(int(diffTime.seconds/60)))


if __name__=='__main__':
    main(*sys.argv)
