#!/usr/bin/env python
# To
# By Yusuke Satoh
# On
import os
import sys
import re
import time
import random
import datetime
import itertools
import pathlib
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
from sklearn.utils import resample
from netCDF4 import Dataset
from tqdm import tqdm
from utiltools import get_region_info

today = datetime.date.today().strftime('%Y%m%d')
plt.rcParams["font.size"] = 20


# --- output options  (Basically, for REUSE=False)
WRITE_NC = True
EXCEL_SUMARRY = True
EXCEL = False  # excel writer takes time...

# --- reuse option to generate figures
REUSE = False
#REUSE = True; reuse_date = 20210419

if REUSE:
    EXCEL = False  # excel writer takes time...
    EXCEL_SUMARRY = False
    WRITE_NC = False
    Figure1 = True
    Figure2 = False
    Figure3 = False
else: 
    Figure1 = False
    Figure2 = False
    Figure3 = False


# --- input from arguments

experiment = sys.argv[1]
sample_type = sys.argv[2]
region_type = sys.argv[3]

digit = int(sys.argv[4])
if digit == 3: K = 1000
elif digit == 4: K = 10000
elif digit == 5: K = 100000
else: raise ValueError

syear_org = int(sys.argv[5])

# --- key parameters of this analysis
tChunk = int(sys.argv[6])
#tChunk = 5
#tChunk = 10
#tChunk = 15
#tChunk = 20
#tChunk = 30

season = sys.argv[7]

case = f'case_{sys.argv[8]}'


# --- fixed parameters
if region_type == 'AR6_regions':
    dict_input_date = {
        ('basic', 'annual', 'Average'):     20210501,
        ('basic', 'DRY', 'Average'):        20210624,
        ('basic', 'WET', 'Average'):        20210630,
        #('basic', 'annual', 'Mean'):        20210707,  # AR6_regions  
        #('basic', 'DRY', 'Mean'):           20210707,  # AR6_regions 
        #('basic', 'WET', 'Mean'):           20210707,  # AR6_regions
        ('basic', 'annual', 'Mean'):        20211114,  # hydroregions  
        ('basic', 'DRY', 'Mean'):           20211114,  # hydroregions 
        ('basic', 'WET', 'Mean'):           20211114,  # hydroregions
        ('co2exp', 'annual', 'Mean'):       20210707,
        ('co2exp', 'DRY', 'Mean'):          20210707,
        ('co2exp', 'WET', 'Mean'):          20210707,
        ('rcp26soc', 'annual', 'Mean'):     20210707,
        ('rcp26soc', 'DRY', 'Mean'):        20210707,
        ('rcp26soc', 'WET', 'Mean'):        20210707,
        ('picontrol', 'annual', 'Average'): 20210505,
        ('picontrol', 'DRY', 'Average'): '',
        ('picontrol', 'WET', 'Average'): '',
        #('basic', 'annual', 'Mean'):        20210722,  # with DryGridMask
        #('basic', 'DRY', 'Mean'):           20210722,  # with DryGridMask  
        #('basic', 'WET', 'Mean'):           20210722,  # with DryGridMask
        #('co2exp', 'DRY', 'Mean'):          20210723,  # with DryGridMask
        #('rcp26soc', 'DRY', 'Mean'):        20210723,  # with DryGridMask
        }
elif region_type == 'HydroBASINS_lev1':
    dict_input_date = {  # date of TFE analysis
        ('basic', 'annual', 'Mean'): 20220102,
        ('basic', 'DRY', 'Mean'):    20220102,
        ('basic', 'WET', 'Mean'):    20220102,
        ('co2exp', 'DRY', 'Mean'):   20220102,
        ('rcp26soc', 'DRY', 'Mean'): 20220102,
        }
elif region_type == 'hydroregions':
    pass  # TODO
else:
    raise ValueError(f'check {region_type}')


threshold = 'max'

block_size = 5  # default
#block_size = 1

trend_type = 'quadratic'  # this is the one!!
#trend_type = 'linear'
#trend_type =  'moving'  # moving average

SHUFFLE = 'None'  # 1st option
#SHUFFLE = 'type1'  # 2nd option
#SHUFFLE = 'type2'

rcps = [sys.argv[9]]
socs = [sys.argv[10]]
co2s = [sys.argv[11]]
gcms = ['hadgem2-es', 'ipsl-cm5a-lr', 'gfdl-esm2m', 'miroc5']  # full members
ghms = ['cwatm', 'h08', 'lpjml', 'matsiro', 'watergap2']  # full members
if experiment == 'co2exp':
    ghms = ['lpjml', 'matsiro']
elif experiment == 'rcp26soc':
    ghms = ['cwatm', 'h08', 'lpjml', 'matsiro']
else:
    pass
"""
# --- basic
rcps = ['rcp85', 'rcp26']
socs = ['2005soc']
co2s = ['co2']
gcms = ['hadgem2-es', 'ipsl-cm5a-lr', 'gfdl-esm2m', 'miroc5']  # full members
ghms = ['cwatm', 'h08', 'lpjml', 'matsiro', 'watergap2']  # full members
# --- for experiment
if experiment == 'picontrol': 
    rcps = ['picontrol', 'rcp85', 'rcp26']
elif experiment == 'co2exp':
    ghms = ['lpjml', 'matsiro']
    co2s = ['co2', '2005co2']
elif experiment == 'rcp26soc': 
    ghms = ['cwatm', 'h08', 'lpjml', 'matsiro']
    rcps = ['rcp26']
    socs = ['2005soc', 'rcp26soc']
"""

tfe_setting = f'{experiment}.abs.05.{sample_type[:3]}.1865-2005.{threshold}.05.historical_histsoc_co2_1861-2005.1865-2099.plt1865-'
drought_setting = f'Q80win15Len30tau4{season}'

print('\n'+'==='*3)
for setting in ['experiment', 'threshold', 'K', 'trend_type', 'tChunk', 'case']:
    print(f'{setting:>15}: {eval(setting)}')
print('==='*3)
for setting in ['tfe_setting', 'drought_setting', 'rcps', 'socs', 'co2s', 'ghms', 'gcms']:
    print(f'{setting:>15}: {eval(setting)}')
print('==='*3+'\n')


if region_type == 'AR6_regions':
    dict_regions = {
         1: 'NWN',   2: 'NEN',  3: 'WNA',   4: 'CNA',   5: 'ENA',
         6: 'NCA',   7: 'SCA',  8: 'CAR',   9: 'NWS',  10: 'NSA',
        11: 'NES',  12: 'SAM', 13: 'SWS',  14: 'SES',  15: 'SSA',
        16: 'NEU',  17: 'WCE', 18: 'EEU',  19: 'MED',  20: 'SAH',
        21: 'WAF',  22: 'CAF', 23: 'NEAF', 24: 'SEAF', 25: 'WSAF',
        26: 'ESAF', 27: 'MDG', 28: 'RAR',  29: 'WSB',  30: 'ESB',
        31: 'RFE',  32: 'WCA', 33: 'ECA',  34: 'TIB',  35: 'EAS',
        36: 'ARP',  37: 'SAS', 38: 'SEA',  39: 'NAU',  40: 'CAU',
        41: 'EAU',  42: 'SAU', 43: 'NZ',
        #19: 'MED'  # for test
        }
elif region_type == 'HydroBASINS_lev1':
        
    region_map, dict_regions = get_region_info(region_type, TEST=False)
    del region_map
    """
    region_excel_path = '/data/rg001/sgec0017/data/mapmask/HydroSHEDS/HydroBASINS/withLakes/Global_merge/hybas_lake____lev02_v1c_merge.xlsx'
    df_region = pd.read_excel(region_excel_path, header=0, usecols=[0,1])
    df_region = df_region.drop_duplicates()
    dict_regions = {}
    n_regions = df_region.shape[0]
    for iregion in range(n_regions):
        region_ID, long_name = df_region.iloc[iregion]
        if region_ID != 0:  # skip Grennland
            dict_regions[region_ID] = long_name
    del dict_regions[57]  # skip Australia_4 (because therea are no grid at 0.5deg...)
    """
elif region_type == 'hydroregions':
    dict_regions = {
        11: 'BOR_Nam',  12: 'NML_Nam',  13: 'NDR_Nam',  14: 'NST_Nam',  25: 'EQT_Sam',  
        26: 'SST_Sam',  27: 'SDR_Sam',  28: 'SML_Sam',  31: 'BOR_Eur',  32: 'NML_Eur',  
        43: 'NDR_Afr',  44: 'NST_Afr',  45: 'EQT_Afr',  46: 'SST_Afr',  47: 'SDR_Afr',  
        48: 'SML_Afr',  511: 'BOR_Asi(WSb)',  512: 'BOR_Asi(ESb)',  52: 'NML_Asi',  531: 'NDR_Asi(MdE)',  
        532: 'NDR_Asi(CAs)',  54: 'NST_Asi',  55: 'EQT_Asi',  66: 'SST_Aus',  67: 'SDR_Aus',  
        68: 'SML_Aus',  
        }
else:
    raise ValueError('something is wrong with region_type...')


main_directory = '/data/rg001/sgec0017/data/figure_box/drought_tpcd_isimip2b'

excel_directory = os.path.join(
    main_directory, 'TPCD.estimator.scns.v2_for_1stRevision',
    region_type, 'all_all', drought_setting,
    tfe_setting)

syear_org, eyear_org = syear_org, 2099
years_org = np.arange(syear_org, eyear_org+1)
if not trend_type == 'moving':
    moving_sample_size = 0
    syear, eyear = syear_org, eyear_org
elif trend_type == 'moving':
    moving_sample_size = 21
    weight = np.ones(moving_sample_size) / moving_sample_size
    window_size_half = int((moving_sample_size-1)/2)
    syear, eyear = syear_org+window_size_half, syear_org-window_size_half
years = np.arange(syear, eyear + 1)

figure_directory = os.path.join(
    main_directory, 'bootstrap', region_type, 
    '{}.{}.{}{}{}.{}.block{}.{}'.format(experiment, sample_type, 
                                        trend_type, str(syear_org)[2:], str(eyear_org)[2:],
                                        K, block_size, season), 
    'tChunk_{:02}'.format(tChunk), 
    today)
    #case, today)
if not os.path.isdir(figure_directory): os.makedirs(figure_directory)
pathlib.Path(os.path.join(figure_directory, f'originalfile__{tfe_setting}')).touch()

if REUSE:
    reuse_data_directory = os.path.join(
        main_directory, 'bootstrap', region_type,
        '{}.{}.samplesize_{}.{}-{}.blocksize_{}.{}'.format(sample_type, trend_type, K, 
                                                           syear_org, eyear_org, 
                                                           block_size, season), 
        'tChunk_{:02}'.format(tChunk), 
        str(reuse_date))
        #case, str(reuse_date))
    #figure_directory = reuse_data_directory  

pT = re.compile('T+')
pF = re.compile('F+')
missing_value = 9999
edges = np.arange(syear, missing_value+3) - 0.5  # edges for histograms


def find_tfe(aSrc, threshold, _years):
    
    syear_tfe_analysis = 2005
    aSrc = aSrc[_years.tolist().index(syear_tfe_analysis):]   # use data after syear_tfe_analysis

    # --- detect unprecedented (find worsen cases)
    #   True:  unprecedented!! = abnormal = more than the hiscorical max
    #   False: normal, within historical range.
    unprecedented_mask = np.ma.make_mask(aSrc > threshold)
    unprecedented = ''.join(['F' if i == b'F' else 'T' for i in unprecedented_mask.astype('S1')])  # 'FFFFFFTTTTTTFFTTTFFFFFFFF...'   'T' means 'unprecedented'
    if len(unprecedented) == 1:
        unprecedented = 'F' * aSrc.shape[0]

    # --- Caution!! At this point, the first item in "unprecedented" is NOT always F in case you do not select the entire historical preiod as the reference period.
    if unprecedented[0] == 'F':
        initial_flug = 'F'
    elif unprecedented[0] == 'T':
        initial_flug = 'T'
    if unprecedented[-1] == 'T':
        might_be_permanent = True
    else:
        might_be_permanent = False

    # --- split the seaquence
    Fs = pF.findall(unprecedented)  # ['FFFFFF', 'FF', 'FFFFFFFF', ...]
    Ts = pT.findall(unprecedented)  # ['TTTTTT', 'TTT, ... ]

    if tChunk == 99:
        if might_be_permanent:  # the last step (2099) is unprecedented
            if len(Ts[-1]) >= 15:  # in use
                Ts = ['F' * len(ts)  for ts in Ts[:-1]] + [Ts[-1]]  # ['FFFFFF', 'FFF', .... 'FFF', 'TTTTTTTTTTTTTTT']
            else:  # reject
                Ts = ['F' * len(ts)  for ts in Ts]  # ['FFFFFF', 'FFF', ...]
        else:  # might_be_permenant = False. So, reject all Ts.
            Ts = ['F' * len(ts)  for ts in Ts]  # ['FFFFFF', 'FFF', ...]
    else:  # remove shorter unprecedented period than tChunk
        Ts = ['F' * len(ts) if len(ts) < tChunk else ts for ts in Ts]  # ['TTTTTT', 'FFF', ...]
        #Ts = ['F' * len(ts) if len(ts) <= tChunk else ts for ts in Ts]  # ['TTTTTT', 'FFF', ...]

    # --- merge
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

    # --- find tfe
    # --- get all non-future years False.
    #_unprecedented_mask = [False] * (years.index(years_analysis_historical[-1]) + 1) + unprecedented_mask[years.index(years_analysis_future[0]):]
    _unprecedented_mask = unprecedented_mask
    if True in _unprecedented_mask:
        # --- find an index of the first True in the future period
        tfe = int(range(syear_tfe_analysis, 2099+1)[_unprecedented_mask.index(True)])
    else:
        #tfe = 1e+20
        tfe = missing_value

    # --- make a mask to calculate the potential population affected under the unprecedented drought condition.
    #     the mask maskes out NON unprecedented grids. -->  If NOT recordbreaking, then this mask is True.
    #     shorter recordbreaking condition than tChunk are removed here.
    unexceptional_mask = np.full(aSrc.shape, False)
    unexceptional_mask[np.where(np.array(unprecedented_mask) == False)] = True  # [True ... False False False ...]

    _Ts = [len(ts) for ts in Ts if 'T' in ts]
    #n_recordbreaking_event = len(_Ts)
    n_recordbreaking_year = sum(_Ts)

    return tfe, unexceptional_mask, n_recordbreaking_year  # for each member


def random_sample(src):
    random.seed(int(time.time()))
    #return random.sample(list(src), src.shape[0])
    src_input = src
    src_output = random.sample(list(src), src.shape[0])
    print(f'random_sample: {src_input} ---> {src_output}')
    return src_output


def write_nc(src, topic, rcp, soc, co2, region):

    out_file = f'{topic}.{region}.nc4'
    outDir = os.path.join(figure_directory, f'{rcp}_{soc}_{co2}')
    if not os.path.isdir(outDir): os.makedirs(outDir)
    outPath = os.path.join(outDir, out_file)

    rootgrp = Dataset(outPath, 'w', format='NETCDF4')
    import time
    rootgrp.title = f'bootstrap resampling output: {topic} for {region}'
    rootgrp.description = 'input'.format()
    rootgrp.history = 'Created ' + time.ctime(time.time())
    rootgrp.source = f'ISI-MIP2b TFE analysis: {drought_setting}/{tfe_setting}'
    rootgrp.institution = 'NIES'
    rootgrp.contact = 'satoh.yusuke@nies.go.jp'

    if topic == 'resampled_TS':
        print(f'--- {topic} ---')
        print(f'src.shape: {src.shape}')
        print(f'ensemble_member: {len(ghms)*len(gcms)}')
        print(f'resampled_size: {K+1}')
        print(f'year: {years.shape[0]}')
        # --- create dimension
        ensemble_member = rootgrp.createDimension('ensemble_member', size=(len(ghms)*len(gcms)))
        resampled_size = rootgrp.createDimension('resampled_size', size=(K+1))
        year = rootgrp.createDimension('year', size=years.shape[0])
        # --- create variable
        _ensemble_members = rootgrp.createVariable('ensemble_member', 'int', ('ensemble_member'))
        _resampled_sizes = rootgrp.createVariable('resampled_size', 'int', ('resampled_size'))
        _years = rootgrp.createVariable('year', 'int', ('year'))
        srcs = rootgrp.createVariable(topic, 'f4', ('ensemble_member', 'resampled_size', 'year'),
                                      zlib=True, complevel=5, fill_value=np.float32(1e+20),
                                      chunksizes=(len(ghms)*len(gcms), K+1, years.shape[0]))
        srcs.long_name = 'bootstrap-resampled time series of regional average FDD ({}-{})'.format(syear, eyear) # longName
        srcs.standard_name = 'resampled time series' # standard_name
        srcs.missing_value = np.float32(1e+20)
        srcs.unit = 'day/year'  # unitName
        srcs.members = f'ghms: {ghms}, gcms: {gcms}'
        # --- allocate data
        _years[:] = years
        _resampled_sizes[:] = range(K+1)
        _ensemble_members[:] = range(len(ghms)*len(gcms))
        srcs[:] = src

    elif topic == 'unrecordbreaking_mask':
        years_future = range(2005, 2099+1)
        print(f'--- {topic} ---')
        print(f'src.shape: {src.shape}')
        print(f'ensemble_member: {len(ghms)*len(gcms)}')
        print(f'resampled_size: {K+1}')
        print(f'year: {len(years_future)}')
        # --- create dimension
        ensemble_member = rootgrp.createDimension('ensemble_member', size=(len(ghms)*len(gcms)))
        resampled_size = rootgrp.createDimension('resampled_size', size=(K+1))
        year = rootgrp.createDimension('year', size=len(years_future))
        # --- create variable
        _ensemble_members = rootgrp.createVariable('ensemble_member', 'int', ('ensemble_member'))
        _resampled_sizes = rootgrp.createVariable('resampled_size', 'int', ('resampled_size'))
        _years = rootgrp.createVariable('year', 'int', ('year'))
        srcs = rootgrp.createVariable(topic, 'f4', ('ensemble_member', 'resampled_size', 'year'),
                                      zlib=True, complevel=5, fill_value=np.float32(1e+20),
                                      chunksizes=(len(ghms)*len(gcms), K+1, len(years_future)))
        srcs.long_name = f'time series of unrecordbreaking mask (2005-{eyear})' # longName
        srcs.standard_name = 'unrecordbreaking mask' # standard_name
        srcs.missing_value = np.float32(1e+20)
        srcs.members = f'ghms: {ghms}, gcms: {gcms}'
        srcs.memo = 'Mask Treu: unrecord breaking, False: recordbreaking'
        # --- allocate data
        _years[:] = years_future
        _resampled_sizes[:] = range(K+1)
        _ensemble_members[:] = range(len(ghms)*len(gcms))
        srcs[:] = src

    elif topic == 'resampled_tfe' or topic == 'n_recordbreaking_years':
        print(f'--- {topic} ---')
        print(f'src.shape: {src.shape}')
        print(f'ensemble_member: {len(ghms)*len(gcms)}')
        print(f'resampled_size: {K+1}')
        # --- create dimension
        resampled_size = rootgrp.createDimension('resampled_size', size=(K+1))
        ensemble_member = rootgrp.createDimension('ensemble_member', size=(len(ghms)*len(gcms)))
        # --- create variable
        _resampled_sizes = rootgrp.createVariable('resampled_size', 'int', ('resampled_size'))
        _ensemble_members = rootgrp.createVariable('ensemble_member', 'int', ('ensemble_member'))
        srcs = rootgrp.createVariable(topic, 'f4', ('ensemble_member', 'resampled_size'),
                                      zlib=True, complevel=5, fill_value=np.float32(1e+20),
                                      chunksizes=(len(ghms)*len(gcms), K+1))
        if topic == 'resampled_tfe':
            srcs.long_name = 'tfe from bootstrap-resampled time series' # longName
            srcs.standard_name = 'resampled tfe' # standard_name
            srcs.missing_value = np.float32(1e+20)
            srcs.unit = 'year'  # unitName
            srcs.members = f'ghms: {ghms}, gcms: {gcms}'
        elif topic == 'n_recordbreaking_years':
            srcs.long_name = f'the number of recordbreaking years during {syear}-{eyear}' # longName
            srcs.standard_name = 'N of recordbreaking years' # standard_name
            srcs.missing_value = np.float32(1e+20)
            srcs.unit = 'years'  # unitName
            srcs.members = f'ghms: {ghms}, gcms: {gcms}'
        # --- allocate data
        _resampled_sizes[:] = range(K+1)
        _ensemble_members[:] = range(len(ghms)*len(gcms))
        srcs[:] = src

    # --- close
    rootgrp.close()
    print(f'\nsave: {outPath}')


# --------------------------------------------------------------------------------------------------------------------------
def main():
    print(f'start this process at {datetime.datetime.now()}')
    start_time_overall = time.time()

    # --- set seed for the random sampling to reproduce it repeatedly.
    SEED = K * int(case[-1])
    np.random.seed(SEED)
    if SHUFFLE == 'type1' or SHUFFLE == 'type2': 
        random.seed(SEED)

    # regions
    regions = [f'{region_id:02}.{region_name}' for region_id, region_name in dict_regions.items()]
    print(f'regions: {regions}')

    for rcp, soc, co2 in itertools.product(rcps, socs, co2s):

        # --- load input 
        excel_path = os.path.join(excel_directory, str(dict_input_date[(experiment, season, sample_type)]), 
                                  f'TPCDs.nDayTot.005_{threshold}.{drought_setting}.all.{sample_type}.{rcp}_{soc}_{co2}.xlsx')
        df_input = pd.read_excel(excel_path, sheet_name='data', index_col=0)
        print(f'df_input:\n{df_input}')
        indices = list(df_input.index)
        columns = list(df_input.columns)

        # --- generate an excel file for summary
        summary_excel_file = f'bootstrap_summary_{rcp}.{soc}.{co2}.xlsx'
        summary_excel_path = os.path.join(figure_directory, summary_excel_file)
        writer_summary = pd.ExcelWriter(summary_excel_path)
        df_summary = pd.DataFrame(index=regions, columns=['median', 'Q10', 'Q25', 'Q75', 'Q90'])

        for region in regions:
            print('\n>>> process for {}'.format(region))
            start_time_region = time.time()

            if Figure1:
                fig1 = plt.figure(num=1, figsize=(45, 16))
                fig1.suptitle(f'{region}  (K={K})')
                plt.subplots_adjust(left=0.03, bottom=0.05, right=0.99, top=0.95, wspace=0.075)  #, hspace=None)
                gs1 = gridspec.GridSpec(3, 4)
                ax1_1 = plt.subplot(gs1[:2, 0]); ax1_1.set_title('(1) original time series')
                ax1_2 = plt.subplot(gs1[0, 1]);  ax1_2.set_title('(2) trend')
                ax1_3 = plt.subplot(gs1[1, 1]);  ax1_3.set_title('(3) without trend')
                ax1_4 = plt.subplot(gs1[0, 2]);  ax1_4.set_title('(4) trend (cwatm-hadgem)')
                ax1_5 = plt.subplot(gs1[1, 2]);  ax1_5.set_title('(5) without trend (cwatm-hadgem)')
                ax1_6 = plt.subplot(gs1[0, 3]);  ax1_6.set_title('(6) resampled (cwatm-hadgem)')
                ax1_7 = plt.subplot(gs1[1, 3]);  ax1_7.set_title('(7) resampled with trend (cwatm-hadgem)')
                ax1_8 = plt.subplot(gs1[2, :2]); ax1_8.set_title('(8) frequency of tfe')
                ax1_9 = plt.subplot(gs1[2, 2:]); ax1_9.set_title('(9) frequency of tfe')
            if Figure2:
                fig2 = plt.figure(num=2, figsize=(45, 30))
                fig2.suptitle('PDF of TFE')
                plt.subplots_adjust(left=0.05, bottom=0.05, right=0.99, top=0.95, wspace=0.1, hspace=0.2)
                gs2 = gridspec.GridSpec(5, 4)
            if Figure3:
                fig3 = plt.figure(num=3, figsize=(45, 30))
                fig3.suptitle('CDF of TFE')
                plt.subplots_adjust(left=0.05, bottom=0.05, right=0.99, top=0.95, wspace=0.1, hspace=0.2)
                gs3 = gridspec.GridSpec(5, 4)

            if EXCEL:  # optional, but not recommended because this option gets the process very slow...
                excel_file = f'bootstrap_{rcp}.{soc}.{co2}_{region}.xlsx'
                output_excel_path = os.path.join(figure_directory, excel_file)
                writer = pd.ExcelWriter(output_excel_path)

            tfes = []  # tfes from full member bootstrap results
            TSs = []   # resampled time series
            if WRITE_NC:
                URB_masks = []  # Unrecordbreaking mask
                N_RBYs = []     # N of recordbreaking years

            for (ighm, ghm), (igcm, gcm) in itertools.product(enumerate(ghms), enumerate(gcms)):
                print(f'\n{region:<6} <<< {4*ighm+igcm+1:02} {ghm} {gcm}')
                label = f'{ghm}_{gcm}'
                index_name = f'{ghm}_{gcm}_{region}'
                idx = indices.index(index_name)
                _tfes = []  # tfes of an ensemble member
                _TSs = []

                if WRITE_NC:
                    _URB_masks = []
                    _N_RBYs = []
                if EXCEL:
                    df_out  = pd.DataFrame(index=list(range(K+1))+['median'], columns=['TFE']+list(years_org))
                    df_out2 = pd.DataFrame(index=list(range(K+1)),            columns=['n_unprecetended_year']+list(years_org))

                # --- original
                src_org = df_input.iloc[idx, columns.index(syear_org):]
                src_org = src_org.replace('--', 0)  # this is for LPJmL in WET under rcp85
                src_org = src_org.astype('float').values  # (nyear_org_full)
                if Figure1:
                    plt.figure(num=1)
                    ax1_1.plot(years_org, src_org, label=label)
                if REUSE:  # reuse tfe info from a netCDF file
                    reuse_ts_nc_path  = os.path.join(reuse_data_directory, f'resampled_TS.{region}.nc4')
                    reuse_tfe_nc_path = os.path.join(reuse_data_directory, f'resampled_tfe.{region}.nc4')
                    tfe               = Dataset(reuse_tfe_nc_path)['resampled_tfe'][len(gcms)*ighm+igcm,0]
                else:  # normal
                    if not trend_type == 'moving':
                        tfe, unrecordbreaking_mask, n_recordbreaking_year = find_tfe(src_org, df_input.loc[index_name, 'threshold'], years_org)
                    else:  # moving
                        tfe, unrecordbreaking_mask, n_recordbreaking_year = find_tfe(src_org[window_size_half:-window_size_half], 
                                                                                     df_input.loc[index_name, 'threshold'], years)
                    tfe_org = tfe
                    if EXCEL:
                        df_out.iloc[0, 0] = tfe
                        df_out.iloc[0, 1:] = src_org
                        df_out2.iloc[0, 0] = n_recordbreaking_year
                        if not trend_type == 'moving':
                            df_out2.iloc[0, 1:] = unrecordbreaking_mask.astype('int')  # 1 if not recordbreaking, 0 is recordbreaking
                        else:  # moving
                            df_out2.iloc[0, 1+window_size_half:] = unrecordbreaking_mask.astype('int')  # 1 if not recordbreaking, 0 is recordbreaking
                _tfes.append(tfe)
                if not trend_type == 'moving':
                    _TSs.append(src_org)
                else: 
                    _TSs.append(src_org[window_size_half:-window_size_half])

                if WRITE_NC:
                    _URB_masks.append(unrecordbreaking_mask)
                    _N_RBYs.append(n_recordbreaking_year)

                # trend and variability
                if trend_type == 'linear':  # --- liner and diff
                    cf1, cf0 = np.polyfit(years, src_org, deg=1)
                    src_trend = cf1 * years + cf0
                    src_diff = src_org - src_trend
                elif trend_type == 'quadratic':
                    cf2, cf1, cf0 = np.polyfit(years, src_org, deg=2)
                    src_trend = cf2 * years**2 + cf1 * years + cf0     # (nyear_org_full) or (nyear_org_full+window)
                    src_diff = src_org - src_trend                     # (nyear_org_full) or (nyear_org_full+window)
                elif trend_type == 'moving':
                    src_trend = np.convolve(src_org, weight, mode='valid')
                    src_diff = src_org[window_size_half:-window_size_half] - src_trend

                if Figure1:
                    ax1_2.plot(years, src_trend)
                    y2_min, y2_max1_ = ax1_2.get_ylim()
                    ax1_3.plot(years, src_diff)
                    y3_min, y3_max1_ = ax1_3.get_ylim()
                    if ighm == 0 and igcm == 0:
                        _src_trend = np.copy(src_trend)
                        _src_diff  = np.copy(src_diff)

                # --- Blockwise Bootstrap --------------------------------------------------------------------------------------------------------
                start_time = time.time()
                if REUSE:
                    _TSs = _TSs + Dataset(reuse_ts_nc_path)['resampled_TS'][len(gcms)*ighm+igcm, 1:].tolist()
                    _TSs_without_trend = _TSs - src_trend
                    _tfes = _tfes + Dataset(reuse_tfe_nc_path)['resampled_tfe'][len(gcms)*ighm+igcm, 1:].tolist()
                else:
                    #for counter in tqdm(range(K)):  # ramdom sampled process...
                    for counter in range(K):  # ramdom sampled process...
                        # resampling...
                        if not block_size == 1:
                            src_diff = src_diff.reshape(-1, block_size)
                        src_resampled = resample(src_diff)  # randam sampling
                        # shuffling...
                        if SHUFFLE == 'type1':  # suffle in each block
                            src_resampled = np.array([random.sample(list(block), block.shape[0]) for block in src_resampled])
                        elif SHUFFLE == 'type2':  # shuffle using different seed
                            src_resampled = np.array([random_sample(block) for block in src_resampled])
                        if not block_size == 1:
                            src_resampled = src_resampled.reshape(-1)
                        # resampled with trend
                        src_resampled_with_trend = src_trend + src_resampled
                        tfe, unrecordbreaking_mask, n_recordbreaking_year = find_tfe(src_resampled_with_trend, df_input.loc[index_name, 'threshold'], years)
                        _tfes.append(tfe)
                        _TSs.append(src_resampled_with_trend)
                        _URB_masks.append(unrecordbreaking_mask)
                        _N_RBYs.append(n_recordbreaking_year)
                        if EXCEL:
                            df_out.iloc[counter+1, 0] = tfe
                            df_out2.iloc[counter, 0] = n_recordbreaking_year
                            if not trend_type == 'moving':
                                df_out.iloc[counter+1, 1:] = src_resampled_with_trend
                                df_out2.iloc[counter, 1:] = unrecordbreaking_mask.astype('int')  # 1 if not recordbreaking, 0 is recordbreaking
                            else:
                                df_out.iloc[counter+1, 1+window_size_half:] = src_resampled_with_trend
                                df_out2.iloc[counter, 1+window_size_half:] = unrecordbreaking_mask.astype('int')  # 1 if not recordbreaking, 0 is recordbreaking
                # ---------------------------------------------------------------------------------------------------------------------------------

                if EXCEL:
                    elapsed_time = time.time() - start_time
                    print(f'Bootstrap took {elapsed_time:.3f} sec')
                    tfe_bootstrap_median = int(df_out.iloc[:-1, 0].median())
                    df_out.loc['median', 'TFE'] = tfe_bootstrap_median
                    print('bootstrap median tfe: {tfe_bootstrap_median}')

                if Figure2:
                    plt.figure(num=2)  # member-specific PDF
                    ax = plt.subplot(gs2[ighm, igcm])
                    ax.set_title(label)
                    N, bins, patches = ax.hist(_tfes, bins=edges, density=True)
                    if REUSE:
                        ax.axvline(x=Dataset(reuse_tfe_nc_path)['resampled_tfe'][len(gcms)*ighm+igcm,0], color='#ff0000', linewidth=2)
                    else:
                        ax.axvline(x=tfe_org, color='#ff0000', linewidth=2)
                    ax.axvline(x=np.median(_tfes), color='#000000', linewidth=2)
                    ax.set_xlim([syear-0.5, eyear+0.5])
                    ax.set_ylim(ymin=0, ymax=max(N[:100].max(), 0.00001))
                if Figure3:
                    plt.figure(num=3)  # member-specific CDF
                    ax = plt.subplot(gs3[ighm, igcm])
                    ax.set_title(label)
                    ax.hist(_tfes, bins=edges, density=True, cumulative=True)
                    if REUSE:
                        ax.axvline(x=Dataset(reuse_tfe_nc_path)['resampled_tfe'][len(gcms)*ighm+igcm,0], color='#ff0000', linewidth=2)
                    else:
                        ax.axvline(x=df_out.iloc[0,0], color='#ff0000', linewidth=2)
                    ax.axvline(x=np.median(_tfes), color='#000000', linewidth=2)
                    ax.set_xlim([syear-0.5, eyear+0.5])
                    ax.set_ylim(ymin=0, ymax=1.0)

                if EXCEL:
                    df_out.to_excel(writer, sheet_name=label)

                tfes.append(_tfes)
                TSs.append(_TSs)
                if WRITE_NC:
                    URB_masks.append(_URB_masks)
                    N_RBYs.append(_N_RBYs)


            # --- after the ensemble member loop
            if WRITE_NC:
                write_nc(np.array(tfes),      'resampled_tfe',          rcp, soc, co2, region)  # (n_member, K+1)
                write_nc(np.array(TSs),       'resampled_TS',           rcp, soc, co2, region)  # (n_member, K+1, years)
                write_nc(np.array(URB_masks), 'unrecordbreaking_mask',  rcp, soc, co2, region)  # (n_member, K+1, years)
                write_nc(np.array(N_RBYs),    'n_recordbreaking_years', rcp, soc, co2, region)  # (n_member, K+1)

            tfes = sum(tfes, [])  # flatten
            full_median_tfe = np.median(tfes)


            # --- stats of the region (based on full member bootstrap samples)
            if EXCEL_SUMARRY:
                df_summary.loc[region, 'median'] = np.median(tfes)
                df_summary.loc[region, 'Q10'] = np.percentile(tfes, 10)
                df_summary.loc[region, 'Q25'] = np.percentile(tfes, 25)
                df_summary.loc[region, 'Q75'] = np.percentile(tfes, 75)
                df_summary.loc[region, 'Q90'] = np.percentile(tfes, 90)

            if Figure1:
                plt.figure(num=1)
                # for the sample of resampling (ighm=0 & igcm=0)
                src_trend = _src_trend
                src_diff = _src_diff
                ax1_4.plot(years, src_trend)
                ax1_4.set_ylim([y2_min, y2_max1_])
                ax1_5.plot(years, src_diff)
                ax1_5.set_ylim([y3_min, y3_max1_])
                for counter in range(K):  # Blockwise Bootstrap
                    if REUSE:
                        src_resampled_with_trend = Dataset(reuse_ts_nc_path)['resampled_TS'][0, counter+1]
                        src_resampled = src_resampled_with_trend - src_trend
                        tfe = Dataset(reuse_tfe_nc_path)['resampled_tfe'][0, counter+1]
                    else:
                        src_resampled_with_trend = TSs[0,counter+1]
                        src_resampled = src_resampled_with_trend - src_trend
                    if counter < 10:
                        ax1_6.plot(years, src_resampled)
                        ax1_7.plot(years, src_resampled_with_trend)
                ax1_1.set_xlim([syear_org-0.5,2130])
                ax1_1.legend()
                for ax in [ax1_2, ax1_3, ax1_4, ax1_5, ax1_6, ax1_7]:
                    ax.set_xlim([syear-0.5, eyear+0.5])
                    
                # member median TFE
                member_tfes = []
                for ghm, gcm in itertools.product(ghms, gcms):
                    index_name = f'{ghm}_{gcm}_{region}'
                    tfe = df_input.loc[index_name, 'TPCD']
                    ax1_8.axvline(x=tfe, color='#606060', linewidth=0.5)  # right gray
                    member_tfes.append(tfe)
                member_median_tfe = np.median(member_tfes)

                # --- histgram
                start_time = time.time()
                ax1_8.hist(tfes, bins=edges, cumulative=False, density=True)
                ax1_8.axvline(x=member_median_tfe, color='#000000', linewidth=1.)
                ax1_8.axvline(x=full_median_tfe, color='#ff0000', linewidth=1.)
                ax1_8.set_xlim([syear-0.5, eyear+0.5])
                elapsed_time = time.time() - start_time
                print(f'PDF took {elapsed_time:.3f} sec')

                # --- histgram (cumulative)
                start_time = time.time()
                ax1_9.hist(tfes, bins=edges, cumulative=True, density=True)
                ax1_9.axvline(x=member_median_tfe, color='#000000', linewidth=1.)
                ax1_8.axvline(x=full_median_tfe, color='#ff0000', linewidth=1.)
                ax1_9.set_xlim([syear-0.5, eyear+0.5])
                ax1_9.set_ylim(ymax=1.0)
                elapsed_time = time.time() - start_time
                print(f'CDF took {elapsed_time:.3f} sec'.format())

                # save Fig1
                print('\nsaving figure1...')
                start_time_savefig = time.time()
                plt.figure(num=1)
                figure_name = f'overview.samplesize{K}.shuffle{SHUFFLE}.blocksize{block_size}.tChunk{tChunk:02}.{rcp}.{soc}.{co2}.{region}.{case}.png'
                figure_path = os.path.join(figure_directory, figure_name)
                plt.savefig(figure_path)
                plt.close(1)
                print(f'savefig: {figure_path}')
                elapsed_time = time.time() - start_time_savefig
                print(f'took {elapsed_time:.3f} sec')

            if Figure2:  # save Fig2
                print('saving figure2...')
                start_time_savefig = time.time()
                plt.figure(num=2)
                figure_name = 'pdf.samplesize{K}.shuffle{SHUFFLE}.blocksize{block_size}.tChunk{tChunk:02}.{rcp}.{soc}.{co2}.{region}.{case}.png'
                figure_path = os.path.join(figure_directory, figure_name)
                plt.savefig(figure_path)
                plt.close(2)
                print(f'savefig: {figure_path}')
                elapsed_time = time.time() - start_time_savefig
                print(f'took {elapsed_time:.3f} sec')

            if Figure3:  # save Fig3
                print('saving figure3...')
                start_time_savefig = time.time()
                plt.figure(num=3)
                figure_name = f'cdf.samplesize{K}.shuffle{SHUFFLE}.blocksize{block_size}.tChunk{tChunk:02}.{rcp}.{soc}.{co2}.{region}.{case}.png'
                figure_path = os.path.join(figure_directory, figure_name)
                plt.savefig(figure_path)
                plt.close(3)
                print(f'savefig: {figure_path}')
                elapsed_time = time.time() - start_time_savefig
                print(f'took {elapsed_time:.3f} sec')

            if EXCEL:
                print('\nclosing Excel file...')
                start_time_closingwriter = time.time()
                writer.close()
                elapsed_time = time.time() - start_time_closingwriter
                print(f'took {elapsed_time:.3f} sec')

                elapsed_time = time.time() - start_time_region
                print(f'\nThis region took {elapsed_time:.3f} sec in total')

        if EXCEL_SUMARRY:
            df_summary.to_excel(writer_summary, sheet_name='summary')
            writer_summary.close()
            elapsed_time = time.time() - start_time_overall
            print(f'\nThis process took {elapsed_time:.3f} sec in total')

    print('\nCongrats!! This process successfully DONE!!  d(^.^)b')


if __name__ == '__main__':
    main()
