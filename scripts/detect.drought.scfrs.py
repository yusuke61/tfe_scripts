#! /usr/bin/python
#--------------------------------------------------------------------
import os,sys,re
import itertools
from copy import copy
from numpy import array, concatenate, arange, ma, mean
from datetime import datetime
from netCDF4 import Dataset
from utiltools import fillinMasked, extractMasked
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Setting >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
## for win, Len, tau, Q, scn, soc, ghm, gcm

#ghms  = ['cwatm', 'h08', 'lpjml','matsiro', 'watergap2']
#ghms  = ['cwatm']
#ghms  = ['h08']
#ghms  = ['matsiro']
#ghms  = ['lpjml']
#ghms  = ['watergap2']
ghms = [sys.argv[1]]

#gcms = ['hadgem2-es','ipsl-cm5a-lr','gfdl-esm2m','miroc5']
#gcms = ['hadgem2-es']
#gcms = ['ipsl-cm5a-lr']
#gcms = ['gfdl-esm2m']
#gcms = ['miroc5']
gcms = [sys.argv[2]]

# target simlation type (future)
sims  = [
          (sys.argv[3], sys.argv[4], sys.argv[5]),
          #('picontrol', '2005soc',  'co2'),
          #('rcp26',     '2005soc',  'co2'),
          #('rcp85',     '2005soc',  'co2'),
          #('rcp26',     'rcp26soc', 'co2'),
          #('rcp26',     '2005soc',  '2005co2'),
          #('rcp85',     '2005soc',  '2005co2'),
          #('rcp26',     'nosoc',    'co2'),
          #('rcp85',     'nosoc',    'co2'),
          ]

#SEASONs = ['ALL', 'DJF', 'MAM', 'JJA', 'SON']      # all
#SEASONs = ['ALL']                                  # annual
#SEASONs = ['DJF', 'MAM', 'JJA', 'SON']              # seasonal
#SEASONs = ['DRY', 'WET']              # seasonal
#SEASONs = ['WET']              # seasonal
SEASONs = [sys.argv[6]]              # seasonal

# target period of the analysis
syear, eyear = 1861, 2099  # main
#syear, eyear = 1971, 2099
#syear, eyear = 1979, 2099  # test

# reference period of drought threshold
ref_syear, ref_eyear = 1861, 2005  # main
#ref_syear, ref_eyear = 1971, 2005
#ref_syear, ref_eyear = 1979, 2013  # test

base_year = 1901

#--- qvalType
autQval = False                  # Default
#autQval = True                   # automatically specify, depending on soc type: 'Qvale_hist_%s'%soc

#--- Qvalue
#Qs = [90, 85, 80, 75, 70, 65]
#Qs = [90,80,70]
#Qs = [80,90,70]
Qs = [80]

#--- window size for Qvalue
#WINs  = [15, 10, 7, 5]
WINs = [15]

#--- minimum drought days
#LENs = [180,90,60,30,14,7]
#LENs = [30,14,60,90,180,7]
#LENs = [1, 30, 60, 90]
LENs = [30]
#LENs = [60, 90]

#--- Pooled duration
TAUs = [4]




#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>#
datDIR     = '/data/rg001/sgec0017/data/isimip2b/out/nc/water_global.processed'
drghtDIR   = '/data/rg001/sgec0017/data/isimip2b.drought'
lndmskPath = '/data/rg001/sgec0017/data/mapmask/ISIMIP2b_landseamask_generic.nc4'

dictIndex = {
    'nDayTot': 'days',
    'dfcTot':  'm3',
    'nEvent':  'times',
    'nOnset':  'times',
    'avrDSL':  'days',
    'maxDSL':  'days',
    }

dictSeason = {
    'ALL': [[[  0,366]       ], 365],
    'DJF': [[[334,366],[0,59]],  90],
    'MAM': [[[ 59,151]       ],  92],
    'JJA': [[[151,243]       ],  92],
    'SON': [[[243,334]       ],  91],
    'DRY': [[                ],  91],
    'WET': [[                ],  91],
    }

years = range(syear, eyear+1)
nyear = len(years)


#------------------------------------------------------------------------------------------
def read_NC(ghm, gcm, scn, soc, co2, iyear, mask):
    #>>>> Read only Land Data

    if iyear <= 2005:  # historical period
        # scn
        if scn == 'picontrol':
            pass
        else:
            scn = 'historical'
        # soc
        if soc == 'nosoc':
            pass
        else:
            if scn == 'picontrol' and ghm == 'h08':
                soc = '2005soc'
            else:
                soc = 'histsoc'
        # co2
        if co2 == '2005co2':
            co2 = 'co2'
    if iyear == syear - 1: iyear = syear
    if iyear == eyear + 1: iyear = eyear
    if iyear > 2099:
        print('Are you sure? {}?'.format(iyear))

    ncFileName = '{}_{}_ewembi_{}_{}_{}_dis_global_daily_{}.nc'.format(ghm, gcm, scn, soc, co2, iyear)
    srcPath = os.path.join(datDIR, ghm, gcm, scn, ncFileName)
    #if not os.path.isfile(srcPath): 
    #    srcPath = srcPath[:-2]+'nc4'

    if not os.path.isfile(srcPath): 
        print('Error!! {} is not exist... Check!!'.format(srcPath))
        sys.exit()
    else:
        ncFile = Dataset(srcPath)
        aSrc = ncFile.variables['dis'][:]
        shp = aSrc.shape
        if shp[0] == 366:  # leap year
            aSrc = concatenate((aSrc[:59], aSrc[60:]), axis=0)
            print('Reading... {} {}'.format(srcPath, shp))  #, type(aSrc)
            print('This is a leap year. Skipped the leap day. {} >> {}'.format(shp, aSrc.shape))
        else:
            print('Reading... {} {}'.format(srcPath, shp))  #, type(aSrc)
        aSrc = extractMasked(aSrc,mask)
        #print('extractMasked {}  {}-{}'.format(aSrc.shape, aSrc.min(), aSrc.max()))
        #if iyear == eyear: draw_MAP(aSrc, ghm, gcm, scn, soc, iyear, mask)

        return aSrc


#------------------------------------------------------------------------------------------------------
def byte2str(src):

    #print('before {}'.format(type(src)))
    #print('before: {}'.format(src))
    #src[where(src==b'T')] = 'T'
    #src[where(src==b'F')] = 'F'

    nlnd, ndy = 67420, 365
    src = [['F' if src[ilnd, idy] == b'F' else 'T' for idy in range(ndy)] for ilnd in range(nlnd)]
    src = array(src)

    #print('after: {}'.format(type(src)))
    #print('after: {} {}'.format(src, array(src).shape))

    return src
    


#------------------------------------------------------------------------------------------------------
def droughtJudge(FLUG,FLUGprv,FLUGnxt, Len, tau, pF, pT, pTF):
    """
    Judge drought days at a grid.
    Retrun a list of bool (True or False). (False: drought)
    """
    #print('\ncheck3 >> FLUG:\n{}\n'.format(FLUG))
    
    # conv from byte type to string type 
    #FLUG = ['F' if flug == b'F' else 'T' for flug in FLUG]

    headFlug = False
    tailFlug = False
    #if not 'F' in FLUG:
    if not 'F' in FLUG:
        return [True for i in FLUG]

    else:
        ## < 1. the Len process: removing shorter droughts >           (ex)  "FLUG" is a list of 'F' and 'T'.
        ft = ''.join(FLUG)                                            # 'TTTTTFTFFTFTFTFFFFTTTFFFFTFTFTTFFTTFFTTTTFFFFFTTTFTTFFFFFTFTFTFFFFTFTFTFFFTFTTFFTTFFFFFTFT' 
        Fs = pF.findall(ft)                                           # ['F', 'FF', 'F', 'F', 'FFFF', 'FFFF', 'F', 'F', 'FF', 'FF', 'FFFFF', 'F', 'FFFFF', 'F', 'F', 'FFFF', 'F', 'F', 'FFF', 'F', 'FF', 'FFFFF', 'F']
        Ts = pT.findall(ft)                                           # ['TTTTT', 'T', 'T', 'T', 'T', 'TTT', 'T', 'T', 'TT', 'TT', 'TTTT', 'TTT', 'TT', 'T', 'T', 'T', 'T', 'T', 'T', 'T', 'TT', 'TT', 'T', 'T']
        first = ft[0]                                                 # 'T'

        ##--- treat the first "F" and the last "F"
        if FLUGprv[-1] == 'F' and FLUG[0] == 'F':
            headFlug = True
            f_prv = pF.findall(''.join(FLUGprv))[-1]                  # 'FFFFF'  <string>
            nHead = len(f_prv)                                        # 5
            Fs[0] = f_prv + Fs[0]                                     # 
        if FLUG[-1] == 'F' and FLUGnxt[0] == 'F':
            tailFlug = True
            f_nxt  = pF.findall(''.join(FLUGnxt))[0]                  # 'FFF'    <string>
            nTail  = len(f_nxt)                                       # 3
            Fs[-1] = Fs[-1] + f_nxt                                   # 

        ##--- remove shorter droughts
        for i in range(len(Fs)):                                      # ex) Len = 3  >>  'FF' and 'F' are removed
            if len(Fs[i]) < Len: Fs[i] = 'T'*len(Fs[i])               # ['T', 'TT', 'T', 'T', 'FFFF', 'FFFF', 'T', 'T', 'TT', 'TT', 'FFFFF', 'T', 'FFFFF', 'T', 'T', 'FFFF', 'T', 'T', 'FFF', 'T', 'TT', 'FFFFF', 'T']

        ##--- marge
        ft = []
        if len(Fs) > len(Ts):  # originally starts from F and ends with F
            Ts.append('')
            for i in range(len(Fs)):
                ft.append(Fs[i])
                ft.append(Ts[i])
        elif len(Fs) < len(Ts): # originally starts from T and ends with T
            Fs.append('')
            for i in range(len(Ts)):
                ft.append(Ts[i])
                ft.append(Fs[i])
        else:                   # originally starts from T or F and ends with another
            if first == 'F':
                for i in range(len(Fs)):
                    ft.append(Fs[i])
                    ft.append(Ts[i])
            elif first == 'T':
                for i in range(len(Ts)):
                    ft.append(Ts[i])
                    ft.append(Fs[i])
        ft = ''.join(ft)                                            # 'TTTTTTTTTTTTTTFFFFTTTFFFFTTTTTTTTTTTTTTTTFFFFFTTTTTTFFFFFTTTTTFFFFTTTTTFFFTTTTTTTTFFFFFTTT'
        if headFlug: ft = ft[nHead:]
        if tailFlug: ft = ft[:-nTail]
        if not len(ft) == 365: print('len(ft) is not 365. Check the process!')
        headFlug = False
        tailFlug = False

        ## < 2. the Tau process: interporating short breaks between droughts >
        Fs = pF.findall(ft)                                          # ['FFFF', 'FFFF', 'FFFFF', 'FFFFF', 'FFFF', 'FFF', 'FFFFF']
        Ts = pT.findall(ft)                                          # ['TTTTTTTTTTTTTT', 'TTT', 'TTTTTTTTTTTTTTTT', 'TTTTTT', 'TTTTT', 'TTTTT', 'TTTTTTTT', 'TTT']
        first = ft[0]                                                # 'T'

        Ts = pT.findall(ft)                                          # ['TTTTTTTTTTTTTT', 'TTT', 'TTTTTTTTTTTTTTTT', 'TTTTTT', 'TTTTT', 'TTTTT', 'TTTTTTTT', 'TTT']
        first = ft[0]                                                # 'T'

        ##--- treat the first T and the last T
        if FLUGprv[-1] == 'T' and ft[0] == 'T':
            headFlug = True
            t_prv = pT.findall(''.join(FLUGprv))[-1]
            nHead = len(t_prv)
            Ts[0] = t_prv + Ts[0]
        if ft[-1] == 'T' and FLUGnxt[0] == 'T':
            tailFlug = True
            t_nxt = pT.findall(''.join(FLUGnxt))[0]
            nTail = len(t_nxt)
            Ts[-1] = Ts[-1] + t_nxt

        ##--- remove shorter breaks
        for i in range(len(Ts)):                                    # ex) tau = 4
            if len(Ts[i]) < tau: Ts[i] = 'F'*len(Ts[i])             # ['TTTTTTTTTTTTTT', 'FFF', 'TTTTTTTTTTTTTTTT', 'TTTTTT', 'TTTTT', 'TTTTT', 'TTTTTTTT', 'FFF']

        ##--- marge
        ft = []
        if   len(Fs) > len(Ts):  # originally starts from F and ends with F
            Ts.append('')
            for i in range(len(Fs)):
                ft.append(Fs[i])
                ft.append(Ts[i])
        elif len(Fs) < len(Ts): # originally starts from T and ends with T
            Fs.append('')
            for i in range(len(Ts)):
                ft.append(Ts[i])
                ft.append(Fs[i])
        else:                   # originally starts from T or F and ends with another
            if   first == 'F':
                for i in range(len(Fs)):
                    ft.append(Fs[i])
                    ft.append(Ts[i])
            elif first == 'T':
                for i in range(len(Ts)):
                    ft.append(Ts[i])
                    ft.append(Fs[i])
        ft = ''.join(ft)                                            # 'TTTTTTTTTTTTTTFFFFFFFFFFFTTTTTTTTTTTTTTTTFFFFFTTTTTTFFFFFTTTTTFFFFTTTTTFFFTTTTTTTTFFFFFFFF'
        if headFlug: ft = ft[nHead:]
        if tailFlug: ft = ft[:-nTail]

        if not len(ft) == 365: 
            print('len(ft) is not 365. Check the process!')

        flug_updated= [False if i=='F' else True for i in list(ft)]
        #print('flug_updated: {}\n'.format(flug_updated))
        #sys.exit()
        #return flug_updated
        return [False if i=='F' else True for i in list(ft)]


#------------------------------------------------------------------------------------------------------
def get_droughtValues(flug, dfc, season, pF, pTF):

    if not ('T' in flug or 'F' in flug):
        print('flug:\n{}'.format(flug))
        print('Terminate this process...')
        sys.exit()

    ##--- Lists of Drought Length
    ft = ''.join(flug)
    Fs = pF.findall(ft)
    DLen = [len(g) for g in Fs]
    
    ## < 3. Get drought values >
    """
     nDayTot [days] : The number of drought days in the year or season.
     dfcTot  [m3]   : Total deficit volume.
     nEvent  [times]: The number of drought event in the year or season. (Multi-year drought could split into yearly)
     nOnset  [times]: The number of onset of derought in the year or season.
     avrDSL  [days] : Average dry spell length of the year or season. 
     maxDSL  [days] : Maximum dry spell lenght of the year or season. 
    """

    nDayTot = sum(DLen)
    nEvent  = len(DLen)
    avrDSL  = mean(DLen)
    nOnset  = len(pTF.findall(ft))
    dfcTot  = dfc.sum() * 86400 * dictSeason[season][1]   # unit conv  [m3/s] >> [m3] per season
    if len(DLen) != 0: maxDSL = max(DLen)
    else             : maxDSL = 0

    if nDayTot > 365: print('Warning Stop : nDayTot > 365...'); sys.exit()

    if season == 'ALL':
        return nDayTot, dfcTot, nEvent, nOnset, avrDSL, maxDSL
    else:
        return nDayTot, dfcTot,  1e+20, nOnset,  1e+20,  1e+20


#------------------------------------------------------------------------------------------
def write_NC(aSrc, var, unit, soc, scn, qvalType, win, Q, Len, tau, ghm, gcm, season, outDir):

    outFile = 'Q{:02}win{:02}_Len{:03}tau{}_{}_{}.nc4'.format(Q,win,Len,tau,season,var)
    outPath = os.path.join(outDir, outFile)

    rootgrp = Dataset(outPath,'w', format='NETCDF4')

    import time
    rootgrp.title = '{} for {} in {}_{}'.format(var,season,scn,soc)
    rootgrp.description = '{} of droughts in {} against win{}&Q{} for events with Len{}&tau{}. \
                           {} faced by {} for {} ({}-{})'\
                          .format(var,season,win,Q,Len,tau,ghm,gcm,scn,syear,eyear)
    rootgrp.history     = 'Created ' + time.ctime(time.time())
    rootgrp.source      = 'ISI-MIP2b: {}_{}_{}_{}'.format(ghm,gcm,scn,soc)
    rootgrp.institution = 'NIES'
    rootgrp.contact     = 'satoh.yusuke@nies.go.jp'

    lon  = rootgrp.createDimension('lon', size=720)
    lat  = rootgrp.createDimension('lat', size=360)
    time = rootgrp.createDimension('time', size=None)

    longitudes               = rootgrp.createVariable('lon', 'f8', ('lon',))
    longitudes.long_name     = 'longitude'
    longitudes.units         = 'degrees east'
    longitudes.standard_name = 'longitude'
    longitudes.axis          = 'X'

    latitudes                = rootgrp.createVariable('lat', 'f8', ('lat',))
    latitudes.long_name      = 'latitude'
    latitudes.units          = 'degrees north'
    latitudes.standard_name  = 'latitude'
    latitudes.axis           = 'Y'

    times                    = rootgrp.createVariable('time', 'f8', ('time',))
    times.units              = 'years since {}-01-01 00:00:00.0'.format(base_year)
    times.calendar           = 'gregorian'

    srcs                     = rootgrp.createVariable(var, 'f4', ('time','lat','lon'),
                                                      fill_value=1e+20, chunksizes=(1, 360, 720), zlib=True)
    srcs.long_name           = '{} in {} against Q{:02}win{:02} based on {}-{} with Len{} and tau{}'.format(
                                var, season, Q, win, ref_syear, ref_eyear, Len, tau)
    srcs.standard_name       = 'Q{:02}win{:02}_Len{:03}tau{}_{}_{}'.format(Q,win,Len,tau,season,var) 
    srcs.units               = unit 

    #=== Allocate data ===
    stime         = syear - 1901
    times[:]      = [stime + n for n in range(nyear)]
    latitudes[:]  = arange(89.75, -90, -0.5)
    longitudes[:] = arange(-179.75, 180, 0.5)
    srcs[:]       = aSrc

    #=== close ===
    rootgrp.close()
    print('write_NC: {}'.format(outPath))


#------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------
def main(*args):
    strTIME = datetime.now()
    print('START detect.drought.py  @{}'.format(strTIME.strftime("%Y-%m-%d %H:%M:%S"), '\n'))

    # Preparation1
    lndmask = Dataset(lndmskPath).variables['LSM'][:][0].mask
    pF      = re.compile('F+')
    pT      = re.compile('T+')
    pTF     = re.compile('TF')

    for win, Len, tau, Q, (scn, soc, co2), ghm, gcm in itertools.product(WINs, LENs, TAUs, Qs, sims, ghms, gcms):
        strTime = datetime.now()

        #if autQval: qvalType = 'Qvalue_historical_{}'.format(soc)
        #else:
        #    if soc == 'nosoc': qvalType = 'Qvalue_historical_nosoc'
        #    else:              qvalType = 'Qvalue_historical_histsoc'
        #print('autQval : {}'.format(autQval))
        #print('qvalType: {}\n'.format(qvalType))
        qvalType = 'Qvalue_historical_histsoc_co2'  # use consistent reference!!
        outDir    = os.path.join(drghtDIR, ghm, gcm, '{}.{}_{}'.format(qvalType, ref_syear, ref_eyear), 
                                 'droughtStats.{}_{}_{}.{}_{}'.format(scn, soc, co2, syear, eyear))
        if not os.path.isdir(outDir): os.makedirs(outDir)

        #checkFile = 'Q{:02}win{:02}_Len{:03}tau{}_ALL_nDayTot.nc4'.format(Q, win, Len, tau)
        #checkPath = os.path.join(outDir, checkFile)
        #if os.path.exists(checkPath):
        #    print('OK, {} already exists. Skip!\n'.format(checkPath))
        #    break

        else:  # main process
            print('\nJob started !!!  {} {} {} {} {} with Q{} win{} for Len{} tau{}'.format(ghm, gcm, scn, soc, co2, Q, win, Len, tau))
            print(strTIME.strftime("%Y-%m-%d %H:%M:%S"))

            qvalfile = 'Q{:02}win{:02}.nc4'.format(Q, win)
            qvalPath = os.path.join(drghtDIR, ghm, gcm, '{}.{}_{}'.format(qvalType, ref_syear, ref_eyear), 
                                    'Qvalues', qvalfile)
            aQval    = extractMasked(Dataset(qvalPath).variables['Qvalue'][:], lndmask).T             # (67420, 365)
            print('Read: {} {}'.format(qvalPath, aQval.T.shape))

            aOUT = []
            for i, iyear in enumerate(years):
                strtime = datetime.now()

                # Read flow Data
                if i == 0:
                    aRFLOWprv = read_NC(ghm, gcm, scn, soc, co2, iyear-1, lndmask).T   # (67420, 365)
                    aRFLOW    = read_NC(ghm, gcm, scn, soc, co2, iyear,   lndmask).T   # (67420, 365)
                    aRFLOWnxt = read_NC(ghm, gcm, scn, soc, co2, iyear+1, lndmask).T   # (67420, 365)
                else:
                    aRFLOWprv = copy(aRFLOW)                                           # (67420, 365)
                    aRFLOW    = copy(aRFLOWnxt)                                        # (67420, 365)
                    aRFLOWnxt = read_NC(ghm, gcm, scn, soc, co2, iyear+1, lndmask).T   # (67420, 365)

                # If aQval>=aRFLOW (aQval-aRFLOW>=0), then it is under drought condition
                # aQval-aRFLOW<0 is maskedout.
                aDFC = ma.masked_less(aQval - aRFLOW, 0)                               # (67420, 365) <maskedarray> [m3/s]

                ## If drought, drough flug is F.    (ex) [...,F,F,F,T,T,F,F,T,T,T,T,T,T,F,F,F,F,...]
                aDFLUG    = aDFC.mask.astype('S1')                                     # (67420, 365) list
                aDFLUGprv = ma.masked_less(aQval-aRFLOWprv, 0).mask.astype('S1')       # (67420, 365) list
                aDFLUGnxt = ma.masked_less(aQval-aRFLOWnxt, 0).mask.astype('S1')       # (67420, 365) list
                # for python3 ---
                aDFLUG    = byte2str(aDFLUG)     
                aDFLUGprv = byte2str(aDFLUGprv)     
                aDFLUGnxt = byte2str(aDFLUGnxt)     
                # ---------------
                if aDFLUG.shape    == (): aDFLUG    = array(['F' for i in range(365)]).astype('str')
                if aDFLUGprv.shape == (): aDFLUGprv = array(['F' for i in range(365)]).astype('str')
                if aDFLUGnxt.shape == (): aDFLUGnxt = array(['F' for i in range(365)]).astype('str')

                # Drought judgement
                # Updating drough flug  line132
                aMask = array([droughtJudge(flug, flugprv, flugnxt, Len, tau, pF, pT, pTF)
                                        for flug, flugprv, flugnxt in list(itertools.zip_longest(aDFLUG, aDFLUGprv, aDFLUGnxt))]) #  (67420, 365) <bool>
                del aDFLUG, aDFLUGprv, aDFLUGnxt
                del aRFLOWprv

                ##-- Update aDFC with the updated drought flug
                aDFC = ma.masked_array(aDFC, mask=aMask)        # (67420, 365)

                ##-- Get yearly or seasonal drought values
                aDFLUG = aMask.astype('S1')                     # (67420, 365) <bool> --> array of 'F'&'T'
                aDFLUG = byte2str(aDFLUG)
                del aMask
                aOut = []  # seasonal list for a year
                for season in SEASONs:
                    # extract data for a target period
                    if season == 'DRY' or season == 'WET':
                        # load sDOY
                        sDOY_directory = '/data/rg001/sgec0017/data/figure_box/drought_tpcd_isimip2b/find_sDOY_for_DRYandWET_season'
                        sDOY_file = 'startDOY_{}season_{}_{}_historical_histsoc_co2.nc4'.format(season, ghm, gcm)
                        sDOY_path = os.path.join(sDOY_directory, sDOY_file)
                        sDOY = Dataset(sDOY_path)['sDOY'][:]  # (360, 720)
                        sDOY = extractMasked(sDOY, lndmask)   # (67420)
                        DFlug, aDfc = [], []
                        # the land grid cell loop (each grid has different dry/wet season)
                        print(' landgridcell loop for {}...'.format(season))
                        for ilnd in range(aDFLUG.shape[0]): 
                            sdoy = int(sDOY[ilnd])
                            # One season is 91 days. If sDOY > 275(=365-91+1), some days in January and Feburuary are used.
                            if sdoy > 275:
                                edoy = 91-(365-sdoy+1)
                                DFlug.append(concatenate((aDFLUG[ilnd,sdoy:], aDFLUG[ilnd,:edoy+1])))
                                aDfc.append(concatenate((   aDFC[ilnd,sdoy:],   aDFC[ilnd,:edoy+1])))
                            else:
                                DFlug.append(aDFLUG[ilnd,sdoy:sdoy+91])
                                aDfc.append(aDFC[ilnd,sdoy:sdoy+91])
                        DFlug = array(DFlug)
                        aDfc = array(aDfc)
                    elif season == 'DJF':
                        index_str1, index_end1 = dictSeason[season][0][0]
                        index_str2, index_end2 = dictSeason[season][0][1]
                        DFlug = concatenate((aDFLUG[:, index_str1:index_end1+1], aDFLUG[:, index_str2:index_end2+1]), axis=1)
                        aDfc  = concatenate((aDFC[:, index_str1:index_end1+1], aDFC[:, index_str2:index_end2+1]), axis=1)
                    else:
                        index_str, index_end = dictSeason[season][0][0]
                        DFlug = aDFLUG[:, index_str:index_end+1]
                        aDfc  = aDFC[:, index_str:index_end+1]
                    print('{} {}'.format(season, DFlug.shape))
                    aOut.append(array([get_droughtValues(flug, dfc, season, pF, pTF) for flug, dfc in zip(DFlug, aDfc)]).T)  # get_droughtValues  (6,67420)
                                                                                                                             # aOut             (5,6,67420)
                aOUT.append(aOut)
                del aDFC, aDFLUG, DFlug, aDfc

                endtime  = datetime.now()
                difftime = endtime - strtime
                print('took {} sec in total.'.format(difftime.seconds))

            aOUT = array(aOUT)                                                       # aOUT (nYear,5,6,  67420)
            #aOUT = ma.masked_equal(fillinMasked(aOUT, lndmask),-999).filled(1e+20)  #      (nYear,5,6,360,720)
            print('\nwriteout...')
            output_variables = ['nDayTot', 'dfcTot', 'nEvent', 'nOnset', 'avrDSL', 'maxDSL']
            for (i, season), (j, index) in itertools.product(enumerate(SEASONs), enumerate(output_variables)):
                print('  {} {}'.format(season, index))
                write_NC(ma.masked_equal(fillinMasked(aOUT[:,i,j], lndmask),-999).filled(1e+20),                  # (nYear,67420) -> (nYear,360,720)
                         index, dictIndex[index], soc, scn, qvalType, win, Q, Len, tau, ghm, gcm, season, outDir)

            endTime  = datetime.now()
            diffTime = endTime - strTime
            #print('end @', endTime.strftime("%Y-%m-%d %H:%M:%S"))
            print('took {} min in total.\n\n'.format(int(diffTime.seconds/60)))

    endTIME  = datetime.now()
    diffTIME = endTIME - strTIME
    #print('END    @{}'.format(endTIME.strftime("%Y-%m-%d %H:%M:%S")))
    print('The whole process took {} min in total.\n\n'.format(int(diffTIME.seconds/60)))
    print('This process successfully finished!!  d(^o^)b')


if __name__=='__main__':
    main(*sys.argv)



