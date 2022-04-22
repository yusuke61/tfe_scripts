#!/home/gec/sgec0017/local/anaconda3/bin/python

import os,sys
import numpy as np
from itertools import product
from datetime import datetime
from utiltools import fillinMasked, extractMasked
from numpy import array, ma, zeros, concatenate, arange
from netCDF4 import Dataset


###!!! Edit Here !!!#######################################################################
#ghms = ['cwatm', 'h08', 'lpjml', 'matsiro', 'watergap2']
#ghms = ['cwatm']
#ghms = ['h08']
#ghms = ['lpjml']
#ghms = ['matsiro']
#ghms = ['watergap2']
ghms = [sys.argv[1]]

#gcms = ['hadgem2-es','ipsl-cm5a-lr','gfdl-esm2m','miroc5']
#gcms = ['hadgem2-es']
#gcms = ['ipsl-cm5a-lr']
#gcms = ['gfdl-esm2m']
#gcms = ['miroc5']
gcms = [sys.argv[2]]

scenarios = ['historical']
socs = ['histsoc']
co2s = ['co2']

#syear, eyear = 1861, 2005  # main
#syear, eyear = 1971, 2005
syear, eyear = 1979, 2013  # test


#wins = [15, 10,  7,  5]
#thrshs = [90, 85, 80, 75, 70, 65]
wins = [15]
thrshs = [80]
#thrshs = [90]
#thrshs = [80, 90]


###########################################################################################
datDIR     = '/data/rg001/sgec0017/data/isimip2b/out/nc/water_global.processed'
outDIR     = '/data/rg001/sgec0017/data/isimip2b.drought'
lndmskPath = '/data/rg001/sgec0017/data/mapmask/ISIMIP2b_landseamask_generic.nc4'


#------------------------------------------------------------------------------------------
def draw_MAP(aSrc, ghm, gcm, scn, soc, year, lndmask):

    print('Drawing map: %s_%s_%s_%s_%i'.format(ghm, gcm, scn, soc, year))
    minVal = aSrc.min()
    maxVal = aSrc.max()

    figure()
    subplots_adjust(left=0.02, bottom=0.02, right=0.98, top=0.95)
    title('dis [m3/s]: {}_{}_{}_{}_{} ({:.2}-{:.2})'.format(ghm, gcm, scn, soc, year, minVal, maxVal))
    imshow(aSrc.mean(0))
    tick_params(labelbottom='off', labelleft='off')
    colorbar(orientation='horizontal', pad=0.001, aspect=35)
    figName = 'dis_{}_{}_{}_{}_{}.png'.format(ghm, gcm, scn, soc, year)
    figDir  = os.path.join(outDIR, 'checkFigBox')
    if not os.path.isdir(figDir): os.makedirs(figDir)
    figPath = os.path.join(figDir, figName)
    savefig(figPath)
    print(figPath)
    close()


#------------------------------------------------------------------------------------------
def getINDXs(thrsh, win, syear, eyear):
    #>>>> get INDX for thrsh
    nYear = eyear - syear + 1
    nData = (win * 2 + 1) * nYear
    indx = int(nData * (1-thrsh/100.))
    return indx


#------------------------------------------------------------------------------------------
def read_NC(scn, soc, ghm, gcm, year, co2, mask):

    #>>>> Read only Land Data
    if year > 2005: scn, soc = 'rcp85', '2005soc'
    ncFileName = '{}_{}_ewembi_{}_{}_{}_dis_global_daily_{}.nc'.format(ghm, gcm, scn, soc, co2, year)
    srcPath = os.path.join(datDIR, ghm, gcm, scn, ncFileName)
    if not os.path.isfile(srcPath): 
        srcPath = srcPath[:-2]+'nc4'

    if not os.path.isfile(srcPath): 
        print('Error!! {} is not exist... Check!!'.format(srcPath))
        sys.exit()
    else:
        ncFile = Dataset(srcPath)
        aSrc = ncFile.variables['dis'][:]
        shp = aSrc.shape
        if shp[0] == 366:
            aSrc = concatenate((aSrc[:59], aSrc[60:]), axis=0)
            print('Reading... {} {}'.format(srcPath, shp))  #, type(aSrc)
            print('This is a leap year. Skipped the leap day. {} >> {}'.format(shp, aSrc.shape))
        else:
            print('Reading... {} {}'.format(srcPath, shp))  #, type(aSrc)

        #if year == eyear: draw_MAP(aSrc, ghm, gcm, scn, soc, year, mask)
        return extractMasked(aSrc, mask)  # Extract only land  (nday, nland)


#------------------------------------------------------------------------------------------
def getQvalue(gridData, indx):
    #>>>> Sort and Select
    gridData = gridData.ravel()                # ((1+win*2), nYear) -> ((1+win*2)*nYear)
    gridData.sort()
    return gridData[indx]


#------------------------------------------------------------------------------------------
def write_NC(aSrc, ghm, gcm, scn, soc, qvalType, win):
    print('writeout...')

    for ithrsh, thrsh in enumerate(thrshs):

        outFile = 'Q{:02}win{:02}.nc4'.format(thrsh, win)
        outDir = os.path.join(outDIR, ghm, gcm, '{}.{}_{}'.format(qvalType, syear, eyear), 'Qvalues')
        if not os.path.isdir(outDir): os.makedirs(outDir)
        outPath = os.path.join(outDir, outFile)

        rootgrp = Dataset(outPath, 'w', format='NETCDF4')
        import time
        rootgrp.title = 'Q{} for drought analysis'.format(win)
        rootgrp.description = 'Qvalue for drought analysis with win={} and thrsh={} ({}-{}): {} faced by {}'.format(win,thrsh,syear,eyear,ghm,gcm)
        rootgrp.history = 'Created ' + time.ctime(time.time())
        rootgrp.source = 'ISI-MIP FT: {}_{}_{}'.format(ghm,gcm,scn)
        rootgrp.institution = 'NIES'
        rootgrp.contact = 'satoh.yusuke@nies.go.jp'

        lon = rootgrp.createDimension('lon', size=720)
        lat = rootgrp.createDimension('lat', size=360)
        time = rootgrp.createDimension('time', size=None)

        longitudes = rootgrp.createVariable('lon', 'f8', ('lon',))
        longitudes.long_name = 'longitude'
        longitudes.units = 'degrees east'
        longitudes.standard_name = 'longitude'
        longitudes.axis = 'X'

        latitudes = rootgrp.createVariable('lat', 'f8', ('lat',))
        latitudes.long_name = 'latitude'
        latitudes.units = 'degrees north'
        latitudes.standard_name = 'latitude'
        latitudes.axis = 'Y'

        times = rootgrp.createVariable('time', 'f8', ('time',))
        times.units = 'days since 1901-01-01 00:00:00.0'
        times.calendar = 'gregorian'

        srcs = rootgrp.createVariable('Qvalue', 'f4', ('time', 'lat', 'lon'), 
                                      zlib=True, complevel=5, fill_value=np.float32(1e+20), 
                                      chunksizes=(1, 360, 720))
        srcs.long_name = 'Q{:02}win{:02} based on {}-{}'.format(thrsh, win, syear, eyear) # longName
        srcs.standard_name = 'Q{:02}win{:02}'.format(thrsh, win) # standard_name
        srcs.missing_value = np.float32(1e+20)
        srcs.unit = 'm3/s'  # unitName
        if eyear > 2005:
            srcs.memo = 'rcp85 was referred during 2006-2013 because of the ISIMIP2b bias-correction method.'

        #=== Allocate data ===
        stime = (datetime(eyear, 1, 1)-datetime(eyear, 1, 1)).days
        times[:] = [stime + n for n in range(365)]
        latitudes[:] = arange(89.75, -90, -0.5)
        longitudes[:] = arange(-179.75, 180, 0.5)
        srcs[:] = aSrc[ithrsh]

        #=== close ===
        rootgrp.close()
        print('save: {}'.format(outPath))


#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
def main(*args):
    print(args)
    strTIME = datetime.now()
    print('Start mkQvalue.py  @', strTIME.strftime("%Y-%m-%d %H:%M:%S"), '\n')

    lndmask = ma.masked_not_equal(Dataset(lndmskPath).variables['LSM'][:][0], 1.0).mask

    for co2, scn, soc, ghm, gcm, win in product(co2s, scenarios, socs, ghms, gcms, wins):
        strTime = datetime.now()
        print('Process : {}, {}, {}, {}, win:{}, thrshs:{},  {}-{}'.format(scn, soc, ghm, gcm, win, ','.join([str(i) for i in thrshs]), syear, eyear))
        print('Start  @', strTime.strftime("%Y-%m-%d %H:%M:%S"))

        qvalType = 'Qvalue_{}_{}_{}'.format(scn, soc, co2)

        INDXs = [getINDXs(thrsh, win, syear, eyear) for thrsh in thrshs]

        aSRC = array([read_NC(scn, soc, ghm, gcm, year, co2, lndmask) for year in range(syear, eyear+1)])  # (nYear, 365, 67420)
        print(aSRC.shape, '<< %.2f GB\n'%(aSRC.shape[0]*aSRC.shape[1]*aSRC.shape[2]*4*1e-9))

        print('Processing addHF...')
        aHead = concatenate((aSRC[-1, -win:, ...][None, ...], aSRC[:-1, -win:, ...]), axis=0)
        aFoot = concatenate((aSRC[1:, :win, ...], aSRC[0, :win, ...][None, ...]), axis=0)
        aSRC = concatenate((aHead, aSRC, aFoot), axis=1)
        del aHead, aFoot

        print('Get Qvalue...')
        aQval = array([[[getQvalue(aSrc[idy-win:idy+win+1, :], indx) for idy in range(win, win+365)]
                                                                     for aSrc in aSRC.T]
                                                                     for indx in INDXs])                # (nthrsh, 67420, 365)
        del aSRC

        aQval = aQval.transpose((0, 2, 1))                                                              # (nthrsh, 365, 67420)
        print('fillinMasked...')
        aQval = ma.masked_equal(fillinMasked(aQval, lndmask), -999).filled(1e+20)                       # (nthrsh, 365, 360, 720)

        write_NC(aQval, ghm, gcm, scn, soc, qvalType, win)
        del aQval

        endTime  = datetime.now()
        diffTime = endTime - strTime
        print('End    @', endTime.strftime("%Y-%m-%d %H:%M:%S"))
        print('took {} min. '.format(int(diffTime.seconds/60)))

    endTIME = datetime.now()
    diffTIME = endTIME - strTIME
    print("\nJob complated !!! d(^.^)")
    print('took {} min in total.'.format(int(diffTIME.seconds/60)))


if __name__=='__main__':
    main(*sys.argv)

