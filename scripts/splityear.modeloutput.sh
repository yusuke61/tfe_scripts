#!/bin/bash

#=== Settings ======================================================================
SPLITYEAR=ON
#SPLITYEAR=OFF  # to check this script...


period=$1;          scns=$2;         socs=$3;       co2=$4
#period=historical; scns=historical; socs=histsoc;  co2=co2
#period=historical; scns=picontrol;  socs=histsoc;  co2=co2
#period=future;     scns=rcp26;      socs=2005soc;  co2=co2
#period=future;     scns=rcp85;      socs=2005soc;  co2=co2
#period=future;     scns=picontrol;  socs=2005soc;  co2=co2
#period=future;     scns=rcp26;      socs=rcp26soc; co2=co2
#period=future;     scns=rcp85;      socs=rcp85soc; co2=co2
#period=future;     scns=rcp26;      socs=2005soc;  co2=2005co2  # for co2 experiment
#period=future;     scns=rcp85;      socs=2005soc;  co2=2005co2  # for co2 experiment

if [ ${co2} = 2005co2 ]; then
  #ghms="lpjml matsiro"
  ghms="matsiro"
elif [ ${socs} = rcp26soc ]; then
  #ghms="cwatm h08 lpjml matsiro"  # selected ghms for the analysis
  ghms="matsiro"  # selected ghms for the analysis
else
  #ghms="clm45 clm50 cwatm dbh h08 jules_w1 lpjml matsiro mpi-hm orchidee pcr-globwb watergap2"  # full
  #ghms="cwatm h08 lpjml matsiro watergap2"  # selected ghms for the analysis
  #ghms="cwatm"
  #ghms="h08"
  #ghms="lpjml"
  ghms="matsiro"
  #ghms="watergap2"
fi

#gcms="hadgem2-es ipsl-cm5a-lr gfdl-esm2m miroc5"
gcms="gfdl-esm2m"
 
dt='daily'

variables="dis"

inDir=/data/rg001/sgec0017/data/isimip2b/out/nc/water_global
outDir=/data/rg001/sgec0017/data/isimip2b/out/nc/water_global.processed


#=== Main ============================================================================
for scn in $scns ; do
  for soc in $socs; do
    for ghm in $ghms; do
      for gcm in $gcms; do
        for variable in $variables; do
  
          if [ ${dt} = daily ]; then
            if [ ${period} = 'pre-industrial' ]; then
              PRDs="1661_1670 1671_1680 1681_1690 1691_1700 1701_1710 1711_1720 1721_1730 1731_1740 1741_1750 1751_1760 1761_1770 1771_1780 1781_1790 1791_1800 1801_1810 1811_1820 1821_1830 1831_1840 1841_1850 1851_1860"
            elif [ ${period} = historical ]; then
              PRDs="1861_1870 1871_1880 1881_1890 1891_1900 1901_1910 1911_1920 1921_1930 1931_1940 1941_1950 1951_1960 1961_1970 1971_1980 1981_1990 1991_2000 2001_2005"
            elif [ ${period} = future ]; then
              PRDs="2006_2010 2011_2020 2021_2030 2031_2040 2041_2050 2051_2060 2061_2070 2071_2080 2081_2090 2091_2099"
            else
              echo something is wrong. Check dt or scn!!
            fi
          elif [ ${dt} = monthly ]; then
            if [ ${period} = 'pre-industrial' ]; then
              PRDs="1661_1860"
            elif [ ${period} = historical ]; then
              PRDs="1861_2005"
            elif [ ${period} = future ]; then
              PEDs="2006_2099" 
            else
              echo something is wrong. Check dt or scn!!; exit
            fi

          fi

          
          indir=${inDir}/${ghm}/${gcm}/${period}
          outdir=${outDir}/${ghm}/${gcm}/${scn}
          if [ ! -d ${outdir} ]; then
            mkdir ${outdir}
          fi
 
          _soc=${soc}
          if [ ${period} = historical -a ${scns} = picontrol -a ${ghm} = h08 ]; then
            _soc=2005soc
          fi

          title=${ghm}_${gcm}_ewembi_${scn}_${_soc}_${co2}_${variable}_global_${dt}_

          for PRD in $PRDs ; do
  
            if   [ -f ${indir}/${title}${PRD}.nc4 ]; then
              infile=${indir}/${title}${PRD}.nc4
            elif [ -L ${indir}/${title}${PRD}.nc4 ]; then  # symbolic link
              infile=${indir}/${title}${PRD}.nc4
            else
                echo '\n!!!!!! infile is NOT found !!!!!!'
                ls -l ${indir}/${title}${PRD}.nc4
                echo 'stop processing...'
                exit
            fi

            ls -l ${infile}
            if [ ${SPLITYEAR} = ON ]; then
              cdo splityear ${infile} ${outdir}/${title}
            fi
  
          done # PRD
  
        done # variable
      done # gcm
    done # ghm
  done # soc
done # scn

echo "----- This is the end of the process :)"

exit
