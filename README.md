# tfe_scripts
This repository archives scripts used in Satoh et al. 2022 NatCommns.

Scripts and their functions are listed below. 
1. splityear.modeloutput.sh			              <-- splits 10-year chunk of daily river discharge data into the annual chunk.
2. mkQvalue.scfrs.py   			                  <-- generates drought detection threshold
3. detect.drought.scfrs.py   			            <-- detects drought conditions and estimates the number of drought days per year
4. TPCD.estimator.scns.final.py   		        <-- estimates the time-series of regional average drought frequency
5. estimate_effective_non-corr_time.py   	    <-- finds available block size for the bootstrap resampling and makes SupFig 23&24
6. bootstrap.py   				                    <-- performs blocksize bootstrap resampling to generate a larger ensemble member dataset
7. basicMaps.RCPxFuture.1Dcolorbar.py   	    <-- makes Fig1a, SupFig4, and SupFig6-9
8. draw_plot_regional_avrFDD_timeseries.py   	<-- makes Fig1b and SupFig5
9. draw_tfe_map_from_bootstrap.py   		      <-- makes Fig2a, and SupFig15&16
10. draw_CDF_from_bootstrap.py   		          <-- makes Fig2b, SupFig11, and SupFig15&16
11. draw_accumulated_UDyears.py   		        <-- makes Fig3
12. draw_dgmt_map.py   			                  <-- makes Fig4a
13. draw_dgmt_CDF_from_bootstrap.py   		    <-- makes Fig4b and SupFig17
14. find_sDOY_for_DRYandWET_season.maps.py   	<-- defines dry and wet seasons and makes SupFig2
15. basicMaps.RCPxFuture_PER.py   		        <-- makes SupFig10
16. draw_heatmap.py   			                  <-- makes SupFig12, 18, and 19
17. draw_archiver_TPCDs.py   			            <-- makes SupFig13&14
18. estimate_average_Qvalue_for_reviewer3.py  <-- makes SupFig20
19. find_the_best_fittingcurve.v2.py   		     <-- makes SupFig21&22

Note
- these scripts are not completely cleaned up so that a part of the scripts contain processes, such as methodological or parameter sensitivity tests.
- Several script names and codes include the term, TPCD. TFE (Time of the First Emergence of unprecedented regional drought condition) was called TPCD (Timing of Perception Change regarding Drought) in the early stage of the study.

