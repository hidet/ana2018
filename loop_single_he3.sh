#!/bin/bash

# He3 runs
#run=(160 161)
run=(162     163     165     166     167     170     171     172
     175     176     177     178     179     180     181     182     183     186
     187     188     189     190     191     192     193     194     195     196
     197     198     201     202     203     204     205     206     207     208
     209     210     212     215     216     217     218     222     223     224
     225     226     227     228     231     232     233     239     240     241
     242     259     260     261     263     264     266     267     268     269
     270     271     272     275     276     277     278     279     281     282
     283     284     286     287     290     291     292     293     294     295
     296     297     298     299     300     301)

## He4 runs
#run=(320     321     327     328     329     330     331     332     333     334
#     335     336     337     338     340     341     345     346     347     348
#     349     350     351     352     353     354     355     356     999     359
#     363     364     366     367     368     369     370     371     372     373
#     374     375     376     377     381     382     383     384     385     386
#     387     389     395     396     397     398     399     400     401     402
#     403     406     407     408     410     411     421     422     423     424)



for r in "${run[@]}"
do
    ## beamON, forceNew=true, summaryNew=true, calibNew=true, exttrigNew=true, grptrigNew=true, dump ROOT=true, ROOTEXT=true, ROOTGRP=true
    #python run_heates_single.py $r -fsceg -REG --beam=on --sprmc=on --jbrsc=on --pre=150 --post=550
    
    ## beamOFF, forceNew=true, summaryNew=true, calibNew=true, exttrigNew=true, grptrigNew=true, dump ROOT=true, ROOTEXT=true, ROOTGRP=true
    #python run_heates_single.py $r -fsceg -REG --beam=off --sprmc=on --jbrsc=on --pre=150 --post=550
    
    ## beamON, new calibration
    python run_heates_single.py $r -c -REG --beam=on --sprmc=on --jbrsc=on --pre=150 --post=550
    
    ## beamOFF, new calibration
    #python run_heates_single.py $r -c -REG --beam=off --sprmc=on --jbrsc=on --pre=150 --post=550
    
    ## beamON, new external trigger and new calib
    #python run_heates_single.py $r -ce -REG --beam=on --sprmc=on --jbrsc=on --pre=150 --post=550
    
    ## beamOn, new group trigger and new calib
    #python run_heates_single.py $r -cg -REG --beam=on --sprmc=on --jbrsc=on --pre=150 --post=550
    
    echo $r
done
