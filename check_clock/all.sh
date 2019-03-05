#!/bin/sh

fname=beam_tes_runs.csv
#fname=tmp.csv
while read line
do
  log=log/match_${line/ /_}.log
  echo python match_tes_beam.py $line > $log
  bsub -o $log python match_tes_beam.py $line
done < $fname