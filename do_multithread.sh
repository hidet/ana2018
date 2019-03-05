#!/bin/sh 
# version 1.0, 2018/07/11, S.Yamada


echo "start run" 

OUTPUT_FILE="for_multi.txt"
maxnproc=24; 
i=0
numrun=6 # dump root, and cuts 

while read line
do
  echo "do " ${i} ${line}

  rm -f runtmp_${i}.runs
  touch runtmp_${i}.runs
  echo ${line} >> runtmp_${i}.runs

  com="./do_run_heates_single.sh ${numrun} runtmp_${i}.runs"
  echo "com :" ${com}

  # do this job background 
  eval ${com} &


  # get the PID 
  pid[${i}]=$!

  # if PID > maxnproc, 
  i=`expr $i \+ 1`;
  if [ ${maxnproc} -eq ${i} ];then
    wait;
    i=0;
  fi
done < ${OUTPUT_FILE}

echo " finish run "