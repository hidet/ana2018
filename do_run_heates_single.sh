#!/bin/sh 
# version 1.1, inc. multithread, 2018/07/11, S.Yamada
# version 1.2, inc. fix bug, 2018/07/12, R.Hayakawa

echo ".... start do_run_reates_single.sh"

# do run 
run(){
echo "run()"
inputf=$1
for runnum in `cat ${inputf}`
do
date
echo $runnum 
mkdir -p runlog
python run_heates_single.py $runnum ${options} > runlog/E62_${runnum}.log 2>&1 
done
}


###### set option #############################
if [ _$1 = _ ];
then
echo "[ERROR] need to specify option number"
echo "./do_run_heates_single.sh optioni xxx.runs"
echo "0) --dumproot"
echo "1) --doiteach"
echo "2) --doiteach --dumproot"
echo "3) --doiteach --dumprootpulse"
echo "4) --dumprootpulse"
echo "5) --cut_category=sec_pr_mean --cutMax=37"
echo "6) --cut_category=sec_pr_mean --cutMax=37 --dumproot"
exit 1 

else

optioni=$1

if [   "${optioni}" = "0" ];then
options="--dumproot"
elif [ "${optioni}" = "1" ];then
options="--doiteach"
elif [ "${optioni}" = "2" ];then
options="--doiteach --dumproot"
elif [ "${optioni}" = "3" ];then
options="--doiteach --dumprootpulse"
elif [ "${optioni}" = "4" ];then
options="--dumprootpulse"
elif [ "${optioni}" = "5" ];then
options="--cut_category=sec_pr_mean --cutMax=37"
elif [ "${optioni}" = "6" ];then
options="--cut_category=sec_pr_mean --cutMax=37 --dumproot"
else
options=""
fi

fi

echo ".... options is set " ${options}

################################################


###### start run ###############################

if [ _$2 = _ ];
then

echo "it automatically process for He3.runs and He4.runs"
run He3.runs & 
run He4.runs 

else

echo "it processes for input files"
inputfile=$2

echo "inputfile = " $inputfile

run $inputfile

fi


echo ".... end do_run_reates_single.sh"






