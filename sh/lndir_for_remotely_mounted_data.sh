# check the mount dir
mountdir="/jp"
if [ -z "$(ls -A ${mountdir})" ]; then
   echo "mounting ${mountdir}, need password of oper at jp..."
   sshfs oper@10.105.54.109:/heates/tes ${mountdir} -o ro
else
   echo "${mountdir} exists! go ahead..."
fi

mountdir2="/jp2/TMU_2018U"
if [ -z "$(ls -A ${mountdir2})" ]; then
   echo "mounting ${mountdir2}, need password of oper at jp2..."
   sshfs oper@10.105.54.109:/heates2/TMU_2018U ${mountdir2} -o ro
else
   echo "${mountdir2} exists! go ahead..."
fi


#echo "Did you finish the rsync of raw data to jp?"
#read -p "If yes, please hit the eneter key: "

#testrun="Timing"
#testrun="TMU_2018U"
# if testrun is empty "", this works for all data
#testrun=""
#LOCALDIR="/jp/${testrun}"
#LNDIR="/home/heates/data/${testrun}"

# from 2018/06/17 run279
LOCALDIR="/jp2/TMU_2018U"
LNDIR="/home/heates/data/TMU_2018U"


mkdir -p $LNDIR

if [ -d "${LNDIR}" ]; then
    # symbolic link
    echo "creating symbolic links"
    #lndir ${LOCALDIR} ${LNDIR}
    lndir -silent ${LOCALDIR} ${LNDIR}
else
    echo "${LNDIR} does not exit. something is wrong..."
    exit 1
fi
