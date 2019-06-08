import numpy as np
import sys
import os
import h5py
import csv

import pyroot_util as util
util= reload(util)

'''
create all fitted peak positions as a csv file from hdf5 files
'''

#runs=util.runs_he3
runs=util.runs_he4

options,args = util.parser.parse_args()
cut_pre  = options.cut_pre
cut_post = options.cut_post
spill = options.spill

add=""
if cut_pre!=0 or cut_post!=0: add+="_pre%03d_post%03d"%(cut_pre,cut_post)
if spill==-1:
    print "Error: please specify --spill=0 or --spill=1"
    sys.exit(0)
elif spill==0: add+="_spilloff_sprmcon_jbrscon_prime"
elif spill==1: add+="_spillon_sprmcon_jbrscon_prime"

fnameins=["%s/run%04d/run%04d_noi%04d_mass_2018%s"%(util.datadir,run,run,noi,add) for run,noi in zip(runs,util.get_noise_list(runs))]
for fnamein in fnameins:
    if os.path.isfile(fnamein+".hdf5")==False:
        print "Error: file is missing %s"%(fnamein+".hdf5")
        print "forgetting options? --pre=xxx --post=yyy"
        sys.exit(0)

fnameout="%s/pp_run%04d_%04d%s"%(util.outdir,runs[0],runs[-1],add)
print "runs: ", runs
print "input files[0]: ", fnameins[0]+".hdf5"
print "output file:    ", fnameout+".csv"
util.check_continue()
#print "debug exit"
#sys.exit(0)# please remove this exit
# --------------------------------------------------------

fout = open(fnameout+".csv", 'w')
print "%s is opened as recreate"%(fnameout+".csv")
writer = csv.writer(fout, lineterminator='\n')
csvlist=[]
csvlist.append(["run","ch","CrKAlpha","CrKBeta","CoKAlpha","CoKBeta","CuKAlpha"])
cal_attrs=['p_filt_value_dc','p_filt_value_phc']
for j, (run,fnamein) in enumerate(zip(runs,fnameins)):
    fin = h5py.File(fnamein+".hdf5","r")
    chs = fin.keys()
    for ch in chs:
        cal_g = fin[ch]['calibration']
        fvs = cal_g.keys()
        if len(fvs)==2 and fvs[1]==cal_attrs[1]:
            cal = cal_g[fvs[1]]
            print run, ch, fvs[1]
        elif len(fvs)==1 and fvs[0]==cal_attrs[0]:
            cal = cal_g[fvs[0]]
            print run, ch, fvs[0]
        else:
            print "no calibration in %s"%ch
            continue
        phs=cal['ph'][:]
        if len(phs)==5:
            csvlist.append(["%03d"%run,"%03d"%(int(ch[4:])),"%.1f"%phs[0],"%.1f"%phs[1],"%.1f"%phs[2],"%.1f"%phs[3],"%.1f"%phs[4]])
        else:
            print "number of calibration peaks is strange in %s"%ch
            continue

    fin.close()

writer.writerows(csvlist)
fout.close()
    
print "%s has been created."%(fnameout+".csv")
# --------------------------------------------------------
