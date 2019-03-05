#!/usr/bin/env python

""" check_heates_grptrig.py is to check group trigger. 

History: 
2018-08-03 ; ver 1.0; created from ana_check_grptrig.py

"""

__author__ =  'Shinya Yamada (syamada(at)tmu.ac.jp'
__version__ = '1.0'

import matplotlib
matplotlib.use("Agg")
import mass
import monkeypatch # patches mass to avoid crashes
import numpy as np
import pylab as plt
params = {'xtick.labelsize': 10, # x ticks
          'ytick.labelsize': 10, # y ticks
          'legend.fontsize': 8
                    }
plt.rcParams['font.family'] = 'serif' # 
plt.rcParams.update(params)
import os
import sys
import h5py
import datetime
import pandas as pd
from collections import OrderedDict
import optparse
import commands
import khe_ana as KHE
KHE= reload(KHE) # for ipython to reload when it changed

print "[START] " + __file__

ANADIR=os.environ.get("HEATESANADIR","")
DATADIR=os.environ.get("HEATESDATADIR","")
print "ANADIR  = ", ANADIR
print "DATADIR = ", DATADIR

if ANADIR == "" or DATADIR == "":
    print "[ERROR] Set HEATESANADIR and HEATESDATADIR"
    print "e.g., for bash users"
    print '''
    export HEATESHOME="$HOME/work/ana/HEATES/JPARC201807_E62"
    export HEATESANADIR="$HEATESHOME/ana_2018_jparc"
    export HEATESDATADIR="$HEATESHOME/data/TMU_2018U" '''.strip()
    sys.exit()

BADCHS = [3,9,39,77,83,85,111,337,367,375,423]# initially disconnected
BADCHS.extend([117,203,233])# bad channels
BADCHS.extend([5,177])# strange channels

maxchans = 240
# maxchans = 60
kdict = OrderedDict()

# prim_chan_list = [115]
getMax = 100
default_pulse_chan = 129
# cut_category = "energy"

cut_category_list = ["sec_pr_mean"]

if os.path.isdir(DATADIR)==False: 
    print "%s is missing"%DATADIR
    sys.exit()

RUNINFO = ANADIR + "/csv/data_TMU_2018U.csv"
if os.path.exists(RUNINFO)==False: 
    print "%s is missing"%RUNINFO
    sys.exit(0)

GRTINFO= ANADIR + "/csv/grptrig_twocol.txt"
#GRTINFO="../csv/grptrig_twocol.txt"
if os.path.exists(RUNINFO)==False: 
    print "%s is missing"%RUNINFO
    sys.exit(0)
setch=np.loadtxt(GRTINFO, delimiter=',')
    
COLUMN_INFO = ANADIR + "/csv/column_info.csv"

if os.path.exists(COLUMN_INFO)==False:
    print "%s is missing"%COLUMN_INFO
    sys.exit()

usage = u'%prog 278 (basic analysis)'
version = __version__
parser = optparse.OptionParser(usage=usage, version=version)
# setting for options
parser.add_option('--dumproot',  dest='dumproot', action="store_true", help='dump ROOT except for pulses (default=False, not implemented)', default=False)

options,args = parser.parse_args()

#### get options ####
dumproot       = options.dumproot

print "--- [OPTIONS] ----------------------"
print "  (ROOT)                  "
print "    dumproot       = ", options.dumproot
print "------------------------------------"

npar = len(args)
if (npar>=1):
    run_p = str(args[0])# 0001 or run0001
    run_p = run_p.replace("run","")
else:
    print "Error: specify run number of E62 ", args
    sys.exit(0)

df = pd.read_csv(RUNINFO)
run_list = df.iloc[:,0].tolist()
ind = run_list.index(int(run_p))
noise_list = df.iloc[:,1].tolist()
irn = int(noise_list[ind])
ana_list = df.iloc[:,2].tolist()
exttrig_list = df.iloc[:,3].tolist()
grptrig_list = df.iloc[:,4].tolist()
cal_list = df.iloc[:,9].tolist()

ana_target = ana_list[ind]
exttrig = exttrig_list[ind]
grptrig = grptrig_list[ind]

if ana_target=="":
    print "Error: ana is empty"
    sys.exit(0)
elif ana_list[ind]=="noise":
    print "Error: this is noise run, please select pulse run"
    sys.exit(0)

run_p="%04d"%int(run_p)
run_n="%04d"%irn

cal_run = None if str(cal_list[ind]) == "None" else int(cal_list[ind])
analist = [(int(run_p),int(run_n), ana_target, exttrig, grptrig, "calibration source on", cal_run, BADCHS)]


print '############################################\n'
print ' Checking Group Trigger\n'

h5name = DATADIR + "/run%s/run%s_mass.hdf5"%(run_p,run_p)

print "..... h5name = ", h5name

if os.path.isfile(h5name)==False: 
    print "%s is missing. please run basic analysis first."%h5name
    sys.exit()

h5f = h5py.File(h5name,'r')

commands.getoutput('mkdir -p log')
fname="./log/check_grptrig_"+str(run_p)+".log"
fp = open(fname, "w")

for sch in setch:
    cch = sch[0]
    ach = sch[1:]
    counter = {}

    try:

        mask = h5f['chan%d'%cch]['grouptrig'][:]==-1
        counter[cch] = len(mask[mask==True])
        print 'center channel: %d - grptrig count %d' %(cch, counter[cch])
        fp.write("%d, %d, "%(cch,counter[cch]))
    except KeyError:
        print 'Center channel not'
        continue

    lch = []
    for i in ach:
        try:
            mask = h5f['chan%d'%i]['grouptrig'][:]== cch
            counter[i] = len(mask[mask==True])
            fp.write(str(counter[i]) + ", ")
            lch.append(i)
                
        except KeyError:
            print 'There is not %d channel' %i
            continue

            
    for j in lch:
        if counter[cch]==counter[j]:
            continue
            #print 'Around channel: %d - grptrig count %d OK' %(j, counter[j])
        elif counter[j]==0:
            mask = h5f['chan%d'%j]['energy'][:]!=0
            count = len(mask[mask==True])
            #print 'Around channel: %d - Count %d' %(j, count)
        else:
            print 'Around channel: %d - grptrig count %d No' %(j, counter[j])
            fp.write(" ERROR")            

    fp.write("\n")                        

print '\n############################################\n'
fp.close()
print "saved as "+fname+"  make sure no ERROR\n"




