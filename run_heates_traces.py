#!/usr/bin/env python

""" run_heates_traces.py is a basic analysis tools for HEATES project. 

History: 
2019-03-18 ; ver 1.0; made by HT

"""

__author__ =  'H.Tatsuno'
__version__ = '1.0'

import matplotlib
matplotlib.use("Agg")
import matplotlib.cm as cm
import matplotlib.colors as colors
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
import optparse

import khe_ana as KHE
import khe_util as util
import tesmap as tesmap
KHE = reload(KHE) # for ipython to reload when it changed
util = reload(util)
tesmap = reload(tesmap)


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
#BADCHS.extend([263,293,295])# group trig strange? for run381-merged

maxchans = 240
#maxchans = 10

if os.path.isdir(DATADIR)==False: 
    print "%s is missing"%DATADIR
    sys.exit()

RUNINFO="./csv/data_TMU_2018U.csv"
if os.path.exists(RUNINFO)==False: 
    print "%s is missing"%RUNINFO
    sys.exit(0)

#GRTINFO="./csv/grptrig_twocol.txt"
GRTINFO="./csv/grptrig_wo_neighbor.csv"
if os.path.exists(GRTINFO)==False: 
    print "%s is missing"%GRTINFO
    sys.exit(0)
    
COLUMN_INFO = "./csv/column_info.csv"
if os.path.exists(COLUMN_INFO)==False:
    print "%s is missing"%COLUMN_INFO
    sys.exit()


usage = u'(1) %prog 278 (basic analysis), (2) %prog 278'
version = __version__
parser = optparse.OptionParser(usage=usage, version=version)
# setting for options
parser.add_option('-f', '--force',    dest='forceNew',   action='store_true',  help='True to update filter (default=False)',      metavar='FORCE',   default=False)
parser.add_option('-s', '--summary',  dest='summaryNew', action='store_true',  help='True to update summary (default=False)',     metavar='SUMMARY', default=False)
parser.add_option('-c', '--calib',    dest='calibNew',   action='store_true',  help='True to update calibration (default=False)', metavar='CALIB',   default=False)
parser.add_option('-e', '--exttrig',  dest='externTrig', action='store_true',  help='True for calc externTrig (default=False)',   metavar='EXTTRIG', default=False)
parser.add_option('-g', '--grptrig',  dest='groupTrig',  action='store_true',  help='True for calc groupTrig (default=False)',    metavar='GRPTRIG', default=False)
parser.add_option('-d', '--delete',   dest='delete',     action="store_true",  help='True to delete hdf5 file (default=False)',   metavar='DELETE',  default=False)
parser.add_option('-m', '--merged',   dest='merged',     action="store_true",  help='use merged data (default=False)', metavar='MERGED', default=False)
parser.add_option('-R', '--dumproot', dest='dumproot',   action="store_true",  help='dump ROOT except for pulses (default=False)',metavar='DUMPROOT',default=False)
parser.add_option('-E', '--rootext',  dest='rootext',    action='store_true',  help='True ROOT with externTrig (default=False)',  metavar='ROOTEXT', default=False)
parser.add_option('-G', '--rootgrp',  dest='rootgrp',    action='store_true',  help='True ROOT with groupTrig (default=False)',   metavar='ROOTGRP', default=False)
parser.add_option('--beam',  dest='beam', action="store",type=str, help='set beam catecut (default=None, on or off)',default="None")
parser.add_option('--sprmc',  dest='sprmc', action="store",type=str, help='set sprmc catecut (default=None, on or off)',default="None")
parser.add_option('--jbrsc',  dest='jbrsc', action="store",type=str, help='set jbrsc catecut (default=None, on or off)',default="None")
parser.add_option('--pre',  dest='cut_pre', action="store",type=int, help='set cut for pre samples',default=0)
parser.add_option('--post',  dest='cut_post', action="store",type=int, help='set cut for post samples',default=0)
parser.add_option('--hdf5optname',  dest='hdf5optname', action="store", type=str, help='add optional name for hdf5 (default=None)', default=None)

options,args = parser.parse_args()

#### get options ####
FORCE          = options.forceNew
SUMMARY        = options.summaryNew
CALIB          = options.calibNew
EXTTRIG        = options.externTrig
GRPTRIG        = options.groupTrig
DELETE         = options.delete
MERGED         = options.merged
DUMPROOT       = options.dumproot
ROOTEXT        = options.rootext
ROOTGRP        = options.rootgrp
beam           = options.beam
sprmc          = options.sprmc
jbrsc          = options.jbrsc
cut_pre        = options.cut_pre
cut_post       = options.cut_post
hdf5optname    = options.hdf5optname

catecut = {}
catecut["prime"] = "on"
if not beam=="None":
    catecut["beam"] = beam
if not sprmc=="None":
    catecut["sprmc"] = sprmc
if not jbrsc=="None":
    catecut["jbrsc"] = sprmc
if beam=="None" and sprmc=="None":
    if jbrsc=="None":
        catecut=None

print ""
print "--- [OPTIONS] ----------------------"
print "  (standard)"
print "    FORCE          = ", FORCE
print "    SUMMARY        = ", SUMMARY
print "    CALIB          = ", CALIB
print "    EXTTRIG        = ", EXTTRIG
print "    GRPTRIG        = ", GRPTRIG
print "    DELETE         = ", DELETE
print "    MERGED         = ", MERGED
print "    catecut        = ", catecut
print "    cut_pre        = ", cut_pre
print "    cut_post       = ", cut_post
print "    hdf5optname    = ", hdf5optname
print ""
print "  (ROOT)             "
print "    DUMPROOT       = ", DUMPROOT
print "    ROOTEXT        = ", ROOTEXT
print "    ROOTGRP        = ", ROOTGRP
print "------------------------------------"


npar = len(args)
if (npar>=1):
    run_p = str(args[0])# 0001 or run0001
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
run_n="%04d"%irn

ana_target = ana_list[ind]
exttrig = exttrig_list[ind]
grptrig = grptrig_list[ind]

if ana_target=="":
    print "Error: ana is empty"
    sys.exit(0)
elif ana_list[ind]=="noise":
    print "Error: this is noise run, please select pulse run"
    sys.exit(0)

cal_run = None if str(cal_list[ind]) == "None" else int(cal_list[ind])
analist = [(int(run_p),int(run_n), ana_target, exttrig, grptrig, cal_run, BADCHS)]

for pulse_runnum, noise_runnum, target, extflag, grpflag, calibration_runnum, badchan in analist:
    print "..... run, noise, target, exttrig, grptrig, cal_run, = ", pulse_runnum, noise_runnum, target, extflag, grpflag, calibration_runnum
    print "BADCHAN = ", badchan
    if extflag == "off": # when external trigger or group trigger is off, those flags are forced to be false. 
    	orgEXTTRIG = EXTTRIG
        EXTTRIG = False
    if grpflag == "off":
        orgGRPTRIG = GRPTRIG
        GRPTRIG = False
#    if grouptrigmax: 
#        GRTINFO = "./csv/grptrig_singlethread_all.txt"
    # ---------------------------------------------------------
    k = KHE.KHE(pulse_runnum, noise_runnum, maxchans, calibration_runnum, badchan, DATADIR, DELETE, GRTINFO, COLUMN_INFO,
                hdf5optname=hdf5optname, catecut=catecut, target=target, cut_pre=cut_pre, cut_post=cut_post,
                merged=MERGED)
    k.anahide(forceNew=FORCE,summaryNew=SUMMARY,calibNew=CALIB,exttrigNew=EXTTRIG,grptrigNew=GRPTRIG)
    # ---------------------------------------------------------
    #if DUMPROOT:
    #    print "\n [dump ROOT except for pulses]"
    #    k.dump_ROOT_2018(EXTTRIG=ROOTEXT, GRTRIG=ROOTGRP)
    # ---------------------------------------------------------
    if extflag == "off": # when external trigger or group trigger is off, those flags are back to the input values 
    	EXTTRIG = orgEXTTRIG
    if grpflag == "off":
    	GRPTRIG = orgGRPTRIG
    # ---------------------------------------------------------


data = k.data


# common parameters
usechans=k.usechans
ds0 = data.first_good_dataset
util.init_row_timebase(ds0)
timebase       = util.ROW_TIMEBASE * util.NUM_ROWS
row_timebase   = util.ROW_TIMEBASE
number_of_rows = util.NUM_ROWS
row_offset     = util.GLOBAL_PT_OFFSET
nsamples    = ds0.nSamples
npresamples = ds0.nPresamples
timestamp_offset = ds0.timestamp_offset


# clock
import time
start = time.time()

np_each=2000# number of pulses to use for each ch
npls=0# number of total pulses
for ds in data:
    #npls += ds.nPulses
    npls += np_each

trec = np.zeros((npls,),
                dtype = [('ev', int),
                         ('traces', float,nsamples)# becareful for all columns
                         ])

rec = np.zeros((npls,),
               dtype = [('idx', int),
                        ('ind', int),
                        ('ch', int),
                        ('row', int),
                        ('col', int),
                        ('twr', int),
                        ('pix', int),
                        ('good', bool),
                        ('primary', bool),
                        ('beam', bool),
                        ('filt_value', float),
                        ('filt_value_dc', float),
                        ('filt_value_phc', float),
                        ('energy', float),
                        ('energy_dc', float),
                        ('energy_phc', float),
                        ('pretrig_mean', float),
                        ('pretrig_rms', float),
                        ('pulse_average', float),
                        ('pulse_rms', float),
                        ('peak_value', float),
                        ('peak_time', float),
                        ('postpeak_deriv', float),
                        ('rise_time', float),
                        ('timestamp', float),
                        ('rowcount', float),
                        ('sec_pr_mean', float)
               ])

trec['ev'] = np.arange(npls)

ev=0
for ds in data:
    print "channel %d start.... for %.3f"%(ds.channum, (time.time() - start))
    #if ds.column_number != 0: continue
    #npl = ds.nPulses
    npl = np_each
    row = ds.row_number
    twr, pix = tesmap.getTwrPixAdd(row)

    np_good=np.array(ds.good())
    np_prim=np.zeros(shape=np_good.shape,dtype=bool)
    np_grt=np.array(ds.p_grouptrig)
    prim_ind=np.where(np_grt==-1)[0]
    np_prim[prim_ind]=True

    np_filt_value_dc  = np.array(ds.p_filt_value_dc)
    np_filt_value_phc = np.array(ds.p_filt_value_phc)
    row_adjust = -1.*(np.array(ds.p_filt_phase)+np.array(ds.p_shift1))
    
    rec['ind'][ev:ev+npl]            = np.arange(npl)
    rec['ch'][ev:ev+npl]             = np.full(npl,ds.channum)
    rec['row'][ev:ev+npl]            = np.full(npl,row)
    rec['col'][ev:ev+npl]            = np.full(npl,ds.column_number)
    rec['twr'][ev:ev+npl]            = np.full(npl,twr)
    rec['pix'][ev:ev+npl]            = np.full(npl,pix)
    rec['good'][ev:ev+npl]           = np_good[0:npl]
    rec['primary'][ev:ev+npl]        = np_prim[0:npl]
    rec['beam'][ev:ev+npl]           = np.array(ds.beamflag)[0:npl]
    rec['filt_value'][ev:ev+npl]     = np.array(ds.p_filt_value)[0:npl]
    rec['filt_value_dc'][ev:ev+npl]  = np_filt_value_dc[0:npl]
    rec['filt_value_phc'][ev:ev+npl] = np_filt_value_phc[0:npl]
    rec['energy'][ev:ev+npl]         = np.array(ds.p_energy)[0:npl]
    attr="p_filt_value_dc"
    if ds.calibration.has_key(attr):
        np_energy_dc = ds.calibration[attr].ph2energy(np_filt_value_dc)
        rec['energy_dc'][ev:ev+npl]  = np_energy_dc[0:npl]
    attr="p_filt_value_phc"
    if ds.calibration.has_key(attr):
        np_energy_phc = ds.calibration[attr].ph2energy(np_filt_value_phc)
        rec['energy_phc'][ev:ev+npl] = np_energy_phc[0:npl]
    rec['pretrig_mean'][ev:ev+npl]   = np.array(ds.p_pretrig_mean)[0:npl]
    rec['pretrig_rms'][ev:ev+npl]    = np.array(ds.p_pretrig_rms)[0:npl]
    rec['pulse_average'][ev:ev+npl]  = np.array(ds.p_pulse_average)[0:npl]
    rec['pulse_rms'][ev:ev+npl]      = np.array(ds.p_pulse_rms)[0:npl]
    rec['peak_value'][ev:ev+npl]     = np.array(ds.p_peak_value)[0:npl]
    rec['peak_time'][ev:ev+npl]      = np.array(ds.p_peak_index)[0:npl]*timebase
    rec['postpeak_deriv'][ev:ev+npl] = np.array(ds.p_postpeak_deriv)[0:npl]
    rec['rise_time'][ev:ev+npl]      = np.array(ds.p_rise_time)[0:npl]
    rec['timestamp'][ev:ev+npl]      = row_timebase*(np.array(ds.p_rowcount)[0:npl]+(row_adjust[0:npl]*number_of_rows)-row_offset)
    rec['rowcount'][ev:ev+npl]       = row_timebase*(np.array(ds.p_rowcount)[0:npl]-row_offset)
    rec['sec_pr_mean'][ev:ev+npl]    = np.array(ds.sec_pr_mean)[0:npl]
    
    apm = np.array(ds.p_pretrig_mean)[0:npl]
    for i,apmi in enumerate(apm):
        atr = np.array(ds.read_trace(i))
        trec['traces'][ev] = atr-apmi
        ev += 1
    
print "sort start..."
sidx = np.argsort(rec, order='rowcount')
rec  = np.sort(rec, order='rowcount')
rec['idx'] = np.arange(npls)
print "==== End. %.3f"%((time.time() - start))

##### for traces plots ######
x = timebase*1e3*np.arange(nsamples)# ms
tr   = trec['traces']
ev     = sidx# sort index = event ID before sort
idx    = rec['idx']# pulse indecies
col    = rec['col']
ch     = rec['ch']
twr    = rec['twr']
pix    = rec['pix']
g      = rec['good']
prim   = rec['primary']
beam   = rec['beam']
flv    = rec['filt_value']
flvdc  = rec['filt_value_dc']
flvphc = rec['filt_value_phc']
ene    = rec['energy']
enedc  = rec['energy_dc']
enephc = rec['energy_phc']
premean = rec['pretrig_mean']
prerms  = rec['pretrig_rms']
plsave = rec['pulse_average']
plsrms = rec['pulse_rms']
pvalue = rec['peak_value']
ptime  = rec['peak_time']*1e3# ms
rtime  = rec['rise_time']*1e3# ms
sprm   = rec['sec_pr_mean']
#time   = rec['timestamp']*1e3# ms
time   = rec['rowcount']*1e3# ms
time = time - time[0]# common offset

prim_cut = prim==True

save = True
if save:
    import commands
    dir_save = "/Users/tatsuno/Desktop/hogehoge"
    commands.getoutput('mkdir -p ' + dir_save)
    savename = "%s/hogehoge_simultaneous_hits"%(dir_save)


#intvl=25.#ms
#intvl_offset=3000.#ms
#for j in xrange(100):
#    tcut = np.logical_and(time>intvl*j+intvl_offset,time<intvl*(j+1)+intvl_offset)
#    ind2 = idx[np.logical_and(tcut,prim_cut)]
#    plt.figure(figsize=(8.27,11.69))
#
#    ax = plt.subplot(2,1,1, ylim=(-2000,35000))
#    for ii in ind2:
#        plt.plot(x+time[ii],tr[ev[ii],:],color=tesmap.getColColor(ch[ii]))
#        ax.annotate("%d"%idx[ii],xy=(time[ii]+ptime[ii],flv[ii]),\
#                    xytext=(time[ii]+ptime[ii],25000),arrowprops=dict(arrowstyle="->"),size=8)
#    plt.xlabel("Time [ms]")
#    plt.ylabel("Pulse FB - pretrig_mean")
#
#    ax = tesmap.getPlotRectPlane(2,1,2,usechans)
#    for ii in ind2:
#        tesmap.writeTextTESPos(ax,ch[ii],idx[ii],fontsize=8)
#    plt.xlabel("Position x [$\mu$m]")
#    plt.ylabel("Position y [$\mu$m]")
#    plt.tight_layout()
#    plt.show()
#    if save: plt.savefig("%s_%d"%(savename,plt.get_fignums()[-1]),ext="png")


intvl=20.#ms
intvl_offset=3000.#ms
for j in xrange(100):
    tcut = np.logical_and(time>intvl*j+intvl_offset,time<intvl*(j+1)+intvl_offset)
    ind2 = idx[np.logical_and(tcut,prim_cut)]
    plt.figure(figsize=(8.27,11.69))

    ax = plt.subplot(3,1,1, ylim=(-2000,35000))
    for i2 in ind2:
        plt.plot(x+time[i2],tr[ev[i2],:],color=tesmap.getColColor(ch[i2]))
        ax.annotate("%d"%idx[i2],xy=(time[i2]+ptime[i2],flv[i2]),\
                    xytext=(time[i2]+ptime[i2],25000),arrowprops=dict(arrowstyle="->"),size=8)
    plt.xlabel("Time [ms]")
    plt.ylabel("Pulse FB - pretrig_mean")

    ax = plt.subplot(3,1,2, ylim=(-100,400))
    for i2 in ind2:
        plt.plot(x+time[i2],tr[ev[i2],:],color=tesmap.getColColor(ch[i2]))
    ind3 = idx[np.logical_and(tcut,~prim_cut)]
    for i3 in ind3:
        plt.plot(x+time[i3],tr[ev[i3],:],color=tesmap.getColColor(ch[i3]))
    plt.xlabel("Time [ms]")
    plt.ylabel("Pulse FB - pretrig_mean")
    
    ax = tesmap.getPlotRectPlane(3,1,3,usechans)
    for ii in ind2:
        tesmap.writeTextTESPos(ax,ch[ii],idx[ii],fontsize=8)
    plt.xlabel("Position x [$\mu$m]")
    plt.ylabel("Position y [$\mu$m]")
    plt.tight_layout()
    plt.show()
    if save: plt.savefig("%s_%d"%(savename,plt.get_fignums()[-1]),ext="png")

    

# end
