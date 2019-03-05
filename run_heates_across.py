#!/usr/bin/env python

""" run_heates_across.py is a basic analysis tools for HEATES project. 

History: 
2018-07-08 ; ver 1.0; created from ana_galen_across.py

"""

__author__ =  'HEATES'
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
from os import path
import sys
import h5py
import datetime
import pandas as pd
from collections import OrderedDict
import optparse

import khe_ana as KHE
import khe_util as khe
KHE= reload(KHE) # for ipython to reload when it changed
khe= reload(khe) # for ipython to reload when it changed

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

if os.path.isdir(DATADIR)==False: 
    print "%s is missing"%DATADIR
    sys.exit()

usage = u'(1) %prog 278 (basic analysis), (2) %prog 278 --energy_cut --energyMin=6900. --energyMax=6975 --cut_category=sec_pr_mean --cutMax=37 --doiteach --pulse_ana --chans=129 (when plot pulses)'
version = __version__
parser = optparse.OptionParser(usage=usage, version=version)
# setting for options
parser.add_option('-t',  dest='target', action="store", type=str, help='set target He3 or He4 (default=He3)', default="He3")
parser.add_option('-c',  dest='cut', action="store_true", help='use cut flag (default=False)', default=False)

options,args = parser.parse_args()

#### get options ####
target    = options.target
cut    = options.cut

print ""
print "--- [OPTIONS] ----------------------"
print "    target    = ", target
print "    cut       = ", cut
print "------------------------------------"

npar = len(args)
if (npar>=1):
    outputname = str(args[0])#
else:
    print "Error: specify outputname ", args
    sys.exit(0)

if target == "He3":
    He4=False
    infile = open("He3.runs","r")
else:
    He4=True
    infile = open("He4.runs","r")

print "..... target is set ", target

# get He3 and He4 runlist 

p_runs = []
for runs in infile:   
    for i, run in enumerate(runs.split()):
	print "..... ", i, " run = ", run 
        p_runs.append(int(run))

RUNINFO="./csv/data_TMU_2018U.csv"
if os.path.exists(RUNINFO)==False: 
    print "%s is missing"%RUNINFO
    sys.exit(0)

df = pd.read_csv(RUNINFO)
run_list = df.iloc[:,0].tolist()
noise_list = df.iloc[:,1].tolist()

inds = [run_list.index(run) for run in p_runs]    
n_runs = [int(noise_list[ind]) for ind in inds]
nr=len(p_runs)
    
# create figures
khe.hist2d_across_runs_merged(p_runs,n_runs,closeFigs=True)


# -------- for ROOT file --------------
try:
    import ROOT
except ImportError:
    raise ValueError('ERROR: cannot import pyROOT')
ROOT.gROOT.SetBatch(1)

def get_hist(x,y,c,hname):
    timex = x[0,:]
    nbinx = len(timex)
    binwx = timex[1] - timex[0]
    eney  = y[:,0]
    nbiny = len(eney)
    binwy = eney[1] - eney[0]
    counts = c.T
    hene = ROOT.TH2F(hname,hname,nbinx,timex[0],timex[-1]+binwx,nbiny,eney[0],eney[-1]+binwy)
    for ix in xrange(len(timex)-1):
        for iy in xrange(len(eney)-1):
            hene.SetBinContent(ix+1,iy+1,counts[ix,iy])
    return hene


# create ROOT file
#Fe, 3He, 4He
outd = khe.get_across_XYC_merged(p_runs,n_runs,closeFigs=True)
hnames=["hene_fe","hene_he3","hene_he4"]
output_dir="./output/"

basename = "_rows_merged_nr%d_%d_%d_%s.root"%(nr,p_runs[0],p_runs[-1], outputname)
if He4:
    fname="He4" + basename
    if cut:
        fname_cut="cut_He4" + basename
else:
    fname="He3" + basename
    if cut:
        fname_cut="cut_He3" + basename

fout=path.join(output_dir,fname)
fopt = "recreate"
f = ROOT.TFile(fout,fopt)
for i,hn in enumerate(hnames):
    x,y,c = outd.values()[i]
    h = get_hist(x,y,c,hn)
    h.Write()
f.Close()
print "%s is created."%fout

if not cut: 
    sys.exit()

# The following will be performed only when cut flag is True. 

khe.hist2d_across_runs_merged_with_cut(p_runs,n_runs,closeFigs=True)
outd_cut = khe.get_across_XYC_merged_with_cut(p_runs,n_runs,closeFigs=True)
fout2=path.join(output_dir,fname_cut)
fopt = "recreate"
f = ROOT.TFile(fout2,fopt)
for i,hn in enumerate(hnames):
    x,y,c = outd_cut.values()[i]
    h = get_hist(x,y,c,hn)
    h.Write()
f.Close()
print "%s is created."%fout2















































