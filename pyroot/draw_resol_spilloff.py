from ROOT import gROOT
gROOT.SetBatch()

import ROOT
import mass

import h5py as h5
import pylab as plt
import numpy as np
import pandas as pd
import sys
import os
import optparse

import pyroot_util as util
util=reload(util)

def plot_linefit(linename,hoff,hon,title="",savename=None):
    elo,ehi = np.array([-50,50]) + mass.STANDARD_FEATURES[linename]
    hist, bin_edges = util.get_hist_bins(hoff,elo,ehi)
    fitter = util.linefit(hist,bin_edges,linename)
    res = fitter.last_fit_params_dict["resolution"][0]
    tail_frac = fitter.last_fit_params_dict["tail_frac"][0]
    plt.figure()
    ax1=plt.subplot(211)
    fitter.plot(axis=ax1,label=False,color="r",ph_units="eV")
    plt.xlabel("Energy (eV)")
    plt.ylabel("Counts per {:.2f} eV bin".format(bin_edges[1]-bin_edges[0]))
    leg=plt.legend(["beam off","fit: {:.2f} eV, {:.0f}% tail".format(res, tail_frac*100)],loc="best",frameon=False)
    leg.get_frame().set_linewidth(0.0)
    plt.title(title)

    hist, bin_edges = util.get_hist_bins(hon,elo,ehi)
    fitter = util.linefit(hist,bin_edges,linename)
    res = fitter.last_fit_params_dict["resolution"][0]
    tail_frac = fitter.last_fit_params_dict["tail_frac"][0]
    ax2=plt.subplot(212)
    fitter.plot(axis=ax2,label=False,color="b",ph_units="eV")
    plt.xlabel("Energy (eV)")
    plt.ylabel("Counts per {:.2f} eV bin".format(bin_edges[1]-bin_edges[0]))
    leg=plt.legend(["beam on","fit: {:.2f} eV, {:.0f}% tail".format(res, tail_frac*100)],loc="best",frameon=False)
    leg.get_frame().set_linewidth(0.0)
    plt.tight_layout()
    if savename!=None: plt.savefig("./fig/"+savename+".pdf")


datadir="/Users/tatsuno/work/heates/data/TMU_2018U/dumproot/run0397/"
template="spilloff"
rlcuts=[
    "spilloff",
    "pre000_post100_spilloff",
    "pre000_post200_spilloff",
    "pre000_post300_spilloff",
    "pre000_post350_spilloff",
    "pre000_post400_spilloff",
    "pre000_post450_spilloff",
    "pre000_post500_spilloff",
    "pre000_post550_spilloff",
    "pre000_post600_spilloff",
    "pre050_post300_spilloff",
    "pre050_post350_spilloff",
    "pre050_post400_spilloff",
    "pre050_post450_spilloff",
    "pre050_post500_spilloff",
    "pre050_post550_spilloff",
    "pre100_post300_spilloff",
    "pre100_post350_spilloff",
    "pre100_post400_spilloff",
    "pre100_post450_spilloff",
    "pre100_post500_spilloff",
    "pre100_post550_spilloff",
    "pre150_post300_spilloff",
    "pre150_post350_spilloff",
    "pre150_post400_spilloff",
    "pre150_post450_spilloff",
    "pre150_post500_spilloff",
    "pre150_post550_spilloff",
    "pre200_post300_spilloff",
    "pre200_post350_spilloff",
    "pre200_post400_spilloff",
    "pre200_post500_spilloff",
    "pre200_post550_spilloff"]

rlpre=[0,0,0,0,0,0,0,0,0,0,
       50,50,50,50,50,50,
       100,100,100,100,100,100,
       150,150,150,150,150,150,
       200,200,200,200,200]

rlpost=[0,100,200,300,350,400,450,500,550,600,
        300,350,400,450,500,550,
        300,350,400,450,500,550,
        300,350,400,450,500,550,
        300,350,400,500,550]

fnames=["run0397_noi0393_mass_2018_%s_sprmcon_jbrscon_prime"%(rlc) for rlc in rlcuts]

fnameins=[datadir+fn for fn in fnames]
treename="chanall"
files=[ROOT.TFile.Open(fname+".root") for fname in fnameins]
trees=[f.Get(treename) for f in files]

nbin=10000
minx=0.
maxx=10000.
hene_offs = [ROOT.TH1F("hene_off_%s"%(rlc),"hene_off_%s"%(rlc),nbin,minx,maxx) for rlc in rlcuts]
hene_ons  = [ROOT.TH1F("hene_on_%s"%(rlc),"hene_on_%s"%(rlc),nbin,minx,maxx) for rlc in rlcuts]

sprm_m =  -10
sprm_p =   10
jbrs_m = -400
jbrs_p =  350
good     = "good==1"
prime    = "primary==1"
spillon  = "beam==1"
spilloff = "beam==0"
sprmc    = "sec_pr_mean>%d && sec_pr_mean<%d"%(sprm_m,sprm_p)
jbrsc    = "jbr_region_sum>%d && jbr_region_sum<%d"%(jbrs_m,jbrs_p)
cut_base = "%s && %s && %s && %s"%(good,prime,sprmc,jbrsc)

linename="CoKAlpha"

for t,fname,rlc,hoff,hon in zip(trees,fnames,rlcuts,hene_offs,hene_ons):
    # for speed up: create tmp file to save the tree with basic cut
    ftmpname=datadir+fname+"_cut1.root"
    #if (os.path.isfile(ftmpname) is False):
    ftmp = ROOT.TFile.Open(ftmpname,"recreate")
    newt = t.CopyTree(cut_base)
    print "%s has been created with %s"%(ftmpname,treename)
    #else:
    #    ftmp = ROOT.TFile.Open(ftmpname,"read")
    #    newt = ftmp.Get(treename)
    chtmp=0
    nch=0
    for e in newt:
        ch=e.ch
        if chtmp<ch:# ch is changed
            nch+=1
            print "ch%d"%(ch)
            chtmp=ch
        if e.beam==1:
            hon.Fill(e.energy)
        elif e.beam==0:
            hoff.Fill(e.energy)
    plot_linefit(linename,hoff,hon,title="run0397 %d pixels %s %s"%(nch,rlc,template),savename=fname)


fout = ROOT.TFile.Open(datadir+"run0397_%s_hist.root"%(template),"recreate")
for hoff,hon in zip(hene_offs,hene_ons):
    hoff.Write()
    hon.Write()
fout.Close()
