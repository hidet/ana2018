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


def get_hist_bins(h,elo,ehi):
    bin_s=h.FindBin(elo)
    bin_e=h.FindBin(ehi)
    bins=np.arange(bin_s,bin_e+1)
    nbin = bin_e - bin_s + 1
    bw=h.GetBinWidth(bin_s)
    minx=h.GetBinLowEdge(bin_s)
    maxx=h.GetBinLowEdge(bin_e)+bw
    bin_edges=np.arange(minx,maxx+bw,bw)
    hist = [h.GetBinContent(b) for b in bins]
    hist = np.array(hist)
    return hist,bin_edges

def linefit(hist,bin_edges,linename):
    fitter = mass.getfitter(linename)
    fitter.fit(hist,bin_edges,plot=False)
    params = fitter.last_fit_params[:]
    params[fitter.param_meaning["tail_frac"]]=0.25
    params[fitter.param_meaning["dP_dE"]]=1
    params[fitter.param_meaning["resolution"]]=6
    fitter.fit(hist,bin_edges,params,hold=[2],vary_tail=True, plot=False)
    return fitter

def plot_linefit(linename,hoff,hon,title="",savename=None):
    elo,ehi = np.array([-50,50]) + mass.STANDARD_FEATURES[linename]
    hist, bin_edges = get_hist_bins(hoff,elo,ehi)
    fitter = linefit(hist,bin_edges,linename)
    plt.figure()
    ax1=plt.subplot(211)
    fitter.plot(axis=ax1,label=False,color="r",ph_units="eV")
    plt.xlabel("Energy (eV)")
    plt.ylabel("Counts per {:.2f} eV bin".format(bin_edges[1]-bin_edges[0]))
    res = fitter.last_fit_params_dict["resolution"][0]
    tail_frac = fitter.last_fit_params_dict["tail_frac"][0]
    leg=plt.legend(["beam off","fit: {:.2f} eV, {:.0f}% tail".format(res, tail_frac*100)],loc="best",frameon=False)
    leg.get_frame().set_linewidth(0.0)
    plt.title(title)

    hist, bin_edges = get_hist_bins(hon,elo,ehi)
    fitter = linefit(hist,bin_edges,linename)
    ax2=plt.subplot(212)
    fitter.plot(axis=ax2,label=False,color="b",ph_units="eV")
    plt.xlabel("Energy (eV)")
    plt.ylabel("Counts per {:.2f} eV bin".format(bin_edges[1]-bin_edges[0]))
    res = fitter.last_fit_params_dict["resolution"][0]
    tail_frac = fitter.last_fit_params_dict["tail_frac"][0]
    leg=plt.legend(["beam on","fit: {:.2f} eV, {:.0f}% tail".format(res, tail_frac*100)],loc="best",frameon=False)
    leg.get_frame().set_linewidth(0.0)
    plt.tight_layout()
    if savename!=None: plt.savefig(savename+".pdf")



runnum=421
datadir="/Users/tatsuno/work/heates/data/TMU_2018U/dumproot/run%04d/"%runnum

rlcuts=["pre150_post550_"]

#fnames=["run0397_noi0393_mass_2018_%ssprmcon_jbrscon_prime"%(rlc) for rlc in rlcuts]
fnames=["run%04d_noi0419_mass_2018_%sspillon_sprmcon_jbrscon_prime"%(runnum,rlc) for rlc in rlcuts]
#fnames=["run0397_noi0393_mass_2018_%sspilloff_sprmcon_jbrscon_prime"%(rlc) for rlc in rlcuts]

#template=""
template="spillon"
#template="spilloff"

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


linenames=["CrKAlpha","CoKAlpha","CuKAlpha"]
for t,fname,rlc,hoff,hon in zip(trees,fnames,rlcuts,hene_offs,hene_ons):
    # for speed up: create tmp file to save the tree with basic cut
    ftmpname=datadir+fname+"_cut2.root"
    if (os.path.isfile(ftmpname) is False):
        ftmp = ROOT.TFile.Open(ftmpname,"recreate")
        newt = t.CopyTree(cut_base)
        print "%s has been created with %s"%(ftmpname,treename)
    else:
        ftmp = ROOT.TFile.Open(ftmpname,"read")
        newt = ftmp.Get(treename)
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
    for linename in linenames:
        plot_linefit(linename,hoff,hon,title="run%04d %d pixels %s %s %s"%(runnum,nch,rlc,template,linename),savename=fname+"_"+linename)


#fout = ROOT.TFile.Open(datadir+"run0397_%s_hist.root"%(template),"recreate")
#for hoff,hon in zip(hene_offs,hene_ons):
#    hoff.Write()
#    hon.Write()
#fout.Close()
