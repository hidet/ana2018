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
    if savename!=None: plt.savefig("./fig/"+savename+".pdf")


    
datadir="/Users/tatsuno/work/heates/data/TMU_2018U/dumproot/run0397/"

rlcuts=["",
       "pre000_post100_",
       "pre000_post200_",
       "pre000_post300_",
       "pre000_post400_",
       "pre000_post450_",
       "pre000_post500_",
       "pre000_post550_",
       "pre000_post600_",
       "pre050_post550_",
       "pre100_post550_",
       "pre150_post550_",
       "pre200_post550_"]


#fnames=["run0397_noi0393_mass_2018_%ssprmcon_jbrscon_prime"%(rlc) for rlc in rlcuts]
#fnames=["run0397_noi0393_mass_2018_%sspillon_sprmcon_jbrscon_prime"%(rlc) for rlc in rlcuts]
fnames=["run0397_noi0393_mass_2018_%sspilloff_sprmcon_jbrscon_prime"%(rlc) for rlc in rlcuts]

#template=""
#template="spillon"
template="spilloff"

fin = ROOT.TFile.Open(datadir+"run0397_%s_hist.root"%(template),"read")
hene_offs = [fin.Get("hene_off_%s"%(rlc)) for rlc in rlcuts]
hene_ons  = [fin.Get("hene_on_%s"%(rlc)) for rlc in rlcuts]
for fname,rlc,hoff,hon in zip(fnames,rlcuts,hene_offs,hene_ons):
    linenames=["CrKAlpha","CoKAlpha","CuKAlpha"]
    for linename in linenames:
        plot_linefit(linename,hoff,hon,title="run0397 %s %s %s"%(rlc,template,linename),savename=fname+"_"+linename)
