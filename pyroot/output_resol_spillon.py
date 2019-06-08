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

    
datadir="/Users/tatsuno/work/heates/data/TMU_2018U/dumproot/run0397/"
template="spillon"
rlcuts=[
    "spillon",
    "pre000_post100_spillon",
    "pre000_post200_spillon",
    "pre000_post300_spillon",
    "pre000_post350_spillon",
    "pre000_post400_spillon",
    "pre000_post450_spillon",
    "pre000_post500_spillon",
    "pre000_post550_spillon",
    "pre000_post600_spillon",
    "pre050_post350_spillon",
    "pre050_post400_spillon",
    "pre050_post450_spillon",
    "pre050_post500_spillon",
    "pre050_post550_spillon",
    "pre100_post350_spillon",
    "pre100_post400_spillon",
    "pre100_post450_spillon",
    "pre100_post500_spillon",
    "pre100_post550_spillon",
    "pre150_post350_spillon",
    "pre150_post400_spillon",
    "pre150_post450_spillon",
    "pre150_post500_spillon",
    "pre150_post550_spillon",
    "pre200_post350_spillon",
    "pre200_post400_spillon",
    "pre200_post450_spillon",
    "pre200_post500_spillon",
    "pre200_post550_spillon"]

rlpre=[0,0,0,0,0,0,0,0,0,0,
       50,50,50,50,50,
       100,100,100,100,100,
       150,150,150,150,150,
       200,200,200,200,200]

rlpost=[0,100,200,300,350,400,450,500,550,600,
        350,400,450,500,550,
        350,400,450,500,550,
        350,400,450,500,550,
        350,400,450,500,550]

fnames=["run0397_noi0393_mass_2018_%s_sprmcon_jbrscon_prime"%(rlc) for rlc in rlcuts]

fin = ROOT.TFile.Open(datadir+"run0397_%s_hist.root"%(template),"read")
hene_offs = [fin.Get("hene_off_%s"%(rlc)) for rlc in rlcuts]
hene_ons  = [fin.Get("hene_on_%s"%(rlc)) for rlc in rlcuts]
linenames=["CrKAlpha","CoKAlpha","CuKAlpha"]

import csv

for linename in linenames:
    elo,ehi = np.array([-50,50]) + mass.STANDARD_FEATURES[linename]
    
    f_off = open('%s_beamoff_spillon.csv'%(linename), 'w')
    writer_off = csv.writer(f_off, lineterminator='\n')
    print ""
    for fname,rlc,pre,post,h in zip(fnames,rlcuts,rlpre,rlpost,hene_offs):
        # hoff
        hist, bin_edges = util.get_hist_bins(h,elo,ehi)
        fitter = util.linefit(hist,bin_edges,linename)
        res = fitter.last_fit_params_dict["resolution"][0]
        res_err = fitter.last_fit_params_dict["resolution"][1]
        outlist=[pre,post,res,res_err]
        writer_off.writerow(outlist)
        print "%s beamoff %03d %03d %.2f +- %.2f eV"%(linename,pre,post,res,res_err)
    f_off.close()


    f_on = open('%s_beamon_spillon.csv'%(linename), 'w')
    writer_on = csv.writer(f_on, lineterminator='\n')
    print ""
    for fname,rlc,pre,post,h in zip(fnames,rlcuts,rlpre,rlpost,hene_ons):
        # hon
        hist, bin_edges = util.get_hist_bins(h,elo,ehi)
        fitter = util.linefit(hist,bin_edges,linename)
        res = fitter.last_fit_params_dict["resolution"][0]
        res_err = fitter.last_fit_params_dict["resolution"][1]
        outlist=[pre,post,res,res_err]
        writer_on.writerow(outlist)
        print "%s beamon %03d %03d %.2f +- %.2f eV"%(linename,pre,post,res,res_err)
    f_on.close()
    
