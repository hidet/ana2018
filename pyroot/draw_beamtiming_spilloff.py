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

datadir="/Volumes/HEATES_HD/jparc2018_data/root/run0397/"
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

nbin=150
minx=0.
maxx=150.
nbiny=450
miny=6150.
maxy=6600.
hdts  = [ROOT.TH1F("hdt_%s"%(rlc),"hdt_%s"%(rlc),nbin,minx,maxx) for rlc in rlcuts]
hdt2s = [ROOT.TH2F("hdt2_%s"%(rlc),"hdt2_%s"%(rlc),nbin,minx,maxx,nbiny,miny,maxy) for rlc in rlcuts]

sprm_m =  -10
sprm_p =   10
jbrs_m = -400
jbrs_p =  350

for t,fname,rlc,hdt,hdt2 in zip(trees,fnames,rlcuts,hdts,hdt2s):
    print fname
    chtmp=0
    nch=0
    hdt.SetTitle(fname)
    hdt.GetXaxis().SetTitle("Row counts (1 row = 240 ns)")
    hdt.GetYaxis().SetTitle("Counts / 1 row (240 ns)")
    hdt2.SetTitle(fname)
    hdt2.GetXaxis().SetTitle("Row counts (1 row = 240 ns)")
    hdt2.GetYaxis().SetTitle("Energy (eV)")
    for e in t:
        ch=e.ch
        if chtmp<ch:# ch is changed
            nch+=1
            print "ch%d"%(ch)
            chtmp=ch
        if e.beam==1 and e.good==1 and e.primary==1:
            if e.energy>6150 and e.energy<6600:
                hdt.Fill(e.row_next_extrig_nrp)
                hdt2.Fill(e.row_next_extrig_nrp,e.energy)
    c1 = ROOT.TCanvas("c1_%s"%fname,"c1_%s"%fname)
    c1.cd()
    hdt.DrawCopy("hist")
    c1.Update()
    c1.SaveAs("beamtiming_%s.pdf"%fname)
    c1.Close()
    c2 = ROOT.TCanvas("c2_%s"%fname,"c2_%s"%fname)
    c2.cd()
    hdt2.DrawCopy("cont4")
    c2.Update()
    c2.SaveAs("beamtiming_vs_energy_%s.png"%fname)
    c2.Close()

foutname = datadir+"run0397_beamtiming_%s_hist.root"%(template)
fout = ROOT.TFile.Open(foutname,"recreate")
for hdt,hdt2 in zip(hdts,hdt2):
    hdt.Write()
    hdt2.Write()
fout.Close()
print foutname, " is created."
