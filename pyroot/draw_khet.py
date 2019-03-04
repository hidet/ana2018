from ROOT import gROOT
gROOT.SetBatch()

import ROOT
import numpy as np
import pandas as pd
import sys
import os
import optparse

import pyroot_util as util
util=reload(util)

# best timing window 77<=t<81

fnames=["test_on_run0160_0301","test_on_run0320_0424"]
runs=["160_301","320_424"]
tnames=["tree%s"%(r) for r in runs]
files = [ROOT.TFile.Open(util.outdir+fname+".root") for fname in fnames]
trees=[f.Get(tn) for f,tn in zip(files,tnames)]

#sprmc="sprm>-10 && sprm<10"
sprmc="sprm>-5 && sprm<5"
cs1=[76,77,78,79,80,81,82,83]
cuts=["khet>=%d && khet<%d && %s"%(c,c+1,sprmc) for c in cs1]
cs2=[76,78,80,82]
cuts2=["khet>=%d && khet<%d && %s"%(c,c+2,sprmc) for c in cs2]

rebin=2# 1 eV
lx=6190.
hx=6260.
nbin=int((hx-lx)/rebin)
hhe3=[ROOT.TH1F("he3_cut%d_%d"%(i,rebin),
                "he3_cut%d_%d"%(i,rebin),nbin,lx,hx) for i in xrange(len(cuts))]
lx=6430.
hx=6500.
nbin=int((hx-lx)/rebin)
hhe4=[ROOT.TH1F("he4_cut%d_%d"%(i,rebin),
                "he4_cut%d_%d"%(i,rebin),nbin,lx,hx) for i in xrange(len(cuts))]

hs=[]
hs.append(hhe3)
hs.append(hhe4)

maxh=15.*rebin
l_khe4=ROOT.TLine(util.KHE4LA_K,0.,util.KHE4LA_K,maxh)
l_khe4.SetLineColor(3)
l_khe3=ROOT.TLine(util.KHE3LA,0.,util.KHE3LA,maxh)
l_khe3.SetLineColor(3)
be_he=24.6
l_khe4c=ROOT.TLine(util.KHE4LA_K-be_he,0.,util.KHE4LA_K-be_he,maxh)
l_khe4c.SetLineColor(5)
l_khe3c=ROOT.TLine(util.KHE3LA-be_he,0.,util.KHE3LA-be_he,maxh)
l_khe3c.SetLineColor(5)

for tree,h in zip(trees,hs):
    for hi,ci in zip(h,cuts):
        print hi.GetName(), ci
        tree.Draw("ene>>"+hi.GetName(),"%s"%(ci),"GOFF")


cnames=["c%s"%(r) for r in runs]
cvs = [ROOT.TCanvas(cn,cn) for cn in cnames]
cvs[0].cd()
ROOT.gStyle.SetOptStat(0)
ROOT.gPad.SetGrid(0,0);
ROOT.gPad.SetTicks();
bw=hhe3[0].GetBinWidth(1)
hhe3[0].GetXaxis().SetTitle("Energy [eV]")
hhe3[0].GetYaxis().SetTitle("Counts / %.1f eV"%(bw))
hhe3[0].GetYaxis().SetRangeUser(0.,maxh)
hhe3[0].SetTitle("He3 ene_phc with timing cuts")
hhe3[0].Draw("hist")
leg1 = ROOT.TLegend(0.55, 0.60, 0.85, 0.85)
leg1.SetFillColor(0)
leg1.SetFillStyle(0);
leg1.SetTextSize(0.02)
leg1.AddEntry(hhe3[0].GetName(), cuts[0], "e")
for i in xrange(len(hhe3)-1):
    hhe3[i+1].SetLineColor(i+2)
    hhe3[i+1].Draw("histsame")
    leg1.AddEntry(hhe3[i+1].GetName(), cuts[i+1], "e")
leg1.Draw("same")
l_khe3.Draw("same")
l_khe3c.Draw("same")
cvs[0].Update()
cvs[0].SaveAs(util.figdir+"khetcuts_run%s_%d.pdf"%(runs[0],rebin))



cvs[1].cd()
ROOT.gStyle.SetOptStat(0)
ROOT.gPad.SetGrid(0,0);
ROOT.gPad.SetTicks();
bw=hhe4[0].GetBinWidth(1)
hhe4[0].GetXaxis().SetTitle("Energy [eV]")
hhe4[0].GetYaxis().SetTitle("Counts / %.1f eV"%(bw))
hhe4[0].GetYaxis().SetRangeUser(0.,maxh)
hhe4[0].SetTitle("He4 ene_phc with timing cuts")
hhe4[0].Draw("hist")
leg2 = ROOT.TLegend(0.55, 0.60, 0.85, 0.85)
leg2.SetFillColor(0)
leg2.SetFillStyle(0);
leg2.SetTextSize(0.02)
leg2.AddEntry(hhe4[0].GetName(), cuts[0], "e")
for i in xrange(len(hhe4)-1):
    hhe4[i+1].SetLineColor(i+2)
    hhe4[i+1].Draw("histsame")
    leg2.AddEntry(hhe4[i+1].GetName(), cuts[i+1], "e")
leg2.Draw("same")
l_khe4.Draw("same")
l_khe4c.Draw("same")
cvs[1].Update()
cvs[1].SaveAs(util.figdir+"khetcuts_run%s_%d.pdf"%(runs[1],rebin))
