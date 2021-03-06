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

# -----------------------------------------------------------
spilltag="on"
runs=["0160_0301","0320_0424"]
fnames=[util.outdir+"/"+"tree_%s_run%s"%(spilltag,r) for r in runs]
tnames=["tree","tree"]
files = [ROOT.TFile.Open(fname+".root") for fname in fnames]
trees=[f.Get(tn) for f,tn in zip(files,tnames)]
# -----------------------------------------------------------
khetc="khet>=77 && khet<81"
# best timing window 77<=t<81
#sprmc="sec_pr_mean>-10 && sec_pr_mean<10"
start=-20
end=25
rebin=5# [eV]
#dc=10
dc=5
cs_l=np.arange(start,end,dc)
cs_h=np.arange(start+dc,end+dc,dc)
cs_l[0]=-1e6
cs_h[-1]=1e6
cuts=["sec_pr_mean>=%d && sec_pr_mean<%d && %s"%(cl,ch,khetc) for (cl,ch) in zip(cs_l,cs_h)]
# -----------------------------------------------------------
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
# -----------------------------------------------------------
maxh=2.*dc*rebin
l_khe4=ROOT.TLine(util.KHE4LA_K,0.,util.KHE4LA_K,maxh)
l_khe4.SetLineColor(3)
l_khe3=ROOT.TLine(util.KHE3LA,0.,util.KHE3LA,maxh)
l_khe3.SetLineColor(3)
be_he=24.6
l_khe4c=ROOT.TLine(util.KHE4LA_K-be_he,0.,util.KHE4LA_K-be_he,maxh)
l_khe4c.SetLineColor(5)
l_khe3c=ROOT.TLine(util.KHE3LA-be_he,0.,util.KHE3LA-be_he,maxh)
l_khe3c.SetLineColor(5)
# -----------------------------------------------------------
for tree,h in zip(trees,hs):
    for hi,ci in zip(h,cuts):
        print hi.GetName(), ci
        tree.Draw("ene>>"+hi.GetName(),"%s"%(ci),"GOFF")
# -----------------------------------------------------------
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
cvs[0].SaveAs(util.figdir+"/sprmcuts_run%s_%d.pdf"%(runs[0],rebin))



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
cvs[1].SaveAs(util.figdir+"/sprmcuts_run%s_%d.pdf"%(runs[1],rebin))
