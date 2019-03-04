import ROOT
import numpy as np
import pandas as pd
import sys
import os
import optparse

import pyroot_util as util
util=reload(util)



# --------------------------------------------------------
# --- files ---
# HE3: hss[0]
# HE4: hss[1]
# --- hnames ---
# [0] hene_phc_on       (spill on with basic cuts)
# [1] hene_phc_sprmon   (sec peak region mean cuts on)
# [2] hene_phc_sprmoff  (out of sec peak region mean cuts)
# [3] hene_phc_kheton   (with ext trig timing cut)
# [4] hene_phc_khetoff  (out of ext trig timing window)
# --------------------------------------------------------
target=["he3","he4"]
fnames=["hene_phc_sprmon160_301_he3_spline",
        "hene_phc_sprmon320_424_he4_spline"]
files = [ROOT.TFile.Open(util.rootdir+fname+".root") for fname in fnames]
lists = [f.GetListOfKeys() for f in files]
lnks  = [l.FirstLink() for l in lists]
hnamess=[]
for lnk in lnks:
    objn=[]
    while(lnk):
        objn.append(lnk.GetObject().GetName())
        lnk=lnk.Next()
    hnamess.append(objn)
hss=[]
for f,hns in zip(files,hnamess):
    hss.append([f.Get(hname) for hname in hns])

maxh=120.
l_khe4=ROOT.TLine(util.KHE4LA_K,0.,util.KHE4LA_K,maxh)
l_khe4.SetLineColor(3)
l_khe3=ROOT.TLine(util.KHE3LA,0.,util.KHE3LA,maxh)
l_khe3.SetLineColor(3)
l_pic=ROOT.TLine(util.PIC4TO3,0.,util.PIC4TO3,maxh)
l_pic.SetLineColor(5)
for j, (tgt,hs,hnames) in enumerate(zip(target,hss,hnamess)):
    nrebin=8# 1bin=2eV
    fscale=0.01
    leg = ROOT.TLegend(0.25, 0.70, 0.60, 0.85)
    leg.SetFillColor(0)
    leg.SetFillStyle(0);
    leg.SetTextSize(0.03)
    for i,(h,hname) in enumerate(zip(hs,hnames)):
        h.SetLineColor(i+1)
        h.SetLineWidth(1)
        h.Rebin(nrebin)
        if i==3:
            leg.AddEntry(h, hname, "e")
        else:
            h.Scale(fscale)
            leg.AddEntry(h, hname+" scale %.3f"%(fscale), "e")
    
    c1name="c1_%s"%(tgt)
    c1 = ROOT.TCanvas(c1name,c1name)
    c1.cd()
    ROOT.gStyle.SetOptStat(0)
    ROOT.gPad.SetGrid(0,0);
    ROOT.gPad.SetTicks();
    bw=hs[0].GetBinWidth(1)
    hs[0].GetXaxis().SetTitle("Energy [eV]")
    hs[0].GetXaxis().SetRangeUser(6000.,6700.)
    hs[0].GetYaxis().SetTitle("Counts / %.1f eV"%(bw))
    hs[0].GetYaxis().SetRangeUser(0.,maxh)
    hs[0].Draw("hist")
    for i in xrange(len(hs)-1):
        hs[i+1].Draw("histsame")
    leg.Draw("same")
    l_khe3.Draw("same")
    l_khe4.Draw("same")
    l_pic.Draw("same")    
    c1.Update()
    c1.SaveAs(util.figdir+hnames[0]+"_%s.pdf"%(tgt))

    c2name="c2_%s"%(tgt)
    c2 = ROOT.TCanvas(c2name,c2name)
    c2.cd()
    ROOT.gPad.SetGrid(0,0);
    ROOT.gPad.SetTicks();
    ROOT.gPad.SetLogy(1);
    hs[0].GetXaxis().SetRangeUser(5400.,7400.)
    hs[0].GetYaxis().SetRangeUser(1,1e5)
    hs[0].Draw("hist")
    for i in xrange(len(hs)-1):
        hs[i+1].Draw("histsame")
    leg.Draw("same")
    c2.Update()
    c2.SaveAs(util.figdir+hnames[0]+"_%s_log.pdf"%(tgt))


# compares khet histograms
c3name="c3_comp"
c3 = ROOT.TCanvas(c3name,c3name)
c3.cd()
ROOT.gPad.SetGrid(0,0);
ROOT.gPad.SetTicks();
leg2 = ROOT.TLegend(0.25, 0.70, 0.60, 0.85)
hss[0][3].SetLineColor(4)# he3
hss[1][3].SetLineColor(2)# he4
leg2.SetFillColor(0)
leg2.SetFillStyle(0);
leg2.SetTextSize(0.03)
leg2.AddEntry(hss[0][3], hnamess[0][3], "e")
leg2.AddEntry(hss[1][3], hnamess[1][3], "e")
bw=hss[0][3].GetBinWidth(1)
hss[0][3].GetXaxis().SetTitle("Energy [eV]")
hss[0][3].GetYaxis().SetTitle("Counts / %.1f eV"%(bw))
hss[0][3].GetXaxis().SetRangeUser(6000.,6700.)
hss[0][3].GetYaxis().SetRangeUser(0,maxh)
hss[0][3].Draw("hist")
hss[1][3].Draw("histsame")
leg2.Draw("same")
l_khe3.Draw("same")
l_khe4.Draw("same")
l_pic.Draw("same")
c3.Update()
c3.SaveAs(util.figdir+hnames[0]+"_comp.pdf")
