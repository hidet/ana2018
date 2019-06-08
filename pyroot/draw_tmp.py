import ROOT
import numpy as np
import pandas as pd
import sys
import os
import optparse

# --------------------------------------------------------
datadir="/Users/tatsuno/work/heates/data/TMU_2018U/"
outdir="./output/"

parser = optparse.OptionParser()
parser.add_option('--pre',  dest='cut_pre', action="store",type=int, help='set cut for pre samples',default=0)
parser.add_option('--post',  dest='cut_post', action="store",type=int, help='set cut for post samples',default=0)

options,args = parser.parse_args()
cut_pre        = options.cut_pre
cut_post       = options.cut_post


# --------------------------------------------------------


npar = len(args)
if (npar>=1):
    run_p = str(args[0])# 0001 or run0001
else:
    print "Error: specify run number of E62 ", args
    sys.exit(0)

RUNINFO="../csv/data_TMU_2018U.csv"
if os.path.exists(RUNINFO)==False: 
    print "%s is missing"%RUNINFO
    sys.exit(0)

df = pd.read_csv(RUNINFO)
run_list = df.iloc[:,0].tolist()
run=int(run_p)
ind = run_list.index(run)
noise_list = df.iloc[:,1].tolist()
noi = int(noise_list[ind])

# --------------------------------------------------------

outrun="%s/run%04d"%(outdir,run)
if not os.path.isdir(outrun):
    os.makedirs(outrun)    
add=""
if cut_pre!=0 or cut_post!=0:
    add="_pre%03d_post%03d"%(cut_pre,cut_post)
fname="run%04d/run%04d_noi%04d_mass_2018%s"%(run,run,noi,add)

# --------------------------------------------------------

ROOT.gStyle.SetOptStat(0)
f = ROOT.TFile.Open(datadir+fname+".root")

# --------------------------------------------------------

nbin=4000
minx=-2000.
maxx=2000.
h_pre_spillon  = ROOT.TH1F("h_pre_spillon","h_pre_spillon",nbin,minx,maxx)
h_pre_spilloff = ROOT.TH1F("h_pre_spilloff","h_pre_spilloff",nbin,minx,maxx)
c6 = ROOT.TCanvas("c6","c6")
legend = ROOT.TLegend(0.7, 0.68, 0.90, 0.78)
legend.SetFillColor(0)
legend.SetTextSize(0.03)
f.chanall.Draw("pre_region_sum>>h_pre_spilloff","good==1 && primary==1 && beam==0","GOFF")
f.chanall.Draw("pre_region_sum>>h_pre_spillon","good==1 && primary==1 && beam==1","GOFF")
legend.AddEntry(h_pre_spilloff, "spill off", "l")
legend.AddEntry(h_pre_spillon, "spill on", "l")
h_pre_spillon.SetLineColor(2)
h_pre_spilloff.SetTitle(fname+" pre region sum spill off/on")
h_pre_spilloff.GetXaxis().SetTitle("pre_region_sum [ch]")
h_pre_spilloff.GetYaxis().SetTitle("Counts / 1 ch")
h_pre_spilloff.Draw("hist")
h_pre_spillon.Draw("histsame")
ROOT.gPad.SetLogy()
legend.Draw("same")
c6.SaveAs(outdir+fname+"_pre_sum.pdf")



nbin=60
minx=6900.
maxx=6960.
pre_cuts=["beam==0 && pre_region_sum!=-1e6",
          "beam==1 && pre_region_sum!=-1e6",
          "beam==1 && pre_region_sum<-400",
          "beam==1 && pre_region_sum>=-400 && pre_region_sum<-350",
          "beam==1 && pre_region_sum>=-350 && pre_region_sum<0",
          "beam==1 && pre_region_sum>=0 && pre_region_sum<100",
          "beam==1 && pre_region_sum>=100 && pre_region_sum<200",
          "beam==1 && pre_region_sum>=200"]
h_ene_phc = [ROOT.TH1F("h_ene_phc%d"%(j),"h_ene_phc%d"%(j),nbin,minx,maxx) for j in xrange(len(pre_cuts))]


c2 = ROOT.TCanvas("c2","c2")
legend = ROOT.TLegend(0.15, 0.70, 0.5, 0.85)
legend.SetFillColor(0)
legend.SetTextSize(0.02)
for i in xrange(len(pre_cuts)):
    f.chanall.Draw("energy_phc>>h_ene_phc%d"%(i),"good==1 && primary==1 && %s"%(pre_cuts[i]),"GOFF")
    h_ene_phc[i].SetLineColor(i+1)
    if i==0:
        legend.AddEntry(h_ene_phc[i], "spill off, no cut", "l")
        h_ene_phc[i].SetTitle(fname+" energy_phc with pre cuts")
        h_ene_phc[i].GetXaxis().SetTitle("Energy [eV]")
        h_ene_phc[i].GetYaxis().SetTitle("Counts / eV")
        h_ene_phc[i].DrawCopy("hist")
    elif i==1:
        legend.AddEntry(h_ene_phc[i], "spill on, no cut", "l")
        h_ene_phc[i].DrawCopy("histsame")
    else:
        legend.AddEntry(h_ene_phc[i], pre_cuts[i], "l")
        h_ene_phc[i].DrawCopy("histsame")
legend.Draw("same")
c2.SaveAs(outdir+fname+"_energy_phc_pre_cuts.pdf")


c3 = ROOT.TCanvas("c3","c3")
h_ene_phc[0].SetTitle(fname+" energy_phc with pre cuts (normalized)")
h_ene_phc[0].GetXaxis().SetTitle("Energy [eV]")
h_ene_phc[0].GetYaxis().SetTitle("Counts / eV")
h_ene_phc[0].DrawNormalized("hist")
for i in xrange(len(pre_cuts)):
    if i==0: continue
    else:
        h_ene_phc[i].DrawNormalized("histsame")
legend.Draw("same")
c3.SaveAs(outdir+fname+"_energy_phc_pre_cuts_norm.pdf")



nbin=60
minx=6900.
maxx=6960.
h_ene_phc_g = ROOT.TH1F("h_ene_phc_g","h_ene_phc_g",nbin,minx,maxx)
h_ene_phc_b = ROOT.TH1F("h_ene_phc_b","h_ene_phc_b",nbin,minx,maxx)

c4 = ROOT.TCanvas("c4","c4")
legend2 = ROOT.TLegend(0.2, 0.80, 0.5, 0.85)
legend2.SetFillColor(0)
legend2.SetTextSize(0.02)
f.chanall.Draw("energy_phc>>h_ene_phc_g","good==1 && primary==1 && beam==1 && pre_region_sum>-350","GOFF")
f.chanall.Draw("energy_phc>>h_ene_phc_b","good==1 && primary==1 && beam==1 && !(pre_region_sum>-350)","GOFF")
h_ene_phc_g.SetLineColor(2)
legend2.AddEntry(h_ene_phc_g, "pre_region_sum>-350", "l")
legend2.AddEntry(h_ene_phc_b, "other", "l")
h_ene_phc_g.SetTitle(fname+" energy_phc with pre cuts")
h_ene_phc_g.GetXaxis().SetTitle("Energy [eV]")
h_ene_phc_g.GetYaxis().SetTitle("Counts / eV")
h_ene_phc_g.DrawCopy("hist")
h_ene_phc_b.DrawCopy("histsame")
legend2.Draw("same")
c4.SaveAs(outdir+fname+"_energy_phc_pre_cuts_goodbad.pdf")



h_ene_vs_pre = ROOT.TH2F("h_ene_vs_pre","h_ene_vs_jnr",1000,5000,10000,400,-2000,2000)
c5 = ROOT.TCanvas("c5","c5")
f.chanall.Draw("pre_region_sum:energy_phc>>h_ene_vs_pre","good==1 && primary==1 && beam==1","GOFF")
h_ene_vs_pre.SetTitle(fname+" energy_phc vs pre spill on")
h_ene_vs_pre.GetXaxis().SetTitle("Energy [eV]")
h_ene_vs_pre.GetYaxis().SetTitle("pre_region_sum [ch]")
h_ene_vs_pre.DrawCopy("cont4")
ROOT.gPad.SetLogz()
c5.SaveAs(outdir+fname+"_energy_phc_vs_pre.png")
