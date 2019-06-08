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

addon  = add + "_spillon_sprmcon_jbrscon_prime"
addoff = add + "_spilloff_sprmcon_jbrscon_prime"
fnameon="run%04d/run%04d_noi%04d_mass_2018%s"%(run,run,noi,addon)
fnameoff="run%04d/run%04d_noi%04d_mass_2018%s"%(run,run,noi,addoff)

# --------------------------------------------------------

ROOT.gStyle.SetOptStat(0)
foff = ROOT.TFile.Open(datadir+fnameoff+".root")
fon  = ROOT.TFile.Open(datadir+fnameon+".root")

# --------------------------------------------------------


nbin=200
minx=-100.
maxx=100.
h_sec_prm_spillon  = ROOT.TH1F("h_sec_prm_spillon","h_sec_prm_spillon",nbin,minx,maxx)
h_sec_prm_spilloff = ROOT.TH1F("h_sec_prm_spilloff","h_sec_prm_spilloff",nbin,minx,maxx)

nbin=60
minx=6900.
maxx=6960.
sec_prm_cuts=["beam==0 && sec_pr_mean!=-1e6",
              "beam==1 && sec_pr_mean!=-1e6",
              "beam==1 && sec_pr_mean<-10",
              "beam==1 && sec_pr_mean>=-10 && sec_pr_mean<-3",
              "beam==1 && sec_pr_mean>=-3 && sec_pr_mean<3",
              "beam==1 && sec_pr_mean>=3 && sec_pr_mean<6",
              "beam==1 && sec_pr_mean>=6 && sec_pr_mean<10",
              "beam==1 && sec_pr_mean>=10"]
h_ene_phc = [ROOT.TH1F("h_ene_phc%d"%(j),"h_ene_phc%d"%(j),nbin,minx,maxx) for j in xrange(len(sec_prm_cuts))]


c1 = ROOT.TCanvas("c1","c1")
legend = ROOT.TLegend(0.7, 0.68, 0.90, 0.78)
legend.SetFillColor(0)
legend.SetTextSize(0.03)
foff.chanall.Draw("sec_pr_mean>>h_sec_prm_spilloff","good==1 && primary==1 && beam==0","GOFF")
fon.chanall.Draw("sec_pr_mean>>h_sec_prm_spillon","good==1 && primary==1 && beam==1","GOFF")
legend.AddEntry(h_sec_prm_spilloff, "spill off", "l")
legend.AddEntry(h_sec_prm_spillon, "spill on", "l")
h_sec_prm_spillon.SetLineColor(2)
h_sec_prm_spilloff.SetTitle(fnameon+" secondary peak region mean spill off/on")
h_sec_prm_spilloff.GetXaxis().SetTitle("sec_pr_mean [ch]")
h_sec_prm_spilloff.GetYaxis().SetTitle("Counts / ch")
h_sec_prm_spilloff.Draw("hist")
h_sec_prm_spillon.Draw("histsame")
ROOT.gPad.SetLogy()
legend.Draw("same")
c1.SaveAs(outdir+fnameon+"_sec_pr_mean.pdf")


c2 = ROOT.TCanvas("c2","c2")
legend = ROOT.TLegend(0.15, 0.70, 0.5, 0.85)
legend.SetFillColor(0)
legend.SetTextSize(0.02)
for i in xrange(len(sec_prm_cuts)):
    if i==0:
        foff.chanall.Draw("energy_phc>>h_ene_phc%d"%(i),"good==1 && primary==1 && %s"%(sec_prm_cuts[i]),"GOFF")
        h_ene_phc[i].SetLineColor(i+1)
        legend.AddEntry(h_ene_phc[i], "spill off, no cut", "l")
        h_ene_phc[i].SetTitle(fnameon+" energy_phc with sec_prm cuts")
        h_ene_phc[i].GetXaxis().SetTitle("Energy [eV]")
        h_ene_phc[i].GetYaxis().SetTitle("Counts / eV")
        h_ene_phc[i].DrawCopy("hist")
    elif i==1:
        fon.chanall.Draw("energy_phc>>h_ene_phc%d"%(i),"good==1 && primary==1 && %s"%(sec_prm_cuts[i]),"GOFF")
        h_ene_phc[i].SetLineColor(i+1)
        legend.AddEntry(h_ene_phc[i], "spill on, no cut", "l")
        h_ene_phc[i].DrawCopy("histsame")
    else:
        fon.chanall.Draw("energy_phc>>h_ene_phc%d"%(i),"good==1 && primary==1 && %s"%(sec_prm_cuts[i]),"GOFF")
        h_ene_phc[i].SetLineColor(i+1)
        legend.AddEntry(h_ene_phc[i], sec_prm_cuts[i], "l")
        h_ene_phc[i].DrawCopy("histsame")
legend.Draw("same")
c2.SaveAs(outdir+fnameon+"_energy_phc_sec_prm_cuts.pdf")


c3 = ROOT.TCanvas("c3","c3")
h_ene_phc[0].SetTitle(fnameon+" energy_phc with sec_prm cuts (normalized)")
h_ene_phc[0].GetXaxis().SetTitle("Energy [eV]")
h_ene_phc[0].GetYaxis().SetTitle("Counts / eV")
h_ene_phc[0].DrawNormalized("hist")
for i in xrange(len(sec_prm_cuts)):
    if i==0: continue
    else:
        h_ene_phc[i].DrawNormalized("histsame")
legend.Draw("same")
c3.SaveAs(outdir+fnameon+"_energy_phc_sec_prm_cuts_norm.pdf")



nbin=60
minx=6900.
maxx=6960.
h_ene_phc_g = ROOT.TH1F("h_ene_phc_g","h_ene_phc_g",nbin,minx,maxx)
h_ene_phc_b = ROOT.TH1F("h_ene_phc_b","h_ene_phc_b",nbin,minx,maxx)

c4 = ROOT.TCanvas("c4","c4")
legend2 = ROOT.TLegend(0.15, 0.80, 0.5, 0.85)
legend2.SetFillColor(0)
legend2.SetTextSize(0.02)
fon.chanall.Draw("energy_phc>>h_ene_phc_g","good==1 && primary==1 && beam==1 && sec_pr_mean>-10 && sec_pr_mean<=10","GOFF")
fon.chanall.Draw("energy_phc>>h_ene_phc_b","good==1 && primary==1 && beam==1 && !(sec_pr_mean>-10 && sec_pr_mean<=10)","GOFF")
h_ene_phc_g.SetLineColor(2)
legend2.AddEntry(h_ene_phc_g, "sec_pr_mean>-10 && sec_pr_mean<=10", "l")
legend2.AddEntry(h_ene_phc_b, "other", "l")
h_ene_phc_g.SetTitle(fnameon+" energy_phc with sec_prm cuts")
h_ene_phc_g.GetXaxis().SetTitle("Energy [eV]")
h_ene_phc_g.GetYaxis().SetTitle("Counts / eV")
h_ene_phc_g.DrawCopy("hist")
h_ene_phc_b.DrawCopy("histsame")
legend2.Draw("same")
c4.SaveAs(outdir+fnameon+"_energy_phc_sec_prm_cuts_goodbad.pdf")



h_ene_vs_prm = ROOT.TH2F("h_ene_vs_prm","h_ene_vs_prm",1000,5000,10000,30,-50,70)
c5 = ROOT.TCanvas("c5","c5")
fon.chanall.Draw("sec_pr_mean:energy_phc>>h_ene_vs_prm","good==1 && primary==1 && beam==1","GOFF")
h_ene_vs_prm.SetTitle(fnameon+" energy_phc vs sec_prm spill on")
h_ene_vs_prm.GetXaxis().SetTitle("Energy [eV]")
h_ene_vs_prm.GetYaxis().SetTitle("sec_pr_sum [ch]")
h_ene_vs_prm.DrawCopy("cont4")
ROOT.gPad.SetLogz()
c5.SaveAs(outdir+fnameon+"_energy_phc_vs_sec_prm.png")


nbin=4000
minx=-2000.
maxx=2000.
h_jbr_spillon  = ROOT.TH1F("h_jbr_spillon","h_jbr_spillon",nbin,minx,maxx)
h_jbr_spilloff = ROOT.TH1F("h_jbr_spilloff","h_jbr_spilloff",nbin,minx,maxx)
c6 = ROOT.TCanvas("c6","c6")
legend = ROOT.TLegend(0.7, 0.68, 0.90, 0.78)
legend.SetFillColor(0)
legend.SetTextSize(0.03)
foff.chanall.Draw("jbr_region_sum>>h_jbr_spilloff","good==1 && primary==1 && beam==0","GOFF")
fon.chanall.Draw("jbr_region_sum>>h_jbr_spillon","good==1 && primary==1 && beam==1","GOFF")
legend.AddEntry(h_jbr_spilloff, "spill off", "l")
legend.AddEntry(h_jbr_spillon, "spill on", "l")
h_jbr_spillon.SetLineColor(2)
h_jbr_spilloff.SetTitle(fnameon+" just before rising region sum spill off/on")
h_jbr_spilloff.GetXaxis().SetTitle("jbr_region_sum [ch]")
h_jbr_spilloff.GetYaxis().SetTitle("Counts / 1 ch")
h_jbr_spilloff.Draw("hist")
h_jbr_spillon.Draw("histsame")
ROOT.gPad.SetLogy()
legend.Draw("same")
c6.SaveAs(outdir+fnameon+"_jbr_sum.pdf")



nbin=60
minx=6900.
maxx=6960.
jbr_cuts=["beam==0 && jbr_region_sum!=-1e6",
          "beam==1 && jbr_region_sum!=-1e6",
          "beam==1 && jbr_region_sum<-400",
          "beam==1 && jbr_region_sum>=-400 && jbr_region_sum<-200",
          "beam==1 && jbr_region_sum>=-200 && jbr_region_sum<0",
          "beam==1 && jbr_region_sum>=0 && jbr_region_sum<350",
          "beam==1 && jbr_region_sum>=350 && jbr_region_sum<400",
          "beam==1 && jbr_region_sum>=400"]
h_ene_phc_jbr = [ROOT.TH1F("h_ene_phc_jbr%d"%(j),"h_ene_phc_jbr%d"%(j),nbin,minx,maxx) for j in xrange(len(jbr_cuts))]


c7 = ROOT.TCanvas("c7","c7")
legend = ROOT.TLegend(0.15, 0.70, 0.5, 0.85)
legend.SetFillColor(0)
legend.SetTextSize(0.02)
for i in xrange(len(jbr_cuts)):
    if i==0:
        foff.chanall.Draw("energy_phc>>h_ene_phc_jbr%d"%(i),"good==1 && primary==1 && %s"%(jbr_cuts[i]),"GOFF")
        h_ene_phc_jbr[i].SetLineColor(i+1)
        legend.AddEntry(h_ene_phc_jbr[i], "spill off, no cut", "l")
        h_ene_phc_jbr[i].SetTitle(fnameon+" energy_phc with jbr cuts")
        h_ene_phc_jbr[i].GetXaxis().SetTitle("Energy [eV]")
        h_ene_phc_jbr[i].GetYaxis().SetTitle("Counts / eV")
        h_ene_phc_jbr[i].DrawCopy("hist")
    elif i==1:
        fon.chanall.Draw("energy_phc>>h_ene_phc_jbr%d"%(i),"good==1 && primary==1 && %s"%(jbr_cuts[i]),"GOFF")
        h_ene_phc_jbr[i].SetLineColor(i+1)
        legend.AddEntry(h_ene_phc_jbr[i], "spill on, no cut", "l")
        h_ene_phc_jbr[i].DrawCopy("histsame")
    else:
        fon.chanall.Draw("energy_phc>>h_ene_phc_jbr%d"%(i),"good==1 && primary==1 && %s"%(jbr_cuts[i]),"GOFF")
        h_ene_phc_jbr[i].SetLineColor(i+1)
        legend.AddEntry(h_ene_phc_jbr[i], jbr_cuts[i], "l")
        h_ene_phc_jbr[i].DrawCopy("histsame")
legend.Draw("same")
c7.SaveAs(outdir+fnameon+"_energy_phc_jbr_cuts.pdf")


c8 = ROOT.TCanvas("c8","c8")
h_ene_phc_jbr[0].SetTitle(fnameon+" energy_phc with jbr cuts (normalized)")
h_ene_phc_jbr[0].GetXaxis().SetTitle("Energy [eV]")
h_ene_phc_jbr[0].GetYaxis().SetTitle("Counts / eV")
h_ene_phc_jbr[0].DrawNormalized("hist")
for i in xrange(len(jbr_cuts)):
    if i==0: continue
    else:
        h_ene_phc_jbr[i].DrawNormalized("histsame")
legend.Draw("same")
c8.SaveAs(outdir+fnameon+"_energy_phc_jbr_cuts_norm.pdf")



nbin=60
minx=6900.
maxx=6960.
h_ene_phc_jbr_g = ROOT.TH1F("h_ene_phc_jbr_g","h_ene_phc_jbr_g",nbin,minx,maxx)
h_ene_phc_jbr_b = ROOT.TH1F("h_ene_phc_jbr_b","h_ene_phc_jbr_b",nbin,minx,maxx)

c9 = ROOT.TCanvas("c9","c9")
legend2 = ROOT.TLegend(0.2, 0.80, 0.5, 0.85)
legend2.SetFillColor(0)
legend2.SetTextSize(0.02)
fon.chanall.Draw("energy_phc>>h_ene_phc_jbr_g","good==1 && primary==1 && beam==1 && jbr_region_sum>-400 && jbr_region_sum<=350","GOFF")
fon.chanall.Draw("energy_phc>>h_ene_phc_jbr_b","good==1 && primary==1 && beam==1 && !(jbr_region_sum>-400 && jbr_region_sum<=350)","GOFF")
h_ene_phc_jbr_g.SetLineColor(2)
legend2.AddEntry(h_ene_phc_jbr_g, "jbr_region_sum>-400 && jbr_region_sum<=350", "l")
legend2.AddEntry(h_ene_phc_jbr_b, "other", "l")
h_ene_phc_jbr_g.SetTitle(fnameon+" energy_phc with jbr cuts")
h_ene_phc_jbr_g.GetXaxis().SetTitle("Energy [eV]")
h_ene_phc_jbr_g.GetYaxis().SetTitle("Counts / eV")
h_ene_phc_jbr_g.DrawCopy("hist")
h_ene_phc_jbr_b.DrawCopy("histsame")
legend2.Draw("same")
c9.SaveAs(outdir+fnameon+"_energy_phc_jbr_cuts_goodbad.pdf")



h_ene_vs_jbr = ROOT.TH2F("h_ene_vs_jbr","h_ene_vs_jbr",1000,5000,10000,400,-2000,2000)
c10 = ROOT.TCanvas("c10","c10")
fon.chanall.Draw("jbr_region_sum:energy_phc>>h_ene_vs_jbr","good==1 && primary==1 && beam==1","GOFF")
h_ene_vs_jbr.SetTitle(fnameon+" energy_phc vs jbr spill on")
h_ene_vs_jbr.GetXaxis().SetTitle("Energy [eV]")
h_ene_vs_jbr.GetYaxis().SetTitle("jbr_region_sum [ch]")
h_ene_vs_jbr.DrawCopy("cont4")
ROOT.gPad.SetLogz()
c10.SaveAs(outdir+fnameon+"_energy_phc_vs_jbr.png")


c11 = ROOT.TCanvas("c11","c11")
legend2 = ROOT.TLegend(0.2, 0.80, 0.5, 0.85)
legend2.SetFillColor(0)
legend2.SetTextSize(0.02)
fon.chanall.Draw("energy_phc>>h_ene_phc_g","good==1 && primary==1 && beam==1 && sec_pr_mean>-10 && sec_pr_mean<=10 && jbr_region_sum>-400 && jbr_region_sum<=350","GOFF")
fon.chanall.Draw("energy_phc>>h_ene_phc_b","good==1 && primary==1 && beam==1 && !(sec_pr_mean>-10 && sec_pr_mean<=10) && jbr_region_sum>-400 && jbr_region_sum<=350","GOFF")
h_ene_phc_g.SetLineColor(2)
legend2.AddEntry(h_ene_phc_g, "sec_pr_mean>-10 && sec_pr_mean<=10 && jbr_region_sum>-400 && jbr_region_sum<=350", "l")
legend2.AddEntry(h_ene_phc_b, "other && jbr_region_sum>-400 && jbr_region_sum<=350", "l")
h_ene_phc_g.SetTitle(fnameon+" energy_phc with sec_prm and jbr_resion_sum cuts")
h_ene_phc_g.GetXaxis().SetTitle("Energy [eV]")
h_ene_phc_g.GetYaxis().SetTitle("Counts / eV")
h_ene_phc_g.DrawCopy("hist")
h_ene_phc_b.DrawCopy("histsame")
legend2.Draw("same")
c11.SaveAs(outdir+fnameon+"_energy_phc_sec_prm_jbr_cuts_goodbad.pdf")


nbin=4000
minx=-2000.
maxx=2000.
h_pre_spillon  = ROOT.TH1F("h_pre_spillon","h_pre_spillon",nbin,minx,maxx)
h_pre_spilloff = ROOT.TH1F("h_pre_spilloff","h_pre_spilloff",nbin,minx,maxx)
c12 = ROOT.TCanvas("c12","c12")
legend = ROOT.TLegend(0.7, 0.68, 0.90, 0.78)
legend.SetFillColor(0)
legend.SetTextSize(0.03)
foff.chanall.Draw("pre_region_sum>>h_pre_spilloff","good==1 && primary==1 && beam==0","GOFF")
fon.chanall.Draw("pre_region_sum>>h_pre_spillon","good==1 && primary==1 && beam==1","GOFF")
legend.AddEntry(h_pre_spilloff, "spill off", "l")
legend.AddEntry(h_pre_spillon, "spill on", "l")
h_pre_spillon.SetLineColor(2)
h_pre_spilloff.SetTitle(fnameon+" pre region sum spill off/on")
h_pre_spilloff.GetXaxis().SetTitle("pre_region_sum [ch]")
h_pre_spilloff.GetYaxis().SetTitle("Counts / 1 ch")
h_pre_spilloff.Draw("hist")
h_pre_spillon.Draw("histsame")
ROOT.gPad.SetLogy()
legend.Draw("same")
c12.SaveAs(outdir+fnameon+"_pre_sum.pdf")



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
h_ene_phc_pre = [ROOT.TH1F("h_ene_phc_pre%d"%(j),"h_ene_phc_pre%d"%(j),nbin,minx,maxx) for j in xrange(len(pre_cuts))]


c13 = ROOT.TCanvas("c13","c13")
legend = ROOT.TLegend(0.15, 0.70, 0.5, 0.85)
legend.SetFillColor(0)
legend.SetTextSize(0.02)
for i in xrange(len(pre_cuts)):
    if i==0:
        foff.chanall.Draw("energy_phc>>h_ene_phc_pre%d"%(i),"good==1 && primary==1 && %s"%(pre_cuts[i]),"GOFF")
        h_ene_phc_pre[i].SetLineColor(i+1)
        legend.AddEntry(h_ene_phc_pre[i], "spill off, no cut", "l")
        h_ene_phc_pre[i].SetTitle(fnameon+" energy_phc with pre cuts")
        h_ene_phc_pre[i].GetXaxis().SetTitle("Energy [eV]")
        h_ene_phc_pre[i].GetYaxis().SetTitle("Counts / eV")
        h_ene_phc_pre[i].DrawCopy("hist")
    elif i==1:
        fon.chanall.Draw("energy_phc>>h_ene_phc_pre%d"%(i),"good==1 && primary==1 && %s"%(pre_cuts[i]),"GOFF")
        h_ene_phc_pre[i].SetLineColor(i+1)
        legend.AddEntry(h_ene_phc_pre[i], "spill on, no cut", "l")
        h_ene_phc_pre[i].DrawCopy("histsame")
    else:
        fon.chanall.Draw("energy_phc>>h_ene_phc_pre%d"%(i),"good==1 && primary==1 && %s"%(pre_cuts[i]),"GOFF")
        h_ene_phc_pre[i].SetLineColor(i+1)
        legend.AddEntry(h_ene_phc_pre[i], pre_cuts[i], "l")
        h_ene_phc_pre[i].DrawCopy("histsame")
legend.Draw("same")
c13.SaveAs(outdir+fnameon+"_energy_phc_pre_cuts.pdf")


c14 = ROOT.TCanvas("c14","c14")
h_ene_phc_pre[0].SetTitle(fnameon+" energy_phc with pre cuts (normalized)")
h_ene_phc_pre[0].GetXaxis().SetTitle("Energy [eV]")
h_ene_phc_pre[0].GetYaxis().SetTitle("Counts / eV")
h_ene_phc_pre[0].DrawNormalized("hist")
for i in xrange(len(pre_cuts)):
    if i==0: continue
    else:
        h_ene_phc_pre[i].DrawNormalized("histsame")
legend.Draw("same")
c14.SaveAs(outdir+fnameon+"_energy_phc_pre_cuts_norm.pdf")



nbin=60
minx=6900.
maxx=6960.
h_ene_phc_preg = ROOT.TH1F("h_ene_phc_preg","h_ene_phc_preg",nbin,minx,maxx)
h_ene_phc_preb = ROOT.TH1F("h_ene_phc_preb","h_ene_phc_preb",nbin,minx,maxx)

c15 = ROOT.TCanvas("c15","c15")
legend2 = ROOT.TLegend(0.2, 0.80, 0.5, 0.85)
legend2.SetFillColor(0)
legend2.SetTextSize(0.02)
fon.chanall.Draw("energy_phc>>h_ene_phc_preg","good==1 && primary==1 && beam==1 && pre_region_sum>-350","GOFF")
fon.chanall.Draw("energy_phc>>h_ene_phc_preb","good==1 && primary==1 && beam==1 && !(pre_region_sum>-350)","GOFF")
h_ene_phc_preg.SetLineColor(2)
legend2.AddEntry(h_ene_phc_preg, "pre_region_sum>-350", "l")
legend2.AddEntry(h_ene_phc_preb, "other", "l")
h_ene_phc_preg.SetTitle(fnameon+" energy_phc with pre cuts")
h_ene_phc_preg.GetXaxis().SetTitle("Energy [eV]")
h_ene_phc_preg.GetYaxis().SetTitle("Counts / eV")
h_ene_phc_preg.DrawCopy("hist")
h_ene_phc_preb.DrawCopy("histsame")
legend2.Draw("same")
c15.SaveAs(outdir+fnameon+"_energy_phc_pre_cuts_goodbad.pdf")



h_ene_vs_pre = ROOT.TH2F("h_ene_vs_pre","h_ene_vs_pre",1000,5000,10000,400,-2000,2000)
c16 = ROOT.TCanvas("c16","c16")
fon.chanall.Draw("pre_region_sum:energy_phc>>h_ene_vs_pre","good==1 && primary==1 && beam==1","GOFF")
h_ene_vs_pre.SetTitle(fnameon+" energy_phc vs pre spill on")
h_ene_vs_pre.GetXaxis().SetTitle("Energy [eV]")
h_ene_vs_pre.GetYaxis().SetTitle("pre_region_sum [ch]")
h_ene_vs_pre.DrawCopy("cont4")
ROOT.gPad.SetLogz()
c16.SaveAs(outdir+fnameon+"_energy_phc_vs_pre.png")



