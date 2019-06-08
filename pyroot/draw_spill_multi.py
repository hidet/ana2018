import ROOT
import numpy as np
import pandas as pd
import sys
import os
import optparse


runs = [381,382,383,384,385,386,387,389,
        395,396,397,398,399,400,401,402,403,
        406,407,408,410,411,421,422,423,424]

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

RUNINFO="../csv/data_TMU_2018U.csv"
if os.path.exists(RUNINFO)==False: 
    print "%s is missing"%RUNINFO
    sys.exit(0)

df = pd.read_csv(RUNINFO)
run_list = df.iloc[:,0].tolist()
noise_list = df.iloc[:,1].tolist()
inds=[run_list.index(run) for run in runs]
nois=[int(noise_list[ind]) for ind in inds]

# --------------------------------------------------------

add=""
if cut_pre!=0 or cut_post!=0:
    add="_pre%03d_post%03d"%(cut_pre,cut_post)

addon  = add + "_spillon_sprmcon_jbrscon_prime"
addoff = add + "_spilloff_sprmcon_jbrscon_prime"
fnameons=["run%04d/run%04d_noi%04d_mass_2018%s"%(run,run,noi,addon) for run, noi in zip(runs,nois)]
fnameoffs=["run%04d/run%04d_noi%04d_mass_2018%s"%(run,run,noi,addoff) for run, noi in zip(runs,nois)]

# --------------------------------------------------------

ROOT.gStyle.SetOptStat(0)
foffs = [ROOT.TFile.Open(datadir+fnameoff+".root") for fnameoff  in fnameoffs]
fons  = [ROOT.TFile.Open(datadir+fnameon+".root") for fnameon in fnameons]
fout = ROOT.TFile.Open(outdir+"he4_run%04d_%04d.root"%(runs[0],runs[-1]),"recreate")

nbin=15000
minx=0.
maxx=15000.
hene_dc_offs  = [ROOT.TH1F("hene_dc_off%d"%(i),"hene_dc_off%d"%(i),nbin,minx,maxx) for i in runs]
hene_phc_offs = [ROOT.TH1F("hene_phc_off%d"%(i),"hene_phc_off%d"%(i),nbin,minx,maxx) for i in runs]
hene_dc_ons   = [ROOT.TH1F("hene_dc_on%d"%(i),"hene_dc_on%d"%(i),nbin,minx,maxx) for i in runs]
hene_phc_ons  = [ROOT.TH1F("hene_phc_on%d"%(i),"hene_phc_on%d"%(i),nbin,minx,maxx) for i in runs]

good = "good==1"
prime = "primary==1"
beamon = "beam==1"
beamoff = "beam==0"
sprmc = "sec_pr_mean>-10 && sec_pr_mean<10"
jbrsc = "jbr_region_sum>-400 && jbr_region_sum<350"
cut_spilloff_all = "%s && %s && %s && %s && %s"%(good,prime,beamoff,sprmc,jbrsc)
cut_spillon_all = "%s && %s && %s && %s && %s"%(good,prime,beamon,sprmc,jbrsc)

# --------------------------------------------------------
for i, (run,foff,fon) in enumerate(zip(runs,foffs,fons)):
    foff.chanall.Draw("energy_dc>>hene_dc_off%d"%(run),"%s"%(cut_spilloff_all),"GOFF")
    foff.chanall.Draw("energy_phc>>hene_phc_off%d"%(run),"%s"%(cut_spilloff_all),"GOFF")
    fon.chanall.Draw("energy_dc>>hene_dc_on%d"%(run),"%s"%(cut_spillon_all),"GOFF")
    fon.chanall.Draw("energy_phc>>hene_phc_on%d"%(run),"%s"%(cut_spillon_all),"GOFF")

    hene_dc_offs[i].Sumw2()
    hene_phc_offs[i].Sumw2()
    hene_dc_ons[i].Sumw2()
    hene_phc_ons[i].Sumw2()
    
    hene_dc_offs[i].SetLineColor(i+1)
    hene_phc_offs[i].SetLineColor(i+1)
    hene_dc_ons[i].SetLineColor(i+1)
    hene_phc_ons[i].SetLineColor(i+1)

    hene_dc_offs[i].SetTitle("run%04d energy_dc, spill off %s"%(run,cut_spilloff_all))
    hene_dc_offs[i].GetXaxis().SetTitle("Energy [eV]")
    hene_dc_offs[i].GetYaxis().SetTitle("Counts / eV")
    hene_dc_ons[i].SetTitle("run%04d energy_dc, spill on %s"%(run,cut_spillon_all))
    hene_dc_ons[i].GetXaxis().SetTitle("Energy [eV]")
    hene_dc_ons[i].GetYaxis().SetTitle("Counts / eV")
    
    hene_phc_offs[i].SetTitle("run%04d energy_phc, spill off %s"%(run,cut_spilloff_all))
    hene_phc_offs[i].GetXaxis().SetTitle("Energy [eV]")
    hene_phc_offs[i].GetYaxis().SetTitle("Counts / eV")
    hene_phc_ons[i].SetTitle("run%04d energy_phc, spill on %s"%(run,cut_spillon_all))
    hene_phc_ons[i].GetXaxis().SetTitle("Energy [eV]")
    hene_phc_ons[i].GetYaxis().SetTitle("Counts / eV")

    # fit and save figures...
    fout.cd()
    hene_dc_offs[i].Write()
    hene_phc_offs[i].Write()
    hene_dc_ons[i].Write()
    hene_phc_ons[i].Write()


hene_dc_off_list = ROOT.TList()
hene_phc_off_list = ROOT.TList()
hene_dc_on_list = ROOT.TList()
hene_phc_on_list = ROOT.TList()

for i in xrange(len(runs)):
    hene_dc_off_list.Add(hene_dc_offs[i])
    hene_phc_off_list.Add(hene_phc_offs[i])
    hene_dc_on_list.Add(hene_dc_ons[i])
    hene_phc_on_list.Add(hene_phc_ons[i])

hene_dc_off_sum  = ROOT.TH1F("hene_dc_off_sum","hene_dc_off_sum",nbin,minx,maxx)
hene_phc_off_sum = ROOT.TH1F("hene_phc_off_sum","hene_phc_off_sum",nbin,minx,maxx)
hene_dc_on_sum   = ROOT.TH1F("hene_dc_on_sum","hene_dc_on_sum",nbin,minx,maxx)
hene_phc_on_sum  = ROOT.TH1F("hene_phc_on_sum","hene_phc_on_sum",nbin,minx,maxx)
hene_dc_off_sum.Merge(hene_dc_off_list)
hene_phc_off_sum.Merge(hene_phc_off_list)
hene_dc_on_sum.Merge(hene_dc_on_list)
hene_phc_on_sum.Merge(hene_phc_on_list)

hene_dc_off_sum.SetLineColor(1)
hene_dc_off_sum.SetTitle("energy_dc_sum, spill off %s"%(cut_spilloff_all))
hene_dc_off_sum.GetXaxis().SetTitle("Energy [eV]")
hene_dc_off_sum.GetYaxis().SetTitle("Counts / eV")

hene_dc_on_sum.SetLineColor(2)
hene_dc_on_sum.SetTitle("energy_dc_sum, spill on %s"%(cut_spillon_all))
hene_dc_on_sum.GetXaxis().SetTitle("Energy [eV]")
hene_dc_on_sum.GetYaxis().SetTitle("Counts / eV")

legend = ROOT.TLegend(0.7, 0.68, 0.90, 0.78)
legend.SetFillColor(0)
legend.SetTextSize(0.03)
legend.AddEntry(hene_dc_off_sum, "spill off", "l")
legend.AddEntry(hene_dc_on_sum, "spill on", "l")
# fit and save figures...
fout.cd()
hene_dc_off_sum.Write()
hene_dc_on_sum.Write()

hene_phc_off_sum.SetLineColor(1)
hene_phc_off_sum.SetTitle("energy_phc_sum, spill off %s"%(cut_spilloff_all))
hene_phc_off_sum.GetXaxis().SetTitle("Energy [eV]")
hene_phc_off_sum.GetYaxis().SetTitle("Counts / eV")

hene_phc_on_sum.SetLineColor(2)
hene_phc_on_sum.SetTitle("energy_phc_sum, spill on %s"%(cut_spillon_all))
hene_phc_on_sum.GetXaxis().SetTitle("Energy [eV]")
hene_phc_on_sum.GetYaxis().SetTitle("Counts / eV")

legend2 = ROOT.TLegend(0.7, 0.68, 0.90, 0.78)
legend2.SetFillColor(0)
legend2.SetTextSize(0.03)
legend2.AddEntry(hene_phc_off_sum, "spill off", "l")
legend2.AddEntry(hene_phc_on_sum, "spill on", "l")
# fit and save figures...
fout.cd()
hene_phc_off_sum.Write()
hene_phc_on_sum.Write()

c2 = ROOT.TCanvas("c2","c2")
hene_dc_off_sum.DrawCopy("hist")
hene_dc_on_sum.DrawCopy("histsame")
legend.Draw("same")
c2.SaveAs(outdir+"he4_energy_dc_sec_prm_cuts.pdf")


c3 = ROOT.TCanvas("c3","c3")
hene_phc_off_sum.DrawCopy("hist")
hene_phc_on_sum.DrawCopy("histsame")
legend2.Draw("same")
c3.SaveAs(outdir+"he4_energy_phc_sec_prm_cuts.pdf")


fout.Close()
