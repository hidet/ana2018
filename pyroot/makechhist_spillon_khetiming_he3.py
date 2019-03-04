from ROOT import gROOT
gROOT.SetBatch()

import ROOT
import numpy as np
import sys
import os

import pyroot_util as util
util= reload(util)


runs=util.runs_he3
#runs=util.runs_he4
#runs=[320,321]# for test

# -- spill on / off --
#beamflag = 0#      0:spill off, 1:spill on
#htag     = "off"#  used for hist name and fout name
beamflag = 1#     0:spill off, 1:spill on
htag     = "on"#  used for hist name and fout name
# -- cut parameters --
sprm_m   =  -10
sprm_p   =   10
jbrs_m   = -400
jbrs_p   =  350
row_next_thl=70# beam timing cut low
row_next_thh=90# beam timing cut high
ene_thl = 6000
ene_thh = 6700
# --------------------------------------------------------
options,args = util.parser.parse_args()
cut_pre  = options.cut_pre
cut_post = options.cut_post
add=""
if cut_pre!=0 or cut_post!=0: add+="_pre%03d_post%03d"%(cut_pre,cut_post)
if beamflag==0:   add+="_spilloff_sprmcon_jbrscon_prime"# root file name to be read
elif beamflag==1: add+="_spillon_sprmcon_jbrscon_prime"# root file name to be read
good    = "good==1"
prime   = "primary==1"
spill   = "beam==%d"%(beamflag)
sprmc   = "sec_pr_mean>%d && sec_pr_mean<%d"%(sprm_m,sprm_p)
jbrsc   = "jbr_region_sum>%d && jbr_region_sum<%d"%(jbrs_m,jbrs_p)
ckheton  = "(row_next_extrig_nrp>%d && row_next_extrig_nrp<%d)"%(row_next_thl,row_next_thh)
ckhetoff = "(row_next_extrig_nrp>%d || row_next_extrig_nrp<%d)"%(row_next_thh,row_next_thl)
cene     = "energy_phc>%d && energy_phc<%d"%(ene_thl,ene_thh)
cut_1           = "%s && %s && %s && %s"%(good,prime,spill,jbrsc)
cut_all         = "%s && %s && %s && %s && %s"%(good,prime,spill,sprmc,jbrsc)
cut_ene_all     = "%s && %s"%(cene,cut_all)
cut_kheton_all  = "%s && %s"%(ckheton,cut_all)
cut_khetoff_all = "%s && %s"%(ckhetoff,cut_all)
# --------------------------------------------------------
fnameins=["run%04d/run%04d_noi%04d_mass_2018%s"%(run,run,noi,add) for run,noi in zip(runs,util.get_noise_list(runs))]
fnameout="%s_%s_run%04d_%04d"%(util.hpht_phc,htag,runs[0],runs[-1])
for fnamein in fnameins:
    if os.path.isfile(util.datadir+fnamein+".root")==False:
        print "Error: file is missing %s"%(util.datadir+fnamein+".root")
        print "forgetting options? --pre=xxx --post=yyy"
        sys.exit(0)
print "runs: ", runs
print "basic cuts: ", cut_1
print "sprmc: ", sprmc
print "khet: ", ckheton
print "output file: ", util.outdir+fnameout+".root"
util.check_continue()
util.backup_rootfile(util.outdir+fnameout+".root")
print "debug exit"
sys.exit(0)# please remove this exit
# --------------------------------------------------------
chans = np.arange(480)[1::2]
titles = ["%d_ch%d"%(run,ch) for run in runs for ch in chans]
# pulse height histograms for calib
nbin = 25000
minx = 0.
maxx = 25000.
houts = [ROOT.TH1F("%s_%s%d_ch%d"%(util.hpht_phc,htag,run,ch),
                   "%s_%s%d_ch%d"%(util.hpht_phc,htag,run,ch),
                   nbin,minx,maxx) for run in runs for ch in chans]
houts_sprmoffs = [ROOT.TH1F("%s_sprmoff%d_ch%d"%(util.hpht_phc,run,ch),
                            "%s_sprmoff%d_ch%d"%(util.hpht_phc,run,ch),
                            nbin,minx,maxx) for run in runs for ch in chans]
houts_sprmons = [ROOT.TH1F("%s_sprmon%d_ch%d"%(util.hpht_phc,run,ch),
                           "%s_sprmon%d_ch%d"%(util.hpht_phc,run,ch),
                           nbin,minx,maxx) for run in runs for ch in chans]
# pulse height histograms with timing cuts
hout_khetoffs = [ROOT.TH1F("%s_khetoff%d_ch%d"%(util.hpht_phc,run,ch),
                           "%s_khetoff%d_ch%d"%(util.hpht_phc,run,ch),
                           nbin,minx,maxx) for run in runs for ch in chans]
hout_khetons  = [ROOT.TH1F("%s_kheton%d_ch%d"%(util.hpht_phc,run,ch),
                           "%s_kheton%d_ch%d"%(util.hpht_phc,run,ch),
                           nbin,minx,maxx) for run in runs for ch in chans]
# timing histograms
hkhets = [ROOT.TH1F("hkhet%d"%(run),"hkhet%d"%(run),300,-150.,150.) for run in runs]# 1bin=1ch=240ns
hkhet_sum = ROOT.TH1F("hkhet_sum","hkhet_sum",300,-150.,150.)
hkhet_ene_sum = ROOT.TH2F("hkhet_ene_sum","hkhet_ene_sum",150,-150.,150.,140,ene_thl,ene_thh)
fins = [ROOT.TFile.Open(util.datadir+fnamein+".root","read") for fnamein in fnameins]
fout = ROOT.TFile.Open(util.outdir+fnameout+".root","recreate")
print "%s is opened as recreate"%(util.outdir+fnameout+".root")
# --------------------------------------------------------
for j, (run,fin,fnamein) in enumerate(zip(runs,fins,fnameins)):
    # loop for runs
    if fin.GetListOfKeys().Contains(util.tree_name)==False:# no tree, skip this run
        print "run%d, no TTree name %s"%(run,util.tree_name)
        if fin.IsOpen(): fin.Close()
        continue
    if not fout.IsOpen():
        fout  = ROOT.TFile.Open(util.outdir+fnameout+".root","update")
        print "%s is re-opened as update"%(util.outdir+fnameout+".root")
    print "%s was read"%(util.datadir+fnamein+".root")
    chtmp=0
    fin.cd()
    t = fin.Get(util.tree_name)
    # for spped up: create tmp file to save the tree with basic cut
    ftmpname=util.datadir+fnamein+"_cut1.root"
    ftmp = ROOT.TFile.Open(ftmpname,"recreate")
    #newt = t.CopyTree(cut_all)
    newt = t.CopyTree(cut_1)# without cutting sprm
    print "%s has been created with %s"%(ftmpname,util.tree_name)
    fin.cd()
    for e in newt:
        # loop for events (from dump_root)
        ch=e.ch
        if chtmp==ch:
            print "run%d ch%d"%(run,ch)
            chtmp+=1
        elif chtmp<ch: chtmp=ch
        try:    i=titles.index("%d_ch%d"%(run,ch))
        except:
            print "Error: cannot find %d_ch%d "%(run,ch)
            break
        # --- cut conditions to fill ---
        houts[i].Fill(e.filt_value_phc)
        if e.sec_pr_mean<sprm_m or e.sec_pr_mean>sprm_p:
            houts_sprmoffs[i].Fill(e.filt_value_phc)
        elif e.sec_pr_mean>sprm_m and e.sec_pr_mean<sprm_p:
            houts_sprmons[i].Fill(e.filt_value_phc)
            if e.energy_phc>ene_thl and e.energy_phc<ene_thh:
                hkhets[j].Fill(-1.*e.row_next_extrig_nrp)
                hkhets[j].Fill(+1.*e.row_after_extrig_nrp)
                hkhet_sum.Fill(-1.*e.row_next_extrig_nrp)
                hkhet_sum.Fill(+1.*e.row_after_extrig_nrp)
                hkhet_ene_sum.Fill(-1.*e.row_next_extrig_nrp,e.energy_phc)
                hkhet_ene_sum.Fill(+1.*e.row_after_extrig_nrp,e.energy_phc)
            if e.row_next_extrig_nrp>row_next_thl and e.row_next_extrig_nrp<row_next_thh:
                hout_khetons[i].Fill(e.filt_value_phc)
            elif e.row_next_extrig_nrp<row_next_thl or e.row_next_extrig_nrp>row_next_thh:
                hout_khetoffs[i].Fill(e.filt_value_phc)
                        
    for chan in chans:
        i=titles.index("%d_ch%d"%(run,chan))
        if houts[i].GetEntries()==0: continue
        fout.cd()
        houts[i].Sumw2()
        houts[i].SetTitle("run%04d ch%d filt_value_phc, %s %s"%(run,chan,htag,cut_1))
        houts[i].GetXaxis().SetTitle("filt_value_phc [ch]")
        houts[i].GetYaxis().SetTitle("Counts / ch")
        houts[i].Write()
        houts_sprmoffs[i].Sumw2()
        houts_sprmoffs[i].SetTitle("run%04d ch%d filt_value_phc, sprmoff %s"%(run,chan,cut_1))
        houts_sprmoffs[i].GetXaxis().SetTitle("filt_value_phc [ch]")
        houts_sprmoffs[i].GetYaxis().SetTitle("Counts / ch")
        houts_sprmoffs[i].Write()
        houts_sprmons[i].Sumw2()
        houts_sprmons[i].SetTitle("run%04d ch%d filt_value_phc, sprmon %s"%(run,chan,cut_1))
        houts_sprmons[i].GetXaxis().SetTitle("filt_value_phc [ch]")
        houts_sprmons[i].GetYaxis().SetTitle("Counts / ch")
        houts_sprmons[i].Write()
        hout_khetons[i].Sumw2()
        hout_khetons[i].SetTitle("run%04d ch%d filt_value_phc, kheton %s"%(run,chan,cut_kheton_all))
        hout_khetons[i].GetXaxis().SetTitle("filt_value_phc [ch]")
        hout_khetons[i].GetYaxis().SetTitle("Counts / ch")
        hout_khetons[i].Write()
        hout_khetoffs[i].Sumw2()
        hout_khetoffs[i].SetTitle("run%04d ch%d filt_value_phc, khetoff %s"%(run,chan,cut_khetoff_all))
        hout_khetoffs[i].GetXaxis().SetTitle("filt_value_phc [ch]")
        hout_khetoffs[i].GetYaxis().SetTitle("Counts / ch")
        hout_khetoffs[i].Write()

    fout.cd()
    hkhets[j].Sumw2()
    hkhets[j].SetTitle("run%04d row_next_extrig_nrp, %s"%(run,cut_ene_all))
    hkhets[j].GetXaxis().SetTitle("Timing (1ch=240ns)")
    hkhets[j].GetYaxis().SetTitle("Counts / ch")
    hkhets[j].Write()
    
    if fin.IsOpen(): fin.Close()
    if ftmp.IsOpen(): ftmp.Close()

fout.cd()
hkhet_sum.Sumw2()
hkhet_sum.SetTitle("row_next_extrig_nrp, %s"%(cut_ene_all))
hkhet_sum.GetXaxis().SetTitle("Timing (1ch=240ns)")
hkhet_sum.GetYaxis().SetTitle("Counts / ch")
hkhet_sum.Write()

hkhet_ene_sum.SetTitle("timing vs energy_phc, %s"%(cut_ene_all))
hkhet_ene_sum.GetYaxis().SetTitle("Energy [eV]")
hkhet_ene_sum.GetXaxis().SetTitle("Timing (1ch=240ns)")
hkhet_ene_sum.Write()
    
if fout.IsOpen(): fout.Close("R")
print "%s has been created."%(util.outdir+fnameout+".root")
# --------------------------------------------------------
