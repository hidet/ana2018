from ROOT import gROOT
gROOT.SetBatch()

import ROOT
import numpy as np
import sys
import os

import pyroot_util as util
util= reload(util)


#runs=util.runs_he3
runs=util.runs_he4
#runs=[320,321]# for test

# -- spill on / off --
beamflag = 0#      0:spill off, 1:spill on
htag     = "off"#  used for hist name and fout name
#beamflag = 1#     0:spill off, 1:spill on
#htag     = "on"#  used for hist name and fout name
# -- cut parameters --
sprm_m   =  -10
sprm_p   =   10
jbrs_m   = -400
jbrs_p   =  350
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
cut_1   = "%s && %s && %s && %s"%(good,prime,spill,jbrsc)
cut_all = "%s && %s && %s && %s && %s"%(good,prime,spill,sprmc,jbrsc)
# --------------------------------------------------------
fnameins=["%s/run%04d/run%04d_noi%04d_mass_2018%s"%(util.dumprootdir,run,run,noi,add) for run,noi in zip(runs,util.get_noise_list(runs))]
fnameout="%s/%s_%s_run%04d_%04d"%(util.outdir,util.hpht_phc,htag,runs[0],runs[-1])
for fnamein in fnameins:
    if os.path.isfile(fnamein+".root")==False:
        print "Error: file is missing %s"%(fnamein+".root")
        print "forgetting options? --pre=xxx --post=yyy"
        sys.exit(0)
print "runs: ", runs
print "cuts (cut_1): ", cut_1
print "sprmc: ", sprmc
print "input files[0]: ", fnameins[0]+".root"
print "output file:    ", fnameout+".root"
util.check_continue()
util.backup_rootfile(fnameout+".root")
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
fins = [ROOT.TFile.Open(fnamein+".root","read") for fnamein in fnameins]
fout = ROOT.TFile.Open(fnameout+".root","recreate")
print "%s is opened as recreate"%(fnameout+".root")
# --------------------------------------------------------
for j, (run,fin,fnamein) in enumerate(zip(runs,fins,fnameins)):
    # loop for runs
    if fin.GetListOfKeys().Contains(util.tree_name)==False:# no tree, skip this run
        print "run%d, no TTree name %s"%(run,util.tree_name)
        if fin.IsOpen(): fin.Close()
        continue
    if not fout.IsOpen():
        fout  = ROOT.TFile.Open(fnameout+".root","update")
        print "%s is re-opened as update"%(fnameout+".root")
    print "%s was read"%(fnamein+".root")
    chtmp=0
    fin.cd()
    t = fin.Get(util.tree_name)
    # for speed up: create tmp file to save the tree with basic cut
    ftmpname=fnamein+"_cut1.root"
    ftmp = ROOT.TFile.Open(ftmpname,"recreate")
    #newt = t.CopyTree(cut_all)
    newt = t.CopyTree(cut_1)# without cutting sprm
    print "%s has been created with %s"%(ftmpname,util.tree_name)
    fin.cd()
    for e in newt:
        # loop for events (from dump_root)
        ch=e.ch
        if chtmp<ch:# ch is changed
            print "run%d ch%d"%(run,ch)
            chtmp=ch
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

    if fin.IsOpen(): fin.Close()
    if ftmp.IsOpen(): ftmp.Close()

if fout.IsOpen(): fout.Close("R")
print "%s has been created."%(fnameout+".root")
# --------------------------------------------------------
