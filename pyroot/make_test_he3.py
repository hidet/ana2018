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
# use cut_1 tree
fnameins=["%s/run%04d/run%04d_noi%04d_mass_2018%s_cut1"%(util.dumprootdir,run,run,noi,add) for run,noi in zip(runs,util.get_noise_list(runs))]
fnameout="%s/test_%s_run%04d_%04d"%(util.outdir,htag,runs[0],runs[-1])
for fnamein in fnameins:
    if os.path.isfile(fnamein+".root")==False:
        print "Error: file is missing %s"%(fnamein+".root")
        print "forgetting options? --pre=xxx --post=yyy"
        sys.exit(0)
print "runs: ", runs
print "basic cuts: ", cut_1
print "sprmc: ", sprmc
print "khet: ", ckheton
print "output file: ", fnameout+".root"
util.check_continue()
util.backup_rootfile(fnameout+".root")
# --------------------------------------------------------
chans = np.arange(480)[1::2]
titles = ["%d_ch%d"%(run,ch) for run in runs for ch in chans]
# timing histograms
hkhet_sum = ROOT.TH1F("hkhet_sum","hkhet_sum",300,-150.,150.)
hkhet_ene_sum = ROOT.TH2F("hkhet_ene_sum","hkhet_ene_sum",150,-150.,150.,140,ene_thl,ene_thh)
fins = [ROOT.TFile.Open(fnamein+".root","read") for fnamein in fnameins]
fout = ROOT.TFile.Open(fnameout+".root","recreate")
print "%s is opened as recreate"%(fnameout+".root")
tout = ROOT.TTree("tree%04d_%04d"%(runs[0],runs[-1]),"tree%04d_%04d"%(runs[0],runs[-1]))
toutfill=tout.Fill
toutbranch=tout.Branch
bp_ch   = np.zeros(1,dtype=np.intc)
bp_ene  = np.zeros(1,dtype=np.float64)
bp_khet = np.zeros(1,dtype=np.float64)
bp_sprm = np.zeros(1,dtype=np.float64)
toutbranch('ch',          bp_ch,           'ch/I')
toutbranch('ene',         bp_ene,          'energy/D')
toutbranch('khet',        bp_khet,         'timediff/D')
toutbranch('sprm',        bp_sprm,         'sec_sprm/D')
# --------------------------------------------------------
for j, (run,fin) in enumerate(zip(runs,fins)):
    # loop for runs
    if fin.GetListOfKeys().Contains(util.tree_name)==False:# no tree, skip this run
        print "run%d, no TTree name %s"%(run,util.tree_name)
        if fin.IsOpen(): fin.Close()
        continue
    # calibration TSpline
    fcalib=util.pardir+"/calib/sprmon/calib_sprmon%d.root"%(run)
    if os.path.isfile(fcalib)==False:
        print "Error: file is missing %s"%(fcalib)
        sys.exit(0)
    fc=ROOT.TFile.Open(fcalib,"read")
    listfc=fc.GetListOfKeys();
    # re-open fout
    if not fout.IsOpen():
        fout  = ROOT.TFile.Open(fnameout+".root","update")
        print "%s is re-opened as update"%(fnameout+".root")
    print "%s was read"%(fin.GetName())
    
    chtmp=0
    fin.cd()
    t = fin.Get(util.tree_name)
    for e in t:
        # loop for events (from dump_root)
        ch=e.ch
        if chtmp<ch:# ch is changed
            print "run%d ch%d"%(run,ch)
            chtmp=ch
            hname_calib="hpht_phc_sprmon%d_ch%d"%(run,ch)
            sp3name="sp3calib_%s"%(hname_calib)
            try:
                sp3calib=fc.Get(sp3name)
                sp3calibEval=sp3calib.Eval
            except:
                print "cannot find %s"%(sp3name)
                sp3calib=None
        try:    i=titles.index("%d_ch%d"%(run,ch))
        except:
            print "Error: cannot find %d_ch%d "%(run,ch)
            break
        # energy calibration
        ene_phc=0.
        if sp3calib is not None: ene_phc=sp3calibEval(e.filt_value_phc)
        # tree fill
        if e.row_next_extrig_nrp>0. and e.row_next_extrig_nrp<200.:
            bp_ch[0]   = ch
            bp_ene[0]  = ene_phc
            bp_khet[0] = e.row_next_extrig_nrp
            bp_sprm[0] = e.sec_pr_mean
            toutfill()
        # --- cut conditions to fill ---
        if e.sec_pr_mean>sprm_m and e.sec_pr_mean<sprm_p:
            if ene_phc>ene_thl and ene_phc<ene_thh:
                hkhet_sum.Fill(-1.*e.row_next_extrig_nrp)
                hkhet_sum.Fill(+1.*e.row_after_extrig_nrp)
                hkhet_ene_sum.Fill(-1.*e.row_next_extrig_nrp,ene_phc)
                hkhet_ene_sum.Fill(+1.*e.row_after_extrig_nrp,ene_phc)
                        
    if fin.IsOpen(): fin.Close()

fout.cd()
tout.Write()
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
print "%s has been created."%(fnameout+".root")
# --------------------------------------------------------
