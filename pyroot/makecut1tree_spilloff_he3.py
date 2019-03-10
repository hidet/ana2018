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
#runs=[277,278,279,281,282,
#      283,284,286,287,290,291,292,293,294,295,
#      296,297,298,299,300,301]

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
fnameins=["%s/run%04d/run%04d_noi%04d_mass_2018%s"%(util.dumprootdir,run,run,noi,add) for run,noi in zip(runs,util.get_noise_list(runs))]
for fnamein in fnameins:
    if os.path.isfile(fnamein+".root")==False:
        print "Error: file is missing %s"%(fnamein+".root")
        print "forgetting options? --pre=xxx --post=yyy"
        sys.exit(0)
print "runs: ", runs
print "basic cuts: ", cut_1
print "sprmc: ", sprmc
print "khet: ", ckheton
print "input files[0]: ", fnameins[0]+".root"
util.check_continue()
#print "debug exit"
#sys.exit(0)# please remove this exit

fins = [ROOT.TFile.Open(fnamein+".root","read") for fnamein in fnameins]
for j, (run,fin,fnamein) in enumerate(zip(runs,fins,fnameins)):
    # loop for runs
    if fin.GetListOfKeys().Contains(util.tree_name)==False:# no tree, skip this run
        print "run%d, no TTree name %s"%(run,util.tree_name)
        if fin.IsOpen(): fin.Close()
        continue
    print "%s was read"%(fnamein+".root")
    fin.cd()
    t = fin.Get(util.tree_name)
    # for speed up: create tmp file to save the tree with basic cut
    ftmpname=fnamein+"_cut1.root"
    ftmp = ROOT.TFile.Open(ftmpname,"recreate")
    newt = t.CopyTree(cut_1)# without cutting sprm
    print "%s has been created with %s"%(ftmpname,util.tree_name)
    ftmp.Close()
    fin.Close()
