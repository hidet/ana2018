from ROOT import gROOT
gROOT.SetBatch()

import ROOT
import numpy as np
import sys
import os
import array

import pyroot_util as util
util= reload(util)
   
# this class is insipired by
# https://root-forum.cern.ch/t/iteration-over-a-tree-in-pyroot-performance-issue/12264/7
class pytree:
    def __init__(self):
        self._branch={}

    def Init(self, tree, branches):
        for _name in branches:         
            br = tree.GetBranch(_name)     
            n = 1          
            if hasattr(tree, _name):       
                val = getattr(tree, _name) 
                cls_n = val.__class__.__name__           
                # GetBranch() doesn't work well with vector
                if cls_n.find("Py") != -1 and cls_n.find("Buffer") != -1:
                    n = len(val)           
            if br:
                self._branch[_name] = br
                type_char = br.GetListOfLeaves().At(0).GetTypeName()
                tp = util.type_convert_root_to_python(type_char)
                setattr(self, _name, array.array(tp, [0 for i in xrange(n)]))
                tree.SetBranchAddress(_name, getattr(self, _name))
        

def event_check2(ev,ene_thl,ene_thh):
    if ev.primary[0]!=1: return False
    elif ev.good[0]!=1: return False
    elif ev.beam[0]!=1: return False
    elif ev.energy[0]>ene_thl and ev.energy[0]<ene_thh: return True
    else: return False

runs=util.runs_he3
#runs=[160,161,162,163]
#runs=util.runs_he4


# -- cut parameters --
sprm_m   =  -10
sprm_p   =   10
jbrs_m   = -400
jbrs_p   =  350
row_next_thl=70# beam timing cut low
row_next_thh=90# beam timing cut high
ene_thl = 6150
ene_thh = 6600
# --------------------------------------------------------
options,args = util.parser.parse_args()
cut_pre  = options.cut_pre
cut_post = options.cut_post
add_spillon=""
add_spilloff=""
if cut_pre!=0 or cut_post!=0:
    add_spillon  += "_pre%03d_post%03d"%(cut_pre,cut_post)
    add_spilloff += "_pre%03d_post%03d"%(cut_pre,cut_post)
add_spillon  += "_spillon_sprmcon_jbrscon_prime"
add_spilloff += "_spilloff_sprmcon_jbrscon_prime"
good    = "good==1"
prime   = "primary==1"
beam    = "beam==1"
enecut  = "energy>%d && energy<%d"%(ene_thl,ene_thh)
cut_all = "%s && %s && %s && %s"%(good,prime,beam,enecut)
# --------------------------------------------------------
#tmp_dumprootdir="/Volumes/HEATES_HD/jparc2018_data/root"
fnameins_spillon=["%s/run%04d/run%04d_noi%04d_mass_2018%s"%(util.dumprootdir,run,run,noi,add_spillon) for run,noi in zip(runs,util.get_noise_list(runs))]
fnameins_spilloff=["%s/run%04d/run%04d_noi%04d_mass_2018%s"%(util.dumprootdir,run,run,noi,add_spilloff) for run,noi in zip(runs,util.get_noise_list(runs))]
for fnamein_on,fnamein_off in zip(fnameins_spillon,fnameins_spilloff):
    if os.path.isfile(fnamein_on+".root")==False or os.path.isfile(fnamein_off+".root")==False :
        print "Error: file is missing %s"%(fnamein_on+".root")
        print "forgetting options? --pre=150 --post=550"
        sys.exit(0)
print "runs: ", runs
print "cuts: ", cut_all
print "fnameins[0]:", fnameins_spillon[0]+".root"
util.check_continue()
# --------------------------------------------------------
fins_spillon  = [ROOT.TFile.Open(fnamein_on+".root","read") for fnamein_on in fnameins_spillon]
fins_spilloff = [ROOT.TFile.Open(fnamein_off+".root","read") for fnamein_off in fnameins_spilloff]

fout = ROOT.TFile.Open("beamtiming_run%04d_%04d.root"%(runs[0],runs[-1]),"recreate")

unit=6
run_tmp=[]
c=None
branches = ["good","beam","primary","energy","row_next_extrig_nrp"]
for j, (run,fin_on,fin_off) in enumerate(zip(runs,fins_spillon,fins_spilloff)):
    if fin_on.GetListOfKeys().Contains(util.tree_name)==False or fin_off.GetListOfKeys().Contains(util.tree_name)==False:
        print "run%d, no TTree name %s"%(run,util.tree_name)
        if fin_on.IsOpen(): fin_on.Close()
        if fin_off.IsOpen(): fin_off.Close()
        continue
            
    hon_name="h_spillon_%03d"%run
    hon = ROOT.TH1F(hon_name,hon_name,150,0,150)
    hon.SetTitle("run%03d spill on"%run)
    hon.GetXaxis().SetTitle("Row counts (1 row = 240 ns)")
    hon.GetYaxis().SetTitle("Counts / 1 row (240 ns)")
    honFill = hon.Fill
    hoff_name="h_spilloff_%03d"%run
    hoff = ROOT.TH1F(hoff_name,hoff_name,150,0,150)
    hoff.SetTitle("run%03d spill off"%run)
    hoff.GetXaxis().SetTitle("Row counts (1 row = 240 ns)")
    hoff.GetYaxis().SetTitle("Counts / 1 row (240 ns)")
    hoffFill = hoff.Fill

    fin_on.cd();  t_on = fin_on.Get(util.tree_name);
    fin_off.cd(); t_off = fin_off.Get(util.tree_name);
    
    ton = pytree();    ton.Init(t_on,branches);
    toff = pytree();   toff.Init(t_off,branches);

    print "%s is drawing..."%(hon.GetName())
    for entry in xrange(t_on.GetEntries()):
	for br in ton._branch.values():
	    br.GetEntry(entry)
        if event_check2(ton,ene_thl,ene_thh)==False: continue
        for br in toff._branch.values():
	    br.GetEntry(entry)
        #print entry, ton.row_next_extrig_nrp[0]
        honFill(ton.row_next_extrig_nrp[0])
        hoffFill(toff.row_next_extrig_nrp[0])
        
    ind=j%unit
    if ind==0:
        run_tmp=[]
        c = ROOT.TCanvas("c_%d"%run,"c_%d"%run)
        c.Divide(4,3)

    c.cd(ind*2+1)
    hon.DrawCopy("hist")
    c.cd(ind*2+2)
    hoff.DrawCopy("hist")

    run_tmp.append(run)
    if ind==unit-1 or run==runs[-1]:
        if c:
            c.Update()
            add=""
            for r in run_tmp: add+="_%d"%r
            c.SaveAs("beamtiming_runs"+add+".pdf")
            fout.cd()
            c.Write()
            c.Close()

    fout.cd()
    hon.Write()
    hoff.Write()
    t_on.Delete()
    t_off.Delete()
    hon.Delete()
    hoff.Delete()
        
    if fin_on.IsOpen(): fin_on.Close()
    if fin_off.IsOpen(): fin_off.Close()

fout.Close()

