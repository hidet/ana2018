from ROOT import gROOT
gROOT.SetBatch()

import ROOT
import numpy as np
import sys
import os

import pyroot_util as util
util= reload(util)


runs=util.runs_he3
#runs=[160,161]# for test
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
cut_1   = "%s && %s && %s && %s"%(good,prime,spill,jbrsc)
cut_all = "%s && %s && %s && %s && %s"%(good,prime,spill,sprmc,jbrsc)
# --------------------------------------------------------
# use cut_1 tree
fnameins=["%s/run%04d/run%04d_noi%04d_mass_2018%s_cut1"%(util.dumprootdir,run,run,noi,add) for run,noi in zip(runs,util.get_noise_list(runs))]
fnameout="%s/tree_%s_run%04d_%04d"%(util.outdir,htag,runs[0],runs[-1])
for fnamein in fnameins:
    if os.path.isfile(fnamein+".root")==False:
        print "Error: file is missing %s"%(fnamein+".root")
        print "forgetting options? --pre=xxx --post=yyy"
        sys.exit(0)
print "runs: ", runs
print "basic cuts: ", cut_1
print "sprmc: ", sprmc
print "fnameins[0]:", fnameins[0]+".root"
print "output file: ", fnameout+".root"
util.check_continue()
util.backup_rootfile(fnameout+".root")
# --------------------------------------------------------
fins = [ROOT.TFile.Open(fnamein+".root","read") for fnamein in fnameins]
fout = ROOT.TFile.Open(fnameout+".root","recreate")
print "%s is opened as recreate"%(fnameout+".root")
#tout = ROOT.TTree("tree%04d_%04d"%(runs[0],runs[-1]),"tree%04d_%04d"%(runs[0],runs[-1]))
tout = ROOT.TTree("tree","tree")
toutfill=tout.Fill
toutbranch=tout.Branch
bp_run             = np.zeros(1,dtype=np.intc)
bp_ch              = np.zeros(1,dtype=np.intc)
bp_ev              = np.zeros(1,dtype=np.intc)
bp_beam            = np.zeros(1,dtype=bool)
bp_ene             = np.zeros(1,dtype=np.float64)
bp_ene_zero        = np.zeros(1,dtype=np.float64)
bp_khet            = np.zeros(1,dtype=np.float64)
bp_filt_value      = np.zeros(1,dtype=np.float64)
bp_filt_value_dc   = np.zeros(1,dtype=np.float64)
bp_filt_value_phc  = np.zeros(1,dtype=np.float64)
bp_energy_phc      = np.zeros(1,dtype=np.float64)
bp_peak_time       = np.zeros(1,dtype=np.float64)
bp_peak_value      = np.zeros(1,dtype=np.float64)
bp_postpeak_deriv  = np.zeros(1,dtype=np.float64)
bp_pretrig_mean    = np.zeros(1,dtype=np.float64)
bp_pretrig_rms     = np.zeros(1,dtype=np.float64)
bp_pulse_average   = np.zeros(1,dtype=np.float64)
bp_rise_time       = np.zeros(1,dtype=np.float64)
bp_timestamp       = np.zeros(1,dtype=np.float64)
bp_rowmodp         = np.zeros(1,dtype=np.float64)
bp_pre_region_sum  = np.zeros(1,dtype=np.float64)
bp_jbr_region_sum  = np.zeros(1,dtype=np.float64)
bp_sec_pr_mean     = np.zeros(1,dtype=np.float64)
bp_sec_enemean     = np.zeros(1,dtype=np.float64)
toutbranch('run',            bp_run,             'run/I')
toutbranch('ch',             bp_ch,              'channel/I')
toutbranch('ev',             bp_ev,              'eventID/I')
toutbranch('beam',           bp_beam,            'beam/O')
toutbranch('ene',            bp_ene,             'energy/D')
toutbranch('ene_zero',       bp_ene_zero,        'energy_with_zero/D')
toutbranch('khet',           bp_khet,            'beam_time_diff/D')
toutbranch('filt_value',     bp_filt_value,      'filt_value/D')
toutbranch('filt_value_dc',  bp_filt_value_dc,   'filt_value_dc/D')
toutbranch('filt_value_phc', bp_filt_value_phc,  'filt_value_phc/D')
toutbranch('energy_phc',     bp_energy_phc,      'mass_energy/D')
toutbranch('peak_time',      bp_peak_time,       'peak_time/D')
toutbranch('peak_value',     bp_peak_value,      'peak_value/D')
toutbranch('postpeak_deriv', bp_postpeak_deriv,  'postpeak_deriv/D')
toutbranch('pretrig_mean',   bp_pretrig_mean,    'pretrig_mean/D')
toutbranch('pretrig_rms',    bp_pretrig_rms,     'pretrig_rms/D')
toutbranch('pulse_average',  bp_pulse_average,   'pulse_average/D')
toutbranch('rise_time',      bp_rise_time,       'rise_time/D')
toutbranch('timestamp',      bp_timestamp,       'timestamp/D')
toutbranch('rowmodp',        bp_rowmodp,         'rowmodifiedp/D')
toutbranch('pre_region_sum', bp_pre_region_sum,  'pre_region_sum/D')
toutbranch('jbr_region_sum', bp_jbr_region_sum,  'jbr_region_sum/D')
toutbranch('sec_pr_mean',    bp_sec_pr_mean,     'sec_peak_region_mean/D')
toutbranch('sec_enemean',    bp_sec_enemean,     'sec_enemean/D')
# --------------------------------------------------------
for j, (run,fin) in enumerate(zip(runs,fins)):
    # loop for runs
    if fin.GetListOfKeys().Contains(util.tree_name)==False:# no tree, skip this run
        print "run%d, no TTree name %s"%(run,util.tree_name)
        if fin.IsOpen(): fin.Close()
        continue
    # --------------------------------------------------------
    # calibration TSpline
    fcalib=""
    if beamflag==1: fcalib=util.pardir+"/calib/sprmon/calib_sprmon%d.root"%(run)
    else:           fcalib=util.pardir+"/calib/offsprmon/calib_offsprmon%d.root"%(run)
    if os.path.isfile(fcalib)==False:
        print "Error: file is missing %s"%(fcalib)
        sys.exit(0)
    fc=ROOT.TFile.Open(fcalib,"read")
    listfc=fc.GetListOfKeys();
    # --------------------------------------------------------
    # re-open fout
    if not fout.IsOpen():
        fout  = ROOT.TFile.Open(fnameout+".root","update")
        print "%s is re-opened as update"%(fnameout+".root")
    # --------------------------------------------------------
    print "%s was read"%(fin.GetName())
    chtmp=0
    fin.cd()
    t = fin.Get(util.tree_name)
    for e in t:# loop for events (from dump_root)
        ch=e.ch
        # --------------------------------------------------------
        if chtmp<ch:# ch is changed
            print "run%d ch%d"%(run,ch)
            chtmp=ch
            hname_calib=""
            if beamflag==1: hname_calib="hpht_phc_sprmon%d_ch%d"%(run,ch)
            else:           hname_calib="hpht_phc_offsprmon%d_ch%d"%(run,ch)
            sp3name="sp3calib_%s"%(hname_calib)
            sp3name_zero="sp3calib_%s_zero"%(hname_calib)
            try:
                sp3calib=fc.Get(sp3name)
                sp3calibEval=sp3calib.Eval
            except:
                print "cannot find %s"%(sp3name)
                sp3calib=None
            try:
                sp3calib_zero=fc.Get(sp3name_zero)
                sp3calib_zeroEval=sp3calib_zero.Eval
            except:
                print "cannot find %s"%(sp3name_zero)
                sp3calib_zero=None    
        # --------------------------------------------------------
        # energy calibration
        ene_phc=0.
        ene_phc_zero=0.
        if sp3calib is not None:      ene_phc=sp3calibEval(e.filt_value_phc)
        if sp3calib_zero is not None: ene_phc_zero=sp3calib_zeroEval(e.filt_value_phc)
        # --------------------------------------------------------
        # tree fill
        bp_run[0] =            run
        bp_ch[0] =             ch
        bp_ev[0] =             e.ev
        bp_beam[0] =           beamflag
        bp_ene[0] =            ene_phc
        bp_ene_zero[0] =       ene_phc_zero
        bp_khet[0] =           e.row_next_extrig_nrp
        bp_filt_value[0] =     e.filt_value
        bp_filt_value_dc[0] =  e.filt_value_dc
        bp_filt_value_phc[0] = e.filt_value_phc
        bp_energy_phc[0] =     e.energy_phc
        bp_peak_time[0] =      e.peak_time
        bp_peak_value[0] =     e.peak_value
        bp_postpeak_deriv[0] = e.postpeak_deriv
        bp_pretrig_mean[0] =   e.pretrig_mean
        bp_pretrig_rms[0] =    e.pretrig_rms
        bp_pulse_average[0] =  e.pulse_average
        bp_rise_time[0] =      e.rise_time
        bp_timestamp[0] =      e.timestamp
        bp_rowmodp[0] =        e.rowmodp
        bp_pre_region_sum[0] = e.pre_region_sum
        bp_jbr_region_sum[0] = e.jbr_region_sum
        bp_sec_pr_mean[0] =    e.sec_pr_mean
        bp_sec_enemean[0] =    e.sec_enemean
        toutfill()
                        
    if fin.IsOpen(): fin.Close()
    if fc.IsOpen(): fc.Close()

fout.cd()
tout.Write()
if fout.IsOpen(): fout.Close("R")
print "%s has been created."%(fnameout+".root")
# --------------------------------------------------------
