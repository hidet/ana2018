from ROOT import gROOT
gROOT.SetBatch()

import ROOT
import numpy as np
import sys
import os

import pyroot_util as util
util= reload(util)


# -- runs --
runs=util.runs_he3
#runs=util.runs_he4
#runs=[160,161]
#runs=[222]


# --------------------------------------------------------
options,args = util.parser.parse_args()
cut_pre  = options.cut_pre
cut_post = options.cut_post
spill = options.spill

add=""
if cut_pre!=0 or cut_post!=0: add+="_pre%03d_post%03d"%(cut_pre,cut_post)
if spill==-1:
    print "Error: please specify --spill=0 or --spill=1"
    sys.exit(0)
elif spill==0: add+="_spilloff_sprmcon_jbrscon_prime"# root file name to be read
elif spill==1: add+="_spillon_sprmcon_jbrscon_prime"# root file name to be read
jbrs_m   = -400
jbrs_p   =  350
good    = "good==1"
prime   = "primary==1"
jbrsc   = "jbr_region_sum>%d && jbr_region_sum<%d"%(jbrs_m,jbrs_p)
cut_1   = "%s && %s && %s"%(good,prime,jbrsc)
# --------------------------------------------------------
fnameins=["%s/run%04d/run%04d_noi%04d_mass_2018%s"%(util.dumprootdir,run,run,noi,add) for run,noi in zip(runs,util.get_noise_list(runs))]
fnamecalibons=["%s/run%04d/run%04d_calib_onsprmon%s"%(util.dumprootdir,run,run,add) for run in runs]
fnamecaliboffs=["%s/run%04d/run%04d_calib_offsprmon%s"%(util.dumprootdir,run,run,add) for run in runs]
fnameout="%s/tree_run%04d_%04d%s"%(util.dumprootdir,runs[0],runs[-1],add)
for fnamein in fnameins:
    if os.path.isfile(fnamein+".root")==False:
        print "Error: file is missing %s"%(fnamein+".root")
        print "forgetting options? --pre=xxx --post=yyy"
        sys.exit(0)
print "runs: ", runs
print "basic cuts: ", cut_1
print "fnameins[0]:       ", fnameins[0]+".root"
print "fnamecalibons[0]:  ", fnamecalibons[0]+".root"
print "fnamecaliboffs[0]: ", fnamecaliboffs[0]+".root"
print "output file:       ", fnameout+".root"
util.check_continue()
util.backup_rootfile(fnameout+".root")
# --------------------------------------------------------
fins = [ROOT.TFile.Open(fnamein+".root","read") for fnamein in fnameins]
fout = ROOT.TFile.Open(fnameout+".root","recreate")
print "%s is opened as recreate"%(fnameout+".root")
tout = ROOT.TTree("tree","tree")
toutfill=tout.Fill
toutbranch=tout.Branch
bp_run             = np.zeros(1,dtype=np.intc)
bp_ch              = np.zeros(1,dtype=np.intc)
bp_ev              = np.zeros(1,dtype=np.intc)
bp_beam            = np.zeros(1,dtype=bool)
bp_sprmc           = np.zeros(1,dtype=bool)
bp_shift1          = np.zeros(1,dtype=bool)
bp_filt_phase      = np.zeros(1,dtype=np.float32)
bp_filt_value      = np.zeros(1,dtype=np.float32)
bp_filt_value_dc   = np.zeros(1,dtype=np.float32)
bp_filt_value_phc  = np.zeros(1,dtype=np.float32)
bp_energy          = np.zeros(1,dtype=np.float32)
bp_energy_dc       = np.zeros(1,dtype=np.float32)
bp_energy_phc      = np.zeros(1,dtype=np.float32)
bp_peak_time       = np.zeros(1,dtype=np.float32)
bp_peak_value      = np.zeros(1,dtype=np.float32)
bp_postpeak_deriv  = np.zeros(1,dtype=np.float32)
bp_pretrig_mean    = np.zeros(1,dtype=np.float32)
bp_pretrig_rms     = np.zeros(1,dtype=np.float32)
bp_pulse_average   = np.zeros(1,dtype=np.float32)
bp_rise_time       = np.zeros(1,dtype=np.float32)
bp_timestamp       = np.zeros(1,dtype=np.float64)
bp_rowp            = np.zeros(1,dtype=np.float64)
bp_rown            = np.zeros(1,dtype=np.float64)
bp_pre_region_sum  = np.zeros(1,dtype=np.float32)
bp_jbr_region_sum  = np.zeros(1,dtype=np.float32)
bp_sec_pr_mean     = np.zeros(1,dtype=np.float32)
bp_sec_pr_meanmean = np.zeros(1,dtype=np.float32)
bp_sec_pr_maxmean  = np.zeros(1,dtype=np.float32)
bp_sec_enemean     = np.zeros(1,dtype=np.float32)
bp_sec_enemax      = np.zeros(1,dtype=np.float32)
bp_ene             = np.zeros(1,dtype=np.float32)
bp_ene_zero        = np.zeros(1,dtype=np.float32)
bp_dt              = np.zeros(1,dtype=np.float32)
toutbranch('run',            bp_run,             'run/I')
toutbranch('ch',             bp_ch,              'ch/I')
toutbranch('ev',             bp_ev,              'ev/I')
toutbranch('beam',           bp_beam,            'beam/O')
toutbranch('sprmc',          bp_sprmc,           'sprmc/O')
toutbranch('shift1',         bp_shift1,          'shift1/O')
toutbranch('filt_phase',     bp_filt_phase,      'filt_phase/F')
toutbranch('filt_value',     bp_filt_value,      'filt_value/F')
toutbranch('filt_value_dc',  bp_filt_value_dc,   'filt_value_dc/F')
toutbranch('filt_value_phc', bp_filt_value_phc,  'filt_value_phc/F')
toutbranch('energy',         bp_energy,          'energy/F')
toutbranch('energy_dc',      bp_energy_dc,       'energy_dc/F')
toutbranch('energy_phc',     bp_energy_phc,      'energy_phc/F')
toutbranch('peak_time',      bp_peak_time,       'peak_time/F')
toutbranch('peak_value',     bp_peak_value,      'peak_value/F')
toutbranch('postpeak_deriv', bp_postpeak_deriv,  'postpeak_deriv/F')
toutbranch('pretrig_mean',   bp_pretrig_mean,    'pretrig_mean/F')
toutbranch('pretrig_rms',    bp_pretrig_rms,     'pretrig_rms/F')
toutbranch('pulse_average',  bp_pulse_average,   'pulse_average/F')
toutbranch('rise_time',      bp_rise_time,       'rise_time/F')
toutbranch('timestamp',      bp_timestamp,       'timestamp/D')
toutbranch('rowp',           bp_rowp,            'rowp/D')
toutbranch('rown',           bp_rown,            'rown/D')
toutbranch('pre_region_sum', bp_pre_region_sum,  'pre_region_sum/F')
toutbranch('jbr_region_sum', bp_jbr_region_sum,  'jbr_region_sum/F')
toutbranch('sec_pr_mean',    bp_sec_pr_mean,     'sec_pr_mean/F')
toutbranch('sec_pr_meanmean',bp_sec_pr_meanmean, 'sec_pr_meanmean/F')
toutbranch('sec_pr_maxmean', bp_sec_pr_maxmean,  'sec_pr_maxmean/F')
toutbranch('sec_enemean',    bp_sec_enemean,     'sec_enemean/F')
toutbranch('sec_enemax',     bp_sec_enemax,      'sec_enemax/F')
toutbranch('dt',             bp_dt,              'dt/F')
toutbranch('ene',            bp_ene,             'ene/F')
toutbranch('ene_zero',       bp_ene_zero,        'ene_zero/F')

# --------------------------------------------------------
for j, (run,fin,fnamecalibon,fnamecaliboff) in enumerate(zip(runs,fins,fnamecalibons,fnamecaliboffs)):
    # loop for runs
    if fin.GetListOfKeys().Contains(util.tree_name)==False:# no tree, skip this run
        print "run%d, no TTree name %s"%(run,util.tree_name)
        if fin.IsOpen(): fin.Close()
        continue
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
        if util.event_check(e,jbrs_m,jbrs_p)==False: continue
        ch=e.ch
        # --------------------------------------------------------
        if chtmp<ch:# ch is changed
            print "run%d ch%d"%(run,ch)
            chtmp=ch
            # calibration TSpline3
            hname_calib_beamon  = "%s_onsprmon%d_ch%d"%(util.hpht_phc,run,ch)
            hname_calib_beamoff = "%s_offsprmon%d_ch%d"%(util.hpht_phc,run,ch)
            sp3calibon, sp3calibon_zero   = util.GetSpline3CalibDir(fnamecalibon+".root",hname_calib_beamon)
            sp3caliboff, sp3caliboff_zero = util.GetSpline3CalibDir(fnamecaliboff+".root",hname_calib_beamoff)
            if sp3calibon is not None: sp3calibonEval=sp3calibon.Eval
            if sp3calibon_zero is not None: sp3calibon_zeroEval=sp3calibon_zero.Eval
            if sp3caliboff is not None: sp3caliboffEval=sp3caliboff.Eval
            if sp3caliboff_zero is not None: sp3caliboff_zeroEval=sp3caliboff_zero.Eval
        # energy calibration
        ene_phc=0.
        ene_phc_zero=0.
        if e.beam==1: 
            if sp3calibon is not None: ene_phc=sp3calibonEval(e.filt_value_phc)
            if sp3calibon_zero is not None: ene_phc_zero=sp3calibon_zeroEval(e.filt_value_phc)
        elif e.beam==0:
            if sp3caliboff is not None: ene_phc=sp3caliboffEval(e.filt_value_phc)
            if sp3caliboff_zero is not None: ene_phc_zero=sp3caliboff_zeroEval(e.filt_value_phc)
        # --------------------------------------------------------
        # tree fill
        bp_run[0] =            run
        bp_ch[0] =             ch
        bp_ev[0] =             e.ev
        bp_beam[0] =           e.beam
        bp_sprmc[0] =          e.sprmc
        bp_shift1[0] =         e.shift1
        bp_filt_phase[0] =     e.filt_phase
        bp_filt_value[0] =     e.filt_value
        bp_filt_value_dc[0] =  e.filt_value_dc
        bp_filt_value_phc[0] = e.filt_value_phc
        bp_energy[0] =         e.energy
        bp_energy_dc[0] =      e.energy_dc
        bp_energy_phc[0] =     e.energy_phc
        bp_peak_time[0] =      e.peak_time
        bp_peak_value[0] =     e.peak_value
        bp_postpeak_deriv[0] = e.postpeak_deriv
        bp_pretrig_mean[0] =   e.pretrig_mean
        bp_pretrig_rms[0] =    e.pretrig_rms
        bp_pulse_average[0] =  e.pulse_average
        bp_rise_time[0] =      e.rise_time
        bp_timestamp[0] =      e.timestamp
        bp_rowp[0] =           e.rowp
        bp_rown[0] =           e.rown
        bp_pre_region_sum[0] = e.pre_region_sum
        bp_jbr_region_sum[0] = e.jbr_region_sum
        bp_sec_pr_mean[0] =    e.sec_pr_mean
        bp_sec_pr_meanmean[0]= e.sec_pr_meanmean
        bp_sec_pr_maxmean[0] = e.sec_pr_maxmean
        bp_sec_enemean[0] =    e.sec_enemean
        bp_sec_enemax[0] =     e.sec_enemax
        bp_ene[0] =            ene_phc
        bp_ene_zero[0] =       ene_phc_zero
        bp_dt[0] =             e.row_next_extrig_nrp
        
        toutfill()
                        
    if fin.IsOpen(): fin.Close()

fout.cd()
tout.Write()
if fout.IsOpen(): fout.Close("R")
print "%s has been created."%(fnameout+".root")
# --------------------------------------------------------
