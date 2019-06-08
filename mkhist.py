import khe_hdf5 as h
import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys
# -------- for ROOT file --------------
try:
    import ROOT
except ImportError:
    raise ValueError('ERROR: cannot import pyROOT')
ROOT.gROOT.SetBatch(1)

cut_on=[['good',1,None],
        ['beam',1,None],
        ['grouptrig',-1,None]
        ]
cut_off=[['good',1,None],
         ['beam',0,None],
         ['grouptrig',-1,None]
         ]
cut_3he=list(cut_on)
cut_3he.append(['energy',6200,6240])
cut_4he=list(cut_on)
cut_4he.append(['energy',6445,6485])
cut_timing=list(cut_on)
cut_timing.append(['rows_until_next_external_trigger',75,80])
cut_asyn1=list(cut_on)
cut_asyn1.append(['rows_until_next_external_trigger',45,65])
cut_asyn2=list(cut_on)
cut_asyn2.append(['rows_until_next_external_trigger',90,110])
cut_group=list(cut_timing)
cut_group.append(['sec_pr_mean',-999,37])
cut_de=list(cut_timing)
cut_de.append(['rows_until_next_external_trigger_clock3',70,90])
cut_asynde=list(cut_asyn1)
cut_asynde.append(['rows_until_next_external_trigger_clock3',30,80])
cut_asynde.append(['sec_pr_mean',-999,37])

cutlist=[cut_off,cut_on,cut_3he,cut_4he,cut_timing,cut_group,cut_de,cut_asyn1,cut_asyn2,cut_asynde]
cutname=['_off','_beam','_he3','_he4','_timing','_group','_de','_asyn1','_asyn2','_asynde']
chs=range(1,480,2)

def mkhist(infile,outfile,PRINT=False):
    hdf=h5py.File(infile)    
    fopt='recreate'
    f=ROOT.TFile(outfile,fopt)
    for cutarray,add in zip(cutlist,cutname):
        if PRINT: print add,cutarray
        cut=h.get_cut_list(hdf,chs,cutarray)
        hname='energy_spectrum'+add
        bins=[10000,0,10000]
        hist=h.mkroot_1d(hname,bins,hdf,chs,'energy',cut)
        hname='timing'+add
        bins=[1000,0,200]
        hist=h.mkroot_1d(hname,bins,hdf,chs,'rows_until_next_external_trigger',cut)
        hname='timing_with_de'+add
        bins=[1000,0,200]
        hist=h.mkroot_1d(hname,bins,hdf,chs,'rows_until_next_external_trigger_clock3',cut)
        hname='sec_pr_mean'+add
        bins=[10000,-200,800]
        hist=h.mkroot_1d(hname,bins,hdf,chs,'sec_pr_mean',cut)

    f.Close()

if __name__ == '__main__':
    infile = sys.argv[1]
    run = sys.argv[2]
    outfile = 'hist_run' + run + '.root'
    mkhist(infile,outfile)
    
