from ROOT import gROOT
gROOT.SetBatch()

import ROOT
import mass

import h5py as h5
import pylab as plt
import numpy as np
import pandas as pd
import sys
import os
import optparse

import pyroot_util as util
util=reload(util)


def get_np_from_h5(f,ch,attr):
    return np.array(f[ch][attr])

def get_npbool_from_h5(f,ch,attr):
    return np.array(f[ch][attr],dtype=bool)

def get_adjusted_timestamp(f,ch):
    filt_phase = get_np_from_h5(f,ch,'filt_phase')
    shift1 = get_np_from_h5(f,ch,'shift1')
    row_adjust = -1.*(filt_phase+shift1)
    rowcount = get_np_from_h5(f,ch,'rowcount')
    return util.row_timebase*(rowcount+(row_adjust*util.number_of_rows)-util.row_offset)


run='run0397'
ch='chan45'
#ch='chan145'
title="%s %s"%(run,ch)

#no record length cut
fnames = [
    "%s_noi0393_mass_2018_sprmcon_jbrscon_prime.hdf5"%(run),
    "%s_noi0393_mass_2018_spillon_sprmcon_jbrscon_prime.hdf5"%(run),
    "%s_noi0393_mass_2018_spilloff_sprmcon_jbrscon_prime.hdf5"%(run)
]

# almost optimum rl cut
#fnames = [
#    "run0397_noi0393_mass_2018_pre150_post550_sprmcon_jbrscon_prime.hdf5",
#    "run0397_noi0393_mass_2018_pre150_post550_spillon_sprmcon_jbrscon_prime.hdf5",
#    "run0397_noi0393_mass_2018_pre150_post550_spilloff_sprmcon_jbrscon_prime.hdf5"
#]

fs = [h5.File("%s/%s/%s"%(util.datadir,run,fname),"r") for fname in fnames]
print fs[0][ch].keys()

# for waveform drawing (ms)
x=(np.arange(1024)-256)*util.frame_timebase*1e3
tss = [get_adjusted_timestamp(f,ch) for f in fs]
attr='good';  gs = [get_npbool_from_h5(f,ch,attr) for f in fs]
attr='beam';  beams = [get_npbool_from_h5(f,ch,attr) for f in fs]
attr='prime'; primes = [get_npbool_from_h5(f,ch,attr) for f in fs]
attr='sprmc'; sprmcs = [get_npbool_from_h5(f,ch,attr) for f in fs]
attr='jbrsc'; jbrscs = [get_npbool_from_h5(f,ch,attr) for f in fs]

# good = good && primary && sprmc && jbrsc
gs = [np.logical_and(g,prime) for g,prime in zip(gs,primes)]
gs = [np.logical_and(g,sprmc) for g,sprmc in zip(gs,sprmcs)]
gs = [np.logical_and(g,jbrsc) for g,jbrsc in zip(gs,jbrscs)]
gbons  = [np.logical_and(g,beam) for g,beam in zip(gs,beams)]
gboffs = [np.logical_and(g,~beam) for g,beam in zip(gs,beams)]
# can share the timestamp, good and beam on/off
ts    = tss[0]
ts    = ts-ts[0]
g     = gs[0]
gbon  = gbons[0]
gboff = gboffs[0]

attr='energy';          enes    = [get_np_from_h5(f,ch,attr) for f in fs]
attr='average_pulse';   apls    = [get_np_from_h5(f,ch,attr) for f in fs]
attr='pretrig_mean';    ptms    = [get_np_from_h5(f,ch,attr) for f in fs]
ptm = ptms[0]
attr='filt_value';      fvs     = [get_np_from_h5(f,ch,attr) for f in fs]
attr='filt_value_dc';   fvdcs   = [get_np_from_h5(f,ch,attr) for f in fs]
attr='filt_value_phc';  fvphcs  = [get_np_from_h5(f,ch,attr) for f in fs]


plt.close('all')
plt.ion()

plt.figure()
legends=["all","spill on","spill off"]
for apl,legd in zip(apls,legends):
    plt.plot(x,apl,"-",label=legd)
plt.legend()
plt.xlabel('ms')
plt.ylabel('feed back')
plt.title(title)

plt.figure()
plt.plot(x,apls[2]-apls[0],"-",label="spill off - all")
plt.plot(x,apls[1]-apls[0],"-",label="spill on - all")
plt.legend()
plt.xlabel('ms')
plt.ylabel('feed back')
plt.title(title)

plt.figure()
plt.plot(ts[g],fvs[2][g]-fvs[0][g],".",label="spill off - all")
plt.plot(ts[g],fvs[1][g]-fvs[0][g],".",label="spill on - all")
plt.legend()
plt.xlabel('time after run start (sec)')
plt.ylabel('filt value difference')
plt.title(title)

plt.figure()
plt.plot(ts[g],fvdcs[2][g]-fvdcs[0][g],".",label="spill off - all")
plt.plot(ts[g],fvdcs[1][g]-fvdcs[0][g],".",label="spill on - all")
plt.legend()
plt.xlabel('time after run start (sec)')
plt.ylabel('difference of filt value with drift correction')
plt.title(title)

plt.figure()
plt.plot(ts[g],enes[2][g]-enes[0][g],".",label="spill off - all")
plt.plot(ts[g],enes[1][g]-enes[0][g],".",label="spill on - all")
plt.legend()
plt.xlabel('time after run start (sec)')
plt.ylabel('energy difference [eV]')
plt.title(title)


plt.figure()
plt.plot(ts[g],ptm[g],".")
plt.xlabel('time after run start (sec)')
plt.ylabel('pretrigger mean')
plt.title(title)


plt.figure()
plt.plot(fvs[0][gboff],ptm[gboff],".",label="all spill off")
plt.plot(fvs[0][gbon],ptm[gbon], ".",label="all spill on")
plt.xlabel('filt value')
plt.ylabel('pretrigger mean')
plt.legend()
plt.title(title)

plt.figure()
plt.plot(fvdcs[0][gboff],ptm[gboff],".",label="all spill off")
plt.plot(fvdcs[0][gbon],ptm[gbon], ".",label="all spill on")
plt.xlabel('filt value with drift correction')
plt.ylabel('pretrigger mean')
plt.legend()
plt.title(title)


plt.figure()
plt.plot(ts[gboff],fvs[0][gboff],".",label="spill off")
plt.plot(ts[gbon],fvs[0][gbon],".",label="spill on")
plt.legend()
plt.xlabel('time after run start (sec)')
plt.ylabel('filt value')
plt.title(title)


plt.figure()
plt.plot(ts[gboff],fvdcs[0][gboff],".",label="spill off")
plt.plot(ts[gbon],fvdcs[0][gbon],".",label="spill on")
plt.legend()
plt.xlabel('time after run start (sec)')
plt.ylabel('filt value with drift correction')
plt.title(title)






nbins=1000
hist_range=(12000.,13000.)
plt.figure()
for fv in fvs:
    hist, bin_edges = np.histogram(fv[g],nbins,hist_range)
    plt.hist(fv[g],bin_edges,hist_range,histtype='step')
plt.title(title)


nbins=200
hist_range=(6800.,7000.)
plt.figure()
for ene in enes:
    hist, bin_edges = np.histogram(ene[g],nbins,hist_range)
    plt.hist(ene[g],bin_edges,hist_range,histtype='step')
plt.title(title)



linename="CoKAlpha"
fitter = mass.getfitter(linename)
nbins=100
hist_range = (mass.STANDARD_FEATURES[linename]-50,mass.STANDARD_FEATURES[linename]+50)
vtail=False


plt.figure()
ax1=plt.subplot(311)
hist, bin_edges = np.histogram(enes[0][g],bins=nbins,range=hist_range)
fitter.fit(hist,bin_edges[0:-1],axis=ax1,color='red',vary_bg=True,vary_bg_slope=False,vary_tail=vtail)
plt.title(title)
ax2=plt.subplot(312)
hist, bin_edges = np.histogram(enes[0][gbon],bins=nbins,range=hist_range)
fitter.fit(hist,bin_edges[0:-1],axis=ax2,color='blue', vary_bg=True,vary_bg_slope=False,vary_tail=vtail)
plt.title("template_all: beam-on events")
ax3=plt.subplot(313)
hist, bin_edges = np.histogram(enes[0][gboff],bins=nbins,range=hist_range)
fitter.fit(hist,bin_edges[0:-1],axis=ax3,color='orange',vary_bg=True,vary_bg_slope=False,vary_tail=vtail)
plt.title("template_all: beam-off events")
plt.tight_layout()


plt.figure()
ax1=plt.subplot(311)
hist, bin_edges = np.histogram(enes[1][g],bins=nbins,range=hist_range)
fitter.fit(hist,bin_edges[0:-1],axis=ax1,color='red',vary_bg=True,vary_bg_slope=False,vary_tail=vtail)
plt.title(title)
ax2=plt.subplot(312)
hist, bin_edges = np.histogram(enes[1][gbon],bins=nbins,range=hist_range)
fitter.fit(hist,bin_edges[0:-1],axis=ax2,color='blue', vary_bg=True,vary_bg_slope=False,vary_tail=vtail)
plt.title("template_beamon: beam-on events")
ax3=plt.subplot(313)
hist, bin_edges = np.histogram(enes[1][gboff],bins=nbins,range=hist_range)
fitter.fit(hist,bin_edges[0:-1],axis=ax3,color='orange',vary_bg=True,vary_bg_slope=False,vary_tail=vtail)
plt.title("template_beamon: beam-off events")
plt.tight_layout()


plt.figure()
ax1=plt.subplot(311)
hist, bin_edges = np.histogram(enes[2][g],bins=nbins,range=hist_range)
fitter.fit(hist,bin_edges[0:-1],axis=ax1,color='red',vary_bg=True,vary_bg_slope=False,vary_tail=vtail)
plt.title(title)
ax2=plt.subplot(312)
hist, bin_edges = np.histogram(enes[2][gbon],bins=nbins,range=hist_range)
fitter.fit(hist,bin_edges[0:-1],axis=ax2,color='blue', vary_bg=True,vary_bg_slope=False,vary_tail=vtail)
plt.title("template_beamoff: beam-on events")
ax3=plt.subplot(313)
hist, bin_edges = np.histogram(enes[2][gboff],bins=nbins,range=hist_range)
fitter.fit(hist,bin_edges[0:-1],axis=ax3,color='orange',vary_bg=True,vary_bg_slope=False,vary_tail=vtail)
plt.title("template_beamoff: beam-off events")
plt.tight_layout()





for f in fs:
    if (f): f.close()


    
