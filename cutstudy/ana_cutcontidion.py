#!/usr/bin/env python


""" ana_cutcondition.py is study the cut contition 

History: 
2018-06-20 ; ver 1.0; made from ichinohe-kun's code
"""

__author__ =  'Shinya Yamada (syamada(at)tmu.ac.jp'
__version__=  '1.0'

import commands
import heates_h5sum as h5sum
import matplotlib.pyplot as plt
import numpy as np

h5sum = reload(h5sum)

"""
(1) set parameters 
"""
elo=6900 # Co Kalpha
ehi=6975 # Co Kalpha
clo=-100 # min of histogram
chi=800  # max of histogram
PLOT=True
#h5sum.homedir = "/home/heates/ichinohe"
h5sum.homedir = "/home/heates"

# read hdf5 files
#xraybeam = h5sum.h5sum("xrayon-beamon",[354],debug=False,maxchans=240)
xraybeam = h5sum.h5sum("xrayon-beamon",[366],debug=False,maxchans=240)
#xraybeam = h5sum.h5sum("xrayon-beamon",[353,355,356,359],debug=False,maxchans=240)
#xraybeam = h5sum.h5sum("xrayon-beamon",[269],debug=False,maxchans=240)
beam = h5sum.h5sum("xrayoff-beamon",[136,137,138,139],debug=False,maxchans=240)

#check if hdf5 is strange 
#print xray
#print xraybeam
#print beam

"""
(2) define cut selection 
"""
#sel0=(xray.energy_sum>elo)*(xray.energy_sum<ehi)
sel1=(xraybeam.beam_sum==1)*(xraybeam.energy_sum>elo)*(xraybeam.energy_sum<ehi)
sel2=(beam.beam_sum==1)*(beam.energy_sum>elo)*(beam.energy_sum<ehi)
etag = "%3.2f < E < %3.2f" % (elo,ehi)


"""
(3) plot histogram 
"""

# xray and beam vs. beamonly 
h5sum.plothistSN(xraybeam.sec_pr_mean_sum[sel1],     beam.sec_pr_mean_sum[sel2]     ,"xraybeam","beamonly", 'meanmax 12_' + etag , xmin=clo, xmax = chi)
h5sum.plothistSN(xraybeam.sec_pr_maxmean_sum[sel1],  beam.sec_pr_maxmean_sum[sel2]  ,"xraybeam","beamonly", 'maxmean 12_' + etag , xmin=clo, xmax = chi)
h5sum.plothistSN(xraybeam.sec_pr_meanmean_sum[sel1], beam.sec_pr_meanmean_sum[sel2] ,"xraybeam","beamonly", 'meanmean 12_' + etag , xmin=clo, xmax = chi)

h5sum.plothistSN(xraybeam.sec_enemean_sum[sel1], beam.sec_enemean_sum[sel2] ,"xraybeam","beamonly", 'enemean 12_' + etag , xmin=clo, xmax = chi)
h5sum.plothistSN(xraybeam.sec_enemax_sum[sel1], beam.sec_enemax_sum[sel2] ,"xraybeam","beamonly", 'enemax 12_' + etag , xmin=clo, xmax = chi)
h5sum.plothistSN(xraybeam.sec_fdmean_sum[sel1], beam.sec_fdmean_sum[sel2] ,"xraybeam","beamonly", 'fdmean 12_' + etag , xmin=clo, xmax = chi)
h5sum.plothistSN(xraybeam.sec_fdmax_sum[sel1], beam.sec_fdmax_sum[sel2] ,"xraybeam","beamonly", 'fdmax 12_' + etag , xmin=clo, xmax = chi)


"""
(4) calculate S/N
"""

### plot SN from pr 
fig, (ax1,ax2,ax3) = plt.subplots(3,1, figsize=(8,14))
h5sum.calcSN(xraybeam.sec_pr_mean_sum[sel1],     beam.sec_pr_mean_sum[sel2]     ,"xraybeam","beamonly", 'meanmax ' + etag ,   ax1)
h5sum.calcSN(xraybeam.sec_pr_maxmean_sum[sel1],  beam.sec_pr_maxmean_sum[sel2]  ,"xraybeam","beamonly", 'maxmean ' + etag ,   ax2)
h5sum.calcSN(xraybeam.sec_pr_meanmean_sum[sel1], beam.sec_pr_meanmean_sum[sel2] ,"xraybeam","beamonly", 'meanmean' + etag ,   ax3)
plt.tight_layout()
odir = "figures_calcSN"	
commands.getoutput('mkdir -p ' + odir)
plt.savefig(odir + "/xraybeam_beamonly_calcSN.png")


### plot SN from energy 
fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1, figsize=(8,14) )
h5sum.calcSN(xraybeam.sec_enemean_sum[sel1],     beam.sec_enemean_sum[sel2]  ,"xraybeam","beamonly", 'enemean ' + etag ,   ax1)
h5sum.calcSN(xraybeam.sec_enemax_sum[sel1],      beam.sec_enemax_sum[sel2]   ,"xraybeam","beamonly", 'enemax  ' + etag ,   ax2)
h5sum.calcSN(xraybeam.sec_fdmean_sum[sel1],      beam.sec_fdmean_sum[sel2]   ,"xraybeam","beamonly", 'fdmean  ' + etag ,   ax3)
h5sum.calcSN(xraybeam.sec_fdmax_sum[sel1],       beam.sec_fdmax_sum[sel2]    ,"xraybeam","beamonly", 'fdmax   ' + etag ,   ax4)
plt.tight_layout()
odir = "figures_calcSN"	
commands.getoutput('mkdir -p ' + odir)
plt.savefig(odir + "/xraybeam_beamonly_calcSNene.png")


"""
(5) look at how the cut works 
"""

val=37
cut=(xraybeam.sec_pr_mean_sum < val)
#cut=(xraybeam.beam_sum==1)*(xraybeam.sec_pr_mean_sum < val)
xraybeam.plotspec(cut,"meanmaxcut-less-than" +str(val), xmin = elo, xmax = ehi, beam=True)

val=13
#val=69
cut=(xraybeam.sec_enemean_sum < val)
#cut=(xraybeam.beam_sum==1)*(xraybeam.sec_pr_meanmean_sum < val)
xraybeam.plotspec(cut,"enemeancut-less-than" +str(val), xmin = elo, xmax = ehi, beam=True)

if PLOT:
    plt.show()
