#!/usr/bin/env python

import h5py
import numpy as np
from sys import argv

import sys
import os
import sys
import math
import commands

import heatesclock as hc

argvs = sys.argv
argc = len(argvs)

if (argc != 3):
    print ('\n[ERROR] Usage: # python %s runn beamn \n' % argvs[0] )
    quit()
    
runn = int(argvs[1])
runs = "%04d" % runn
beamn = int(argvs[2])
beams = "%05d" % beamn
print "***** start run = ", runs, " beam = ", beams

ANADIR=os.environ.get("HEATESANADIR","")
TESDATADIR=os.environ.get("HEATESDATADIR","")
CLOCKDATADIR=os.environ.get("CLOCKDATADIR","")
print "ANADIR  = ", ANADIR
print "TESDATADIR = ", TESDATADIR
print "CLOCKDATADIR = ", CLOCKDATADIR

testrig = hc.TesTrig(DATADIR=TESDATADIR, runs=runs)
#plot if you wish to make a figure  
#testrig.qlplot(0,1000)

btrig = hc.BeamTrig(DATADIR=CLOCKDATADIR, beams=beams)
#plot if you wish to make a figure  
#btrig.qlplot(0,1000)

print "*  Matching check start  *"

CLOCKOUTDIR=TESDATADIR + "clock"
outfname= TESDATADIR + "/clock/" + btrig.name + "_spill.txt"
print outfname
outfile=open(outfname, "w")

offset = 0
spilloffs=3
for i in np.arange(btrig.nspill):
    # dtimeStill is a segment of the nsp-th spill
    # num is the line number of the nsp-th spill
    dtimeSpill, num = btrig.getSpill(nsp = i, width=5, offset=spilloffs)
    # 
    offset, offsetToBeam = testrig.findPattern(dtimeSpill, num, nsp = i, offset = offset)
    # offset : number of the column found in TES clock file
    # offsetToBeam : the line difference bet. Beam Spill file and TES external trigger file
    diff = int(testrig.tesclock[offset] - btrig.beamclock2[num])

    #    outstr = str(i) + ", " + str(int(testrig.tesclock[offset])) + ", " + str(int(btrig.beamclock2[num])) + ", " + str(diff) + ", " + str(offsetToBeam)
    #    outstr = str(i) + ", " + str(int(num)) + ", " + str(diff)
    print testrig.tesclock[offset],btrig.beamclock2[num],diff
    outstr = '{:>10} {:>10} {:>10} {:>15}'.format(i,int(num),offsetToBeam,diff)
    outfile.write(outstr + "\n")

outfile.close()
    
print "*  Matching check end  *"

