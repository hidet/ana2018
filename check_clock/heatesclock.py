#!/usr/bin/env python


""" clock.py is to define two classes 

History: 
2018-06-16 ; ver 1.0, from scratch
"""

__author__ =  'Shinya Yamada (syamada(at)tmu.ac.jp'
__version__=  '1.0'

import os
import sys
import math
import commands
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, FixedLocator
import matplotlib.colors as cr
import optparse
from matplotlib.font_manager import fontManager, FontProperties
import matplotlib.pylab as plab
import re
from numpy import linalg as LA
import csv
import h5py


params = {'xtick.labelsize': 13, # x ticks
          'ytick.labelsize': 13, # y ticks
          'legend.fontsize': 9          
                    }
plt.rcParams['font.family'] = 'serif' 
plt.rcParams.update(params)


class TesTrig():
    
    def __init__ (self, DATADIR, runs, debug = False):
        self.fullfilename = "%s/run%s/run%s_extern_trig.hdf5" % (DATADIR, runs, runs)
        self.shortfilename = os.path.basename(self.fullfilename)
        self.name = self.shortfilename.replace(".hdf5","").replace(".dat","")
        self.debug = debug
        self.dt = 240. * 1e-9
        
        #extract header information
        if os.path.isfile(self.fullfilename):
            self.h5file = h5py.File(self.fullfilename,'r')
            self.tesclock = np.array(self.h5file['trig_times'])
            self.dtesclock = np.diff(self.tesclock)
            self.dtesclock = self.dtesclock.astype('int64')
            self.time = (self.tesclock - self.tesclock[0]) * self.dt # time(sec)
            
            self.exposure = self.time[-1] - self.time[0]  # (sec)
            self.dtime = self.dtesclock * self.dt  # time(sec)            
            self.detectSpill()
            
        else:
            print "ERROR : %s could not open " % self.fullfilename
            sys.exit()

    def detectSpill(self):
        
        self.spillid = np.where(self.dtime > 2.)[0]
        self.nspill = len(self.spillid)
        print "..... ", self.shortfilename, " # of spill = ", self.nspill, " duration = ", self.exposure
        
        
    def findPattern(self, inputpattern, num, offset = 0, debug = False, nsp = 1):

        """
        inputpattern : input clock patten used for matching 
        num : the line number of the nsp-th spill
        offset : number of offset to look for the pattern from num
        nsp : number of Spill 

        return         
        j : number of the column found in TES clock file
        offsetToBeam : the line difference bet. Beam Spill file and TES external trigger file 
        """

        fixedoffset = 10
        
        iny = inputpattern
        samplelen = len(iny)

        istart = offset - fixedoffset if offset >  fixedoffset else 0
        
        for j in np.arange(istart, len(self.dtesclock)-samplelen,1):
#        for i in np.arange(num - noffset, len(self.dtesclock),1):            
            inx = self.dtesclock[j : j + samplelen]
            flag = np.allclose(iny,inx, rtol=1e-05, atol=2.1) # since two clock could be different.
            
#            print i, inx, iny, flag

            if flag:
                offsetToBeam = j - num
                print "Spill # %06d " % (nsp) , " offset-to-Beam = ", offsetToBeam , " (", j, " - ", num, ")"
                return j, offsetToBeam

        if flag == False:
            print "ERROR, there is no corresponding spill."
            print nsp, num, offset, inputpattern
            return -1, -1

        
    def plotpattern(self, inputpattern, x1 = None, x2 = None):
         """
         plot just to grasp the content
         """
         x = self.time
         dx = self.dtime         
         y  = self.tesclock
         dy = self.dtesclock
         iny = inputpattern

         tag = "plotpattern"
         
         if x1 is not None:
             x = x[x1:]
             dx = dx[x1:]             
             y = y[x1:]
             dy = dy[x1:]
             iny = iny[x1:]
             tag += "_xmin" + str(x1)
            
         if x2 is not None:
             x = x[:x2]
             dx = dx[:x2]                          
             y = y[:x2]
             dy = dy[:x2]
             iny = iny[:x2]             
             tag += "_xmax" + str(x2)
             
         odir = "figures_plotpattern"
         commands.getoutput('mkdir -p ' + odir)

         F = plt.figure(figsize=(12,8))

         plt.title(self.shortfilename)
         ax = plt.subplot(1,1,1)

         tmpx = np.arange(len(dy))
         plt.errorbar(tmpx, dy, fmt = "r.", label="TES clock differential", alpha = 0.3)
         tmpx = np.arange(len(iny))
         plt.errorbar(tmpx, iny, fmt="b.-", label="input patten", alpha = 0.3)
         plt.ylim(0,1e4)
         plt.legend(numpoints=1, frameon=False, loc='best')

         plt.savefig(odir + "/" + self.name + "_" + tag + ".png")
            
            
    def qlplot(self, x1 = None, x2 = None):
         """
         plot just to grasp the content
         """
         x = self.time
         dx = self.dtime         
         y  = self.tesclock
         dy = self.dtesclock

         tag = "xrange"
         
         if x1 is not None:
             x = x[x1:]
             dx = dx[x1:]             
             y = y[x1:]
             dy = dy[x1:]
             tag += "_xmin" + str(x1)
            
         if x2 is not None:
             x = x[:x2]
             dx = dx[:x2]                          
             y = y[:x2]
             dy = dy[:x2]
             tag += "_xmax" + str(x2)
             
         odir = "figures_qlplot"
         commands.getoutput('mkdir -p ' + odir)

         F = plt.figure(figsize=(12,8))

         plt.title(self.shortfilename)
         ax = plt.subplot(4,1,1)
                          
         plt.plot(y, label="TES clock raw value")             
         plt.legend(numpoints=1, frameon=False, loc='best')
         
         ax = plt.subplot(4,1,2)
         plt.plot(dy, label="TES clock differential")
         plt.legend(numpoints=1, frameon=False, loc='best')

         ax = plt.subplot(4,1,3)
         plt.plot(x, label="time (sec)")
         plt.legend(numpoints=1, frameon=False, loc='best')

         ax = plt.subplot(4,1,4)
         plt.plot(dx, label="diff time (sec)")
         plt.legend(numpoints=1, frameon=False, loc='best')

         plt.savefig(odir + "/" + self.name + "_" + tag + ".png")
         
         

class BeamTrig():

    def __init__ (self,DATADIR, beams, debug = False):

        self.fullfilename = '%s/clock_run%s.dat' % (DATADIR, beams)
        self.shortfilename = os.path.basename(self.fullfilename)
        self.name = self.shortfilename.replace(".hdf5","").replace(".dat","")    
        self.debug = debug
        self.dt = 240. * 1e-9
        
        #extract header information
        if os.path.isfile(self.fullfilename):
            self.beamclock = np.loadtxt(self.fullfilename, usecols=(0))
            self.beamclock2 = self.beamclock + self.beamclock
            self.dbeamclock2 = np.diff(self.beamclock2)
            self.dbeamclock2 = self.dbeamclock2.astype('int64')
            self.time = (self.beamclock2 - self.beamclock2[0]) * self.dt # time(sec)
            self.dtime = self.dbeamclock2  * self.dt # time(sec)
            self.exposure = np.sum( self.dtime[ np.where( self.dtime > 0 ) ] ) # (sec)            
            self.detectSpill()
            
        else:
            print "ERROR : %s could not open " % self.fullfilename
            sys.exit()

            
    def detectSpill(self):
        self.spillid = np.insert(np.where(self.dtime < -1.)[0],0,0)
        self.nspill = len(self.spillid)
        print "..... ", self.shortfilename, " # of spill = ", self.nspill, " duration(beam ON) = ", self.exposure

        
    def getSpill(self, nsp = 3, width = 10, offset = 3, debug = False):
        """
        *input*
        nsp : number of the spill 
        width :  width used for patten matching
        offset : offset for pattern matching
        
        *output*
        self.dtimeSpill : a segment of the wave packet 
        self.num : the line number of file when the nsp-th spill 
        """
        
        self.num = self.spillid[nsp]
        self.dtimeSpill = self.dbeamclock2[ self.num + offset : self.num + offset + width ]        

        if debug:
            odir = "figures_getSpill"
            commands.getoutput('mkdir -p ' + odir)

            F = plt.figure(figsize=(12,8))        
            plt.title("safe pattern of " + str(nsp) + " in " +   self.shortfilename)
            ax = plt.subplot(1,1,1)                          
            plt.plot(self.dtimeSpill, label="safe still pattern in dbeamclock2")             
            plt.legend(numpoints=1, frameon=False, loc='best')    
            plt.savefig(odir + "/difftime_" + self.name + "_" + tag + ".png")
            
        return self.dtimeSpill, self.num+offset
    
                    
    def qlplot(self, x1 = None, x2 = None):
         """
         plot just to grasp the content
         """
         y  = self.beamclock2
         dy = self.dbeamclock2
         x = self.time
         dx = self.dtime          
         
         tag = "xrange"
         
         if x1 is not None:
             x = x[x1:]
             dx = dx[x1:]                                       
             y = y[x1:]
             dy = dy[x1:]
             tag += "_xmin" + str(x1)
             
         if x2 is not None:
             x = x[:x2]
             dx = dx[:x2]                                       
             y = y[:x2]
             dy = dy[:x2]
             tag += "_xmax" + str(x2)
             
         odir = "figures_qlplot"
         commands.getoutput('mkdir -p ' + odir)

         F = plt.figure(figsize=(12,8))

         plt.title(self.shortfilename)
         ax = plt.subplot(4,1,1)
                          
         plt.plot(y, label="BEAM clock raw value x 2 ")             
         plt.legend(numpoints=1, frameon=False, loc='best')
         
         ax = plt.subplot(4,1,2)
         plt.plot(dy, label="BEAM clock differential")
         plt.legend(numpoints=1, frameon=False, loc='best')

         ax = plt.subplot(4,1,3)
         plt.plot(x, label="time (sec)")
         plt.legend(numpoints=1, frameon=False, loc='best')

         ax = plt.subplot(4,1,4)
         plt.plot(dx, label="diff time (sec)")
         plt.legend(numpoints=1, frameon=False, loc='best')
         
         plt.savefig(odir + "/" + self.name + "_" + tag + ".png")
         
         
         
