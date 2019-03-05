#!/usr/bin/env python


""" heates_h5sum.py is a class to sum several values 

History: 
2018-06-20 ; ver 1.0; made from ichinohe-kun's code
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
import sys
from numpy import linalg as LA
import csv
import h5py

params = {'xtick.labelsize': 13, # x ticks
                    'ytick.labelsize': 13, # y ticks
                    'legend.fontsize': 9
                              }
plt.rcParams['font.family'] = 'serif' 
plt.rcParams.update(params)


def calcSN(xx,bb,label1,label2,title,ax1,xmin=-50,xmax=800):
    """
    Args:
        xx (numpy array) : input array1 
        bb (numpy array) : input array2
        label1 : label for xx
        label2 : label for bb
        title (string): title for the plot
        ax1 : axies for the plot
        xmin : minimum of the histogram 
        xmax : maximum of the histogram 
    """

    if len(xx) > 0 and len(bb) > 0:
        pass
    else:
        print "ERROR : either is empty."
        print "xx = ", xx
        print "yy = ", yy
        sys.exit(1)
    
    idxs=np.arange(xmin,xmax,1)
    nx=float(len(xx))
    nb=float(len(bb))
    lxs=np.array([np.sum(xx<i)/nx for i in idxs]) # xx for cumulative plot
    lbs=np.array([np.sum(bb<i)/nb for i in idxs]) # bb for cumulative plot
    sns=lxs/lbs # S/N 
    snsmax=np.nanargmax(sns)+xmin
    sn95x = np.min(np.where(lxs>0.95))
    sn95y = sns[sn95x]    
    x95=sn95x+xmin
    
    ax1.set_title(title)
    ax1.plot(idxs,lxs,color='k',label=label1,ls='-')
    ax1.plot(idxs,lbs,color='k',label=label2,ls='--')
    ax1.axvline(x=x95,color='k',linestyle='--',linewidth=0.5)
    xt=ax1.get_xticks()
    xt=np.append(xt,snsmax)
    xt=np.append(xt,x95)
    xtl=xt.tolist()
    xtl[-2]='%d'%(snsmax)
    xtl[-1]='%d'%(x95)
    ax1.set_xticks(xt)
    ax1.set_xticklabels(xtl)
    ax1.set_ylabel('',color='k')
    ax1.tick_params('y',colors='k')
    ax1.legend()
    ax2 = ax1.twinx()
    ax2.plot(idxs,sns,color='r')
    ax2.axhline(y=sn95y,color='r',linestyle='--',linewidth=0.5)
    ax2.set_ylabel(label1 + '/' + label2,color='r')
    ax2.tick_params('y',colors='r')
    ax2.axvline(x=snsmax,color='r',linestyle='--',linewidth=0.5)
    ax2.set_ylim(0.5,5)


    #print "calcSN[%s,%s,%s] sn95x=%03.3f, sn95y=%03.3f" %(title,label1,label2,sn95x,sn95y)
    print "calcSN[%s,%s,%s] sn95x=%03.3f, sn95y=%03.3f" %(title,label1,label2,x95,sn95y)


def plothistSN(xx,bb,label1,label2,title,xmin=0,xmax=1e5,xnbin=10):
    """
    Args:
        xx (numpy array) : input array1 
        bb (numpy array) : input array2
        label1 : label for xx
        label2 : label for bb
        title (string): title for the plot
        ax1 : axies for the plot
        xmin : minimum of the histogram 
        xmax : maximum of the histogram 
    """

    if len(xx) > 0 and len(bb) > 0:
        pass
    else:
        print "ERROR : either is empty."
        print "xx = ", xx
        print "yy = ", yy
        sys.exit(1)
    

    F = plt.figure(figsize=(12,8))
    
    hxx,_,_  = plt.hist(xx,bins=np.arange(xmin,xmax,xnbin),histtype='step',color='red',label=label1)
    hbb,_,_  = plt.hist(bb,bins=np.arange(xmin,xmax,xnbin),histtype='step',color='black',label=label2)    
    plt.title("plot summed parameters")
    plt.yscale('log')
    plt.xlabel(title)
    plt.grid(True)
    plt.legend()
    
    odir = "figures_histSN"
    commands.getoutput('mkdir -p ' + odir)
    plt.savefig(odir + "/" + title.replace(" ","-").replace(".","p") + ".png")
    plt.close()
    print "plothist[%s,%s,%s] counts of xx =%06d, bb = %06d" %(title,label1,label2,len(xx),len(bb))
    
    
    
class h5sum():

    home = "/home/heates/ichinohe"
    #homedir = "/Users/syamada/work/ana/HEATES/JPARC201806_E62/ichinohe"
    
    def __init__ (self, name, runlist, debug = False, maxchans = 240):    
        self.home = "/home/heates/ichinohe"
#        self.home = "/home/heates"        
        self.name = name
        self.runlist = runlist
        self.fullfilename = []
        self.shortfilename = []
        self.sname = []
        
        self.beam_sum=[]
        self.energy_sum=[]        
        self.good_sum=[]        
        self.grouptrig_sum=[]        
        self.sec_pr_max_sum=[]
        self.sec_pr_mean_sum=[]
        self.sec_pr_maxmean_sum=[]
        self.sec_pr_meanmean_sum=[]
        self.sec_pr_maxmin_sum=[]
        self.sec_pr_meanmin_sum=[]
        self.sec_pr_sum_sum=[]

        self.sec_enemean_sum=[]
        self.sec_enemax_sum=[]
        self.sec_fdmean_sum=[]
        self.sec_fdmax_sum=[]
        
        self.debug = debug
        self.name = name
        # copied from ana_galen_single.py, 2018/06/18
        self.BADCHS = [3,9,39,77,83,85,111,337,367,375,423]# initially disconnected
        self.BADCHS.extend([117,203,233])# bad channels
        self.BADCHS.extend([5,177])# strange channels
        self.maxchans = maxchans
        
        for i, runn in enumerate(runlist):

            runs = "%04d" % runn
            self.fullfilename.append(homedir + "/data/TMU_2018U/run%s/run%s_mass.hdf5" % (runs, runs))
            self.shortfilename.append(os.path.basename(self.fullfilename[-1]))
            self.sname.append(self.shortfilename[-1].replace(".hdf5","").replace(".dat",""))
            #extract header information
            if os.path.isfile(self.fullfilename[-1]):
                h5f = h5py.File(self.fullfilename[-1],'r')

                for j, chan in enumerate(h5f):

                    if j > self.maxchans: break # shot cut 
                    
                    if chan in self.BADCHS:
                        print "skip chan = ", chan, " due to bad chan"
                        continue
                    else:
                        cc=h5f[chan]
                        energy=cc['energy'].value
                        try:
                            good=cc['good'].value
                        except:
                            good=None
                        try:
                            beam=cc['beam'].value
                        except:
                            beam=None
                        grouptrig=cc['grouptrig'].value
                        sec_pr_max=cc['sec_pr_max'].value
                        sec_pr_mean=cc['sec_pr_mean'].value
                        sec_pr_maxmean=cc['sec_pr_maxmean'].value
                        sec_pr_meanmean=cc['sec_pr_meanmean'].value
                        sec_pr_maxmin=cc['sec_pr_maxmin'].value
                        sec_pr_meanmin=cc['sec_pr_meanmin'].value
                        sec_pr_sum=cc['sec_pr_sum'].value

                        sec_enemean=cc['sec_enemean'].value                        
                        sec_enemax=cc['sec_enemax'].value                        
                        sec_fdmean=cc['sec_fdmean'].value                        
                        sec_fdmax=cc['sec_fdmax'].value                        

                    try:
                        self.beam_sum.extend(list(beam[grouptrig==-1]))
                    except:
                        self.beam_sum.extend([])

                    try:
                        self.good_sum.extend(list(good[grouptrig==-1]))
                    except:
                        self.good_sum.extend([])

                    self.energy_sum.extend(list(energy[grouptrig==-1]))
                    self.grouptrig_sum.extend(list(grouptrig[grouptrig==-1]))
                    self.sec_pr_max_sum.extend(list(sec_pr_max[grouptrig==-1]))
                    self.sec_pr_mean_sum.extend(list(sec_pr_mean[grouptrig==-1]))
                    self.sec_pr_maxmean_sum.extend(list(sec_pr_maxmean[grouptrig==-1]))
                    self.sec_pr_meanmean_sum.extend(list(sec_pr_meanmean[grouptrig==-1]))
                    self.sec_pr_maxmin_sum.extend(list(sec_pr_maxmin[grouptrig==-1]))
                    self.sec_pr_meanmin_sum.extend(list(sec_pr_meanmin[grouptrig==-1]))
                    self.sec_pr_sum_sum.extend(list(sec_pr_sum[grouptrig==-1]))
                    
                    self.sec_enemean_sum.extend(list(sec_enemean[grouptrig==-1]))                    
                    self.sec_enemax_sum.extend(list(sec_enemax[grouptrig==-1]))                    
                    self.sec_fdmean_sum.extend(list(sec_fdmean[grouptrig==-1]))                    
                    self.sec_fdmax_sum.extend(list(sec_fdmax[grouptrig==-1]))                    

            else:
                print "ERROR. could not open ", self.fullname
                print "quit"
                sys.exit()

        self.beam_sum=np.array(self.beam_sum)
        self.energy_sum=np.array(self.energy_sum)
        self.good_sum=np.array(self.good_sum)
        self.grouptrig_sum=np.array(self.grouptrig_sum)
        self.sec_pr_max_sum=np.array(self.sec_pr_max_sum)
        self.sec_pr_mean_sum=np.array(self.sec_pr_mean_sum)
        self.sec_pr_maxmean_sum=np.array(self.sec_pr_maxmean_sum)
        self.sec_pr_meanmean_sum=np.array(self.sec_pr_meanmean_sum)
        self.sec_pr_maxmin_sum=np.array(self.sec_pr_maxmin_sum)
        self.sec_pr_meanmin_sum=np.array(self.sec_pr_meanmin_sum)
        self.sec_pr_sum_sum=np.array(self.sec_pr_sum_sum)

        self.sec_enemean_sum=np.array(self.sec_enemean_sum)        
        self.sec_enemax_sum=np.array(self.sec_enemax_sum)        
        self.sec_fdmean_sum=np.array(self.sec_fdmean_sum)        
        self.sec_fdmax_sum=np.array(self.sec_fdmax_sum)        
            
       
    def __repr__(self):
        return "%s: %s" %(self.name, self.__dict__)

    
    def plotspec(self, cutarray = None, cname="", beam=False, xmin=0, xmax=1e5, xnbin=10):

        F = plt.figure(figsize=(12,8))
        
        plt.title(self.name)
        ax = plt.subplot(1,1,1)

        if beam:
            beamon=(self.beam_sum==1)
            esum = self.energy_sum[beamon]
        else:
            esum = self.energy_sum
            
        hb,_,_  = plt.hist(esum,bins=np.arange(4000,9000,10),histtype='step',color='black',label="all")
        c67all = np.where( (esum > xmin) & (esum < xmax) )[0]
        
        if len(cutarray) > 0:
            
            if beam:
                ene1=self.energy_sum[ cutarray * beamon]
                ene2=self.energy_sum[~cutarray * beamon]
            else:
                ene1=self.energy_sum[cutarray]
                ene2=self.energy_sum[~cutarray]
                
            c67a1 = np.where( (ene1 > xmin) & (ene1 < xmax) )[0]
            c67a2 = np.where( (ene2 > xmin) & (ene2 < xmax) )[0]
            
            hbm1,_,_= plt.hist(ene1,bins=np.arange(4000,9000,10),histtype='step',color='red',label= cname +  '=TRUE')
            hbnm1,_,_= plt.hist(ene2,bins=np.arange(4000,9000,10),histtype='step',color='pink',label= cname + '=FALSE')

            call = len(c67all)
            c1 = len(c67a1)
            c2 = len(c67a2)
            p1 = 100. * c1/call 
            p2 = 100. * c2/call
            status = "total %7d counts in %d - %d keV, %7d (%02.3f%%) in %s, %7d (%02.3f%%) in %s" % (call, xmin, xmax, c1, p1, cname, c2, p2, "not" + cname)
            print status
            plt.figtext(0.1,0.92,"[" + status + "]")
            
        plt.yscale('log')
        plt.xlabel('energy (eV)')
        plt.grid(True)
        plt.legend()
        
        odir = "figures_h5sum"
        commands.getoutput('mkdir -p ' + odir)
        plt.savefig(odir + "/" + self.name + "_" + cname + ".png")
        plt.close()
