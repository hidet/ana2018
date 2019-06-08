# utilities for pyroot
import sys
import os
import pandas as pd
import numpy as np
import mass
import optparse
import ROOT

number_of_rows=30.
frame_timebase=7.18e-06
row_timebase=frame_timebase/number_of_rows
row_offset=(1024-256)*number_of_rows

datadir="%s"%(os.environ['HEATESDATADIR'])
dumprootdir="%s/dumproot"%(os.environ['HEATESDATADIR'])
rootdir="%s/root/output"%(os.environ['HEATESANADIR'])
pardir="%s/root/par"%(os.environ['HEATESANADIR'])
outdir="./output"
figdir="./fig"
RUNINFO="%s/csv/data_TMU_2018U.csv"%(os.environ['HEATESANADIR'])

tree_name="chanall"

hpht_phc       = "hpht_phc";# important hname
hene_phc       = "hene_phc";
fitTES_tag     = "fitTES";# fit function                            
fitElem_tag    = "fitTES";# fit function                            
fitUser_tag    = "fitTES";# fit function                            
ps_tag         = "peak_search";# peak search                        
sp3meander_tag = "sp3meander";# energy calibration                  
sp3mean_tag    = "sp3mean";# energy calibration                     
gmean_tag      = "gmean";# energy calibration                       
gcalib_tag     = "gcalib";# energy calibration                      
sp3calib_tag   = "sp3calib";# energy calibration                    
hene_tag       = "hene";# energy converted histo                    
calib_tag      = "calib";# output root file    

parser = optparse.OptionParser()
parser.add_option('--pre',  dest='cut_pre', action="store",type=int, help='set cut for pre samples',default=0)
parser.add_option('--post',  dest='cut_post', action="store",type=int, help='set cut for post samples',default=0)
parser.add_option('--spill',  dest='spill', action="store",type=int, help='set spill 0:off 1:on',default=-1)


root_type_names = [
    "Bool_t",
    "Char_t",
    "UChar_t",
    "Short_t",
    "UShort_t",
    "Int_t",
    "UInt_t",
    "Long64_t",
    "ULong64_t",
    "Float_t",
    "Double_t"]

python_types = [
    "B",#       unsigned char   int                 1 (used as boolean)
    "b",#       signed char     int                 1
    "B",#       unsigned char   int                 1
    "h",#       signed short    int                 2
    "H",#       unsigned short  int                 2
    "i",#       signed int      int                 2
    "I",#       unsigned int    long                2
    "l",#       signed long     int                 4
    "L",#       unsigned long   long                4
    "f",#       float           float               4
    "d"]#       double          float               8


def type_convert_root_to_python(n):
    if n in root_type_names:
        idx = root_type_names.index(n)
        return python_types[idx]
    else:
        print "Error: ROOT type is strange %s"%n
        return ""

def get_noise_list(runs):
    if os.path.isfile(RUNINFO)==False: 
        print "%s is missing"%RUNINFO
        sys.exit(0)
    df = pd.read_csv(RUNINFO)
    run_list = df.iloc[:,0].tolist()
    noise_list = df.iloc[:,1].tolist()
    inds=[run_list.index(run) for run in runs]
    nois=[int(noise_list[ind]) for ind in inds]
    return nois

def backup_rootfile(fname):
    if os.path.isfile(fname)==False: return False
    end=fname.find(".root")
    fnamebkp=fname[0:end]
    from datetime import datetime
    dt2=datetime.now().strftime('%Y%m%d%H%M%S')
    fnamebkp+="_%s.root"%(dt2)
    from subprocess import call
    call('cp -p %s %s'%(fname,fnamebkp), shell=True)
    print "backup file: %s has been created"%(fnamebkp)
    return True

def check_continue():
    # raw_input for python2
    # input for python3
    # raw_input returns the empty string for "enter"
    yes = {'YES','Yes','Y','yes','y','ye'}
    no = {'NO','No','N','no','n'}
    ans = raw_input("This will take a long time. Do you want to do? (y/n)")
    if ans in yes:
        return True
    elif ans in no:
        sys.exit(0)
    else:
        sys.stdout.write("Please respond with 'yes' or 'no'")

def get_hist_bins(h,elo,ehi):
    bin_s=h.FindBin(elo)
    bin_e=h.FindBin(ehi)
    bins=np.arange(bin_s,bin_e+1)
    nbin = bin_e - bin_s + 1
    bw=h.GetBinWidth(bin_s)
    minx=h.GetBinLowEdge(bin_s)
    maxx=h.GetBinLowEdge(bin_e)+bw
    bin_edges=np.arange(minx,maxx+bw,bw)
    hist = [h.GetBinContent(b) for b in bins]
    hist = np.array(hist)
    return hist,bin_edges

def linefit(hist,bin_edges,linename):
    fitter = mass.getfitter(linename)
    fitter.fit(hist,bin_edges,plot=False)
    params = fitter.last_fit_params[:]
    params[fitter.param_meaning["tail_frac"]]=0.25
    params[fitter.param_meaning["dP_dE"]]=1
    params[fitter.param_meaning["resolution"]]=6
    fitter.fit(hist,bin_edges,params=params,hold=[2],vary_tail=True, plot=False)
    return fitter


def event_check(ev,jbrs_m,jbrs_p):
    if ev.good!=1: return False
    if ev.primary!=1: return False
    if ev.jbr_region_sum>jbrs_m and ev.jbr_region_sum<jbrs_p: return True
    else: return False

def GetSpline3CalibDir(fname, hname):
    '''
    energy calibration with ROOT.TSpline3
    '''
    f= ROOT.TFile.Open(fname,"read")
    if f.IsZombie():
        print "Warning: cannot open file "+fname
        return None, None
    fd = f.GetDirectory(hname);
    if not fd:# fd could be nullptr
        f.Close()
        return None, None
    sp3_name = "%s_%s"%(sp3calib_tag,hname)
    sp3zero_name = "%s_%s_zero"%(sp3calib_tag,hname)
    sp3 = fd.Get(sp3_name)
    sp3zero = fd.Get(sp3zero_name)
    if not sp3 or not sp3zero:# sp3 or sp3zero could be nullptr
        f.Close()
        return None, None
    return sp3, sp3zero


# --------------------------------------------------------
# He3
nruns_he3=96
runs_he3 = [160,161,162,163,165,166,167,170,171,172,
            175,176,177,178,179,180,181,182,183,186,
            187,188,189,190,191,192,193,194,195,196,
            197,198,201,202,203,204,205,206,207,208,
            209,210,212,215,216,217,218,222,223,224,
            225,226,227,228,231,232,233,239,240,241,
            242,259,260,261,263,264,266,267,268,269,
            270,271,272,275,276,277,278,279,281,282,
            283,284,286,287,290,291,292,293,294,295,
            296,297,298,299,300,301]

# He4
nruns_he4=70
runs_he4 = [320,321,327,328,329,330,331,332,333,334,
            335,336,337,338,340,341,345,346,347,348,
            349,350,351,352,353,354,355,356,999,359,
            363,364,366,367,368,369,370,371,372,373,
            374,375,376,377,381,382,383,384,385,386,
            387,389,395,396,397,398,399,400,401,402,
            403,406,407,408,410,411,421,422,423,424]


# for initial guess
# KHeX
# Relative intensity
#           Batty(1976)  Baird(1983)
# KHeX La : 100 +- 4     100 +- 4
#      Lb :  37 +- 5      26 +- 4
#      Lg :  20 +- 4      18 +- 3
#                          4 +- 3
#
# T.Koike calculation
KHE4LA_K = 6463.46;
KHE4LB_K = 8721.73;
KHE4LG_K = 9766.78;
KHE4LD_K = 10334.43;

KHE3LA = 6224.6;

PIC4TO3 = 6428.5;
PI3HE2TO1 = 10678.0;

# from Koike-san
K_PIC4F3D_0 = 6428.78;
K_PIC4D3P_0 = 6435.66;
K_PIC4F3D_1 = 6428.57;
K_PIC4D3P_1 = 6435.39;
K_PIC4F3D_2 = 6428.46;
K_PIC4D3P_2 = 6435.23;

# from HIRENZAKI-san
H_PIC4D3P_STRONG = 0.7784;
H_PIC4F3D = 6428.742;
H_PIC4D3P = 6436.390;
