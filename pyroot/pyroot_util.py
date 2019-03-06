# utilities for pyroot
import sys
import os
import pandas as pd
import optparse

datadir="%s"%(os.environ['HEATESDATADIR'])
dumprootdir="%s/dumproot"%(os.environ['HEATESDATADIR'])
rootdir="%s/root/output"%(os.environ['HEATESANADIR'])
pardir="%s/root/par"%(os.environ['HEATESANADIR'])
outdir="./output"# ./=pyroot
figdir="./fig"# ./=pyroot
hpht_phc="hpht_phc"# important hname
RUNINFO="%s/csv/data_TMU_2018U.csv"%(os.environ['HEATESANADIR'])
tree_name="chanall"

parser = optparse.OptionParser()
parser.add_option('--pre',  dest='cut_pre', action="store",type=int, help='set cut for pre samples',default=0)
parser.add_option('--post',  dest='cut_post', action="store",type=int, help='set cut for post samples',default=0)

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
