import numpy as np
import h5py
import matplotlib.pylab as plt
import csv
import pandas as pd
import matplotlib.colors as colors
import matplotlib.cm as cm

ANAHOME="/home/heates/"
ANADIR="%s/noda/"%(ANAHOME)# maybe current dir
DATADIR="/a/heates/data/anadata/TMU_2018U/"
RUNINFO="noiseonly.csv"
#RUNINFO="%s/ana/ana_2018_jparc/csv/data_TMU_2018U.csv"%(ANAHOME)

#a = pd.read_csv(RUNINFO)
#b = a.run_noise
#b = a.run_noise[10:]

a = np.loadtxt(RUNINFO, delimiter = ',', dtype = str)
a1 = a[:,0]
a2 = a1[:]
#a2 = a1[:]
print a2
b = a2.astype(np.int32)
#b = a.run_noise[-10:]
#print b

print b

xtime = np.arange(513) * 240e-9 * 30 # 240ns x 30 column
dt = 240e-9*30 # timing resolution of one pixel, 240ns x 30 columns 
fmax = 2 * 1./dt # nyquist frequency
reclen = 513
freq = np.linspace(0,fmax,reclen)


def getnoise(ch, debug = False):

    runid = []
    psd = []
    
    for i, row in enumerate(b):

        run = "{0:04d}".format(row)
        fname = "%s/run%s/run%s_mass.hdf5"%(DATADIR, run, run)

        if debug:
            print i, fname

#        if debug == False:
        try:
            f = h5py.File(fname)        
            noise = f["chan"+str(ch)]["noise_psd"][:]
            runid.append(row)
            psd.append(np.sum( noise[ np.where( freq < 3000)]))
 
#	else:       
        except:
            pass


    runid = np.array(runid)    
    psd = np.array(psd)

    zerocut = np.where(psd > 0)
    runid = runid[zerocut]
    psd = psd[zerocut]

    return runid, psd


F = plt.figure(figsize=(12,8))

ax = plt.subplot(1,1,1)
plt.xlabel("Frequency (Hz)")
#plt.xlabel("Time (s)")
plt.ylabel("integrated noise (f<3kHz)")
plt.yscale("linear")
plt.xscale("linear")
plt.xlabel("run number")
plt.grid(True)

for ch in [1,97,114,313]:
    runid, psd = getnoise(ch)
    plt.errorbar(runid, psd, fmt="o-", label = "ch="+str(ch),alpha = 0.9)
    plt.legend()
    
plt.show()    













