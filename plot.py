import khe_hdf5 as h
import csv
import sys

f = open('He4.runs','r')
runs=map(int,f.readline().split())
hdf=h.khe_hdf5(runs)
#runs2=range(345,358)
#chs=[1,7,11,13]
runs2=[421]
chs=[13]
hdf.plot_average_pulse(runs2,chs,'average_pulse',False)
sys.exit()
hdf.plot_calib_trend(runs,chs,peaks=['CrKAlpha','CoKAlpha','CuKAlpha'])
attrs=['pretrig_mean','filt_value','filt_value_dc','energy']
cutarray=[['good',1,None],
          ['beam',0,None],
          ['grouptrig',-1,None],
          ['energy',6700,7100]]
hdf.plot_trend_multi(runs2,chs,attrs,cutarray)
hdf.plot_average_pulse(runs2,chs,'average_pulse',True)
hdf.plot_average_pulse(runs2,chs,'filters/filt_noconst',False)
hdf.plot_average_pulse(runs2,chs,'filters/filt_noconst',True)


