import khe_hdf5 as h
import h5py
import numpy as np
import matplotlib.pyplot as plt

runs=np.asarray(np.loadtxt("He4.runs"),dtype=np.int32).tolist()
#runs=[401,402,403]
runs2=range(136,140)
hdf=h.khe_hdf5(runs)
ch=1

#BADCHAN =  [3, 9, 39, 77, 83, 85, 111, 337, 367, 375, 423, 117, 203, 233, 5, 177]

cutarray=[['good',1,None],
          ['beam',1,None],
          ['grouptrig',-1,None],
#          ['energy',6450,6480]]
          ['energy',6000,6500]]

cutarray2=[['good',1,None],
          ['beam',1,None],
          ['grouptrig',-1,None],
#          ['energy',6450,6480]]
          ['rows_until_next_external_trigger',75,80]]

chs=range(1,480,2)
bins=np.arange(65,90,0.2)
hdf.plot_hist(runs, chs, attr='rows_until_next_external_trigger_clock2',edges=bins,cutarray=cutarray)
plt.savefig('fig/timing2.png')
plt.close()
hdf.plot_hist(runs, chs, attr='rows_until_next_external_trigger_clock3',edges=bins,cutarray=cutarray)
plt.savefig('fig/timing3.png')
plt.close()

#ebins=np.arange(5000,9000,1)
ebins=np.arange(6300,6500,2)
hdf.plot_hist(runs, chs, attr='energy',edges=ebins,cutarray=cutarray2)
plt.savefig('fig/energy.png')
plt.close()

#hdf.plot_hist(runs, chs, attr='rows_after_spill',edges=bins,cutarray=cutarray)
#plt.savefig('fig/spill.png')

#h.plot_trend_scat([hdf],[ch],'pretrig_mean',cut)
#h.plot_hist2d_cut([hdf],[ch],'energy','pretrig_mean',cut,edges_ene,edges_ptm)
#h.plot_hist_cut([hdf],[ch],'energy',cut,6900,6950,1)

'''
'average_pulse',
'beam',
'calibration',
'cuts',
'energy', 
'filt_phase', 
'filt_phase_corr', 
'filt_value', 
'filt_value_dc', 
'filt_value_phc', 
'filt_value_tdc', 
'filters', 
'good', 
'grouptrig', 
'min_value', 
'noise_autocorr', 
'noise_psd', 
'peak_index', 
'peak_region_max', 
'peak_region_mean', 
'peak_region_sum', 
'peak_value', 
'phase_corrector_n', 
'phase_corrector_x', 
'phase_corrector_y', 
'postpeak_deriv', 
'pretrig_mean', 
'pretrig_rms', 
'promptness', 
'pulse_average', 
'pulse_rms', 
'rise_time', 
'rowcount', 
'rows_after_last_external_trigger_refined', 
'rows_after_last_external_trigger_refined_shifted', 
'rows_until_next_external_trigger_refined', 
'rows_until_next_external_trigger_refined_shifted', 
'sec_enemax', 
'sec_enemean', 
'sec_fdmax', 
'sec_fdmean', 
'sec_pr_max', 
'sec_pr_maxmean', 
'sec_pr_maxmin', 
'sec_pr_mean', 
'sec_pr_meanmean', 
'sec_pr_meanmin', 
'sec_pr_sum', 
'shift1', 
'sprmc', 
'timestamp'
'''
