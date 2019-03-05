import khe_hdf5 as h
import h5py
import numpy as np

hdf = h5py.File("/Users/hashi/local/data/TMU_2018U/run0400_mass.hdf5",'r')
print hdf['chan1'].keys()
ch=1

cutarray=[['good',1,None],
          ['beam',1,None],
          ['grouptrig',-1,None],
          ['energy',6900,7000]]

cut=h.get_cut(hdf,ch,cutarray)

val1=h.get_attr_np(hdf, ch, 'grouptrig')
val2=h.get_attr_np(hdf, ch, 'energy', cut)

edges_ene=np.arange(6900,6950,1)
edges_ptm=np.arange(8880,8950,1)
h.plot_trend_scat([hdf],[ch],'pretrig_mean',cut)
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
