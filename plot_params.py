import khe_hdf5 as h
import h5py
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

runs=np.asarray(np.loadtxt("He4.runs"),dtype=np.int32).tolist()
runs3=range(136,139)
runs2=range(345,358)
#runs2=range(132,133)
hdf=h.khe_hdf5(runs2+runs3)
ch=1
chs=[1]

cutarray_on=[['good',1,None],
             ['beam',1,None],
             ['grouptrig',-1,None]]


cutarray_off=[['good',1,None],
             ['beam',0,None],
             ['grouptrig',-1,None]]

ecut= [['energy',6900,6975]]

energy_edges=range(5000,9000,20)
params=[['filt_phase',[0.1*x-1 for x in range(20)]],
        ['filt_phase_corr', [0.1*x-1 for x in range(20)]],
#        ['filt_value',  50],
#        ['filt_value_dc',  50],
#        ['filt_value_phc',  50],
#        ['filt_value_tdc',  50],
        ['min_value',  range(8750,9000,4)],
        ['peak_index',  range(285,305,1)],
        ['peak_region_max',  50],
        ['peak_region_mean',  50],
        ['peak_region_sum',  50],
        ['peak_value',  50],
        ['postpeak_deriv',  range(0,50)],
        ['pretrig_mean',  range(8800,9100,4)],
        ['pretrig_rms',  range(0,100)],
        ['promptness',  [0.01*x+0.3 for x in range(20)]],
        ['pulse_average',  50],
        ['pulse_rms',  50],
        ['rise_time',  [2e-6*x+1.5e-4 for x in range(25)]]]

pdfname="fig/plot_params"+h.get_name(runs2+runs3,chs)+".pdf"
print "printing...", pdfname
with PdfPages(pdfname) as pdf:
    for p in params:
        print p[0]
        fig=plt.figure()
        ax1=fig.add_subplot(2,3,1)
        ax1=hdf.plot_hist2d(runs2,chs,'energy',p[0],energy_edges,p[1],cutarray_on,ax1)   
        ax2=fig.add_subplot(2,3,2,sharex=ax1,sharey=ax1)
        ax2=hdf.plot_hist2d(runs2,chs,'energy',p[0],energy_edges,p[1],cutarray_off,ax2)   
        ax2=fig.add_subplot(2,3,3)
        ax2=hdf.plot_hist2d(runs3,chs,'energy',p[0],energy_edges,100,cutarray_on,ax2)   

        ax3=fig.add_subplot(2,3,4)
        ax3,H3=hdf.plot_hist(runs2,chs,p[0],p[1],cutarray_on,ax3)   
        ax3,H3=hdf.plot_hist(runs2,chs,p[0],H3[1],cutarray_off,ax3)   
        ax3,H3=hdf.plot_hist(runs3,chs,p[0],H3[1],cutarray_on,ax3)   

        ax4=fig.add_subplot(2,3,5)
        ax4,H4=hdf.plot_hist(runs2,chs,p[0],p[1],cutarray_on+ecut,ax4)   
        ax4,H4=hdf.plot_hist(runs2,chs,p[0],H4[1],cutarray_off+ecut,ax4)   
        ax4,H4=hdf.plot_hist(runs3,chs,p[0],H4[1],cutarray_on+ecut,ax4)   
        ax3.set_yscale('log')
        ax4.set_yscale('log')
        ax3.set_ylim(bottom=0.1)
        ax4.set_ylim(bottom=0.1)
        ax1.set_title(p[0])
        pdf.savefig()
        plt.close()

#cut=h.get_cut(hdf,ch,cutarray)
#val1=h.get_attr_np(hdf, ch, 'grouptrig')
#val2=h.get_attr_np(hdf, ch, 'energy', cut)
#edges_ene=np.arange(6900,6950,1)
#edges_ptm=np.arange(8880,8950,1)
#h.plot_trend_scat([hdf],[ch],'pretrig_mean',cut)

'''
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
