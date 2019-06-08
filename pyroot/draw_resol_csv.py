import pylab as plt
import pandas as pd

csvdir="./csv/"

template="spilloff"
#template="spillon"

linenames=["CrKAlpha",
           "CoKAlpha",
           "CuKAlpha"]

fnames_beamon  = [csvdir+"%s_beamon_%s.csv"%(l,template) for l in linenames]
fnames_beamoff = [csvdir+"%s_beamoff_%s.csv"%(l,template) for l in linenames]

df_beamon  = [pd.read_csv(f,header=None) for f in fnames_beamon]
df_beamoff = [pd.read_csv(f,header=None) for f in fnames_beamoff]

for dfon,dfoff in zip(df_beamon,df_beamoff):
    dfon.columns=["pre","post","resol","err"]
    dfoff.columns=["pre","post","resol","err"]


def plt_df(df,pre,ax=None,label=''):
    df_tmp=df.query('pre == %d'%pre)
    x=df_tmp.iloc[:,1]
    y=df_tmp.iloc[:,2]
    yerr=df_tmp.iloc[:,3]
    if ax:
        ax.errorbar(x,y,yerr=yerr,fmt='--o',label=label)
    
pres=[0,50,100,150,200]

plt.ion()
plt.close('all')

for dfon,dfoff,line in zip(df_beamon,df_beamoff,linenames):
    fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(6,8))
    ax = axs[0]
    for pre in pres:
        plt_df(dfon,pre,ax,label="pre=%d"%(pre))
    ax.set_ylim(5.0,8.0)
    ax.set_xlabel("post cut")
    ax.set_ylabel("FWHM (eV)")
    ax.set_title('run 397 beam on %s template %s'%(line,template))
    ax.legend(loc='upper left')

    ax = axs[1]
    for pre in pres:
        plt_df(dfoff,pre,ax,label="pre=%d"%(pre))
    ax.set_ylim(4.5,6.5)
    ax.set_xlabel("post cut")
    ax.set_ylabel("FWHM (eV)")
    ax.set_title('run 397 beam off %s template %s'%(line,template))
    ax.legend(loc='upper left')

    plt.tight_layout()
    plt.show()
    plt.savefig("./fig/run397_resol_pre_post_cuts_%s_%s.pdf"%(line,template))

