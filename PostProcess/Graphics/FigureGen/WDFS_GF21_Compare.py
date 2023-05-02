
# run from .../Tables/GentileFusillo21

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def makePlot():
    df=pd.read_table('WDFS_GF21_matched.dat', sep='\s+')
    plt.figure(figsize=(7, 4.5))
    plt.rc('font',family='serif')
    plt.rc('font',serif='Times New Roman')
    plt.rc('ytick', labelsize=15)

    f, (ax1, ax2, ax3) = plt.subplots(3,1)
    ax1.errorbar(df['logg_med'].to_numpy(), df['loggH'].to_numpy(), xerr=df['logg_sigma'].to_numpy(), yerr=df['e_loggH'].to_numpy(), marker='o', ls='',elinewidth=0.5, capsize=0.5, capthick=0.5, markersize=3)
    ax1.plot([7.0, 9.5],[7.0,9.5], ls='--')
    ax1.set_xlabel('logg (WDFS)')
    f, (ax1, ax2, ax3) = plt.subplots(3,1)
    f.subplots_adjust(hspace=0.5)
    ax1.set_xlabel('logg (WDFS)')
    f.subplots_adjust(hspace=0.6)
    ax1.set_ylabel('logg (GF21)')
    ax1.errorbar(df['logg_med'].to_numpy(), df['loggH'].to_numpy(), xerr=df['logg_sigma'].to_numpy(), yerr=df['e_loggH'].to_numpy(), marker='o', ls='',elinewidth=0.5, capsize=0.5, capthick=0.5, markersize=3)
    ax1.plot([7.0, 9.5],[7.0,9.5], ls='--')
    ax1.set_xlim(7.0,8.2)
    ax1.set_ylim(7.0,8.2)
    ax2.errorbar(df['teff_med'].to_numpy(), df['TeffH'].to_numpy(), xerr=df['teff_sigma'].to_numpy(), yerr=df['e_TeffH'].to_numpy(), marker='o', ls='',elinewidth=0.5, capsize=0.5, capthick=0.5, markersize=3)
    ax2.plot([20000, 70000],[20000,70000], ls='--')
    ax2.set_xlabel('Teff (WDFS)')
    ax2.set_ylabel('Teff (GF21)')
    ax3.errorbar(df['Av_med'].to_numpy(), df['meanAV'].to_numpy(), xerr=df['Av_sigma'].to_numpy(), yerr=df['AVrange'].to_numpy()*0.5, marker='o', ls='',elinewidth=0.5, capsize=0.5, capthick=0.5, markersize=3)
    ax3.set_ylabel('Av (GF21)')
    ax3.set_xlabel('Av (WDFS)')
    ax3.plot([0,0.5],[0,0.5], ls='--')

    plt.show()
