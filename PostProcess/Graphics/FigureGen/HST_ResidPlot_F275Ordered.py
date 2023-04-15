import numpy as np

import pandas as pd
import matplotlib.pyplot as plt

def wghtedStats(data, errs):
    wghtedMean = np.average(data, weights=1./errs**2)
    wghtedErr = np.sqrt(np.average((data-wghtedMean)**2, weights=1./errs**2))
    return wghtedMean, wghtedErr
    
def setLims(data, upper, lower, upperReplace, lowerReplace):
    idxUpper = np.where(data>upper)
    boolUpper = np.zeros_like(data)
    boolUpper[idxUpper] = True
    data[idxUpper] = upperReplace
    
    idxLower = np.where(data<lower)
    boolLower = np.zeros_like(data)
    boolLower[idxLower] = True
    data[idxLower] = lowerReplace
    
    return boolUpper, boolLower

def plotResids(xResids, resids, errs, upperLimit, lowerLimit, clr):
    mask = np.isnan(resids)
    idsort = np.where(mask==False)
    residsGood = resids[idsort]
    xResidsGood = xResids[idsort]
    errsGood = errs[idsort]
    upZap, lowZap = setLims(residsGood, upperLimit, lowerLimit, upperLimit*0.8, lowerLimit*0.8)
    plt.errorbar(xResidsGood, residsGood, yerr=errsGood, uplims=lowZap, lolims=upZap, marker='o', ls='',elinewidth=0.5, capsize=0.5, capthick=0.5, color=clr, mec=clr, markersize=3)

upperLimit = 0.03
lowerLimit = -0.03
plt.ion()
#def DECAM_ResidPlot(upperLimit, lowerLimit):
#plt.style.use('plot_style.txt')
plt.figure(figsize=(7, 4.5))
plt.rc('font',family='serif')
plt.rc('font',serif='Times New Roman')
plt.rc('ytick', labelsize=15)

df = pd.read_table('HST_resids.dat', sep='\s+')

F275 = df['F275W'].to_numpy()
idsort = np.argsort(F275)

span275=len(F275)

x275 = np.arange(span275)
x336 = x275 + span275*1.25
x475 = x336 + span275*1.25
x625 = x475 + span275*1.25
x775 = x625 + span275*1.25
x160 = x775 + span275*1.25

F275Resid=df['rF275W_med'].to_numpy()
F336Resid=df['rF336W_med'].to_numpy()
F475Resid=df['rF475W_med'].to_numpy()
F625Resid=df['rF625W_med'].to_numpy()
F775Resid=df['rF775W_med'].to_numpy()
F160Resid=df['rF160W_med'].to_numpy()

F275Err=df['dF275W'].to_numpy()
F336Err=df['dF336W'].to_numpy()
F475Err=df['dF475W'].to_numpy()
F625Err=df['dF625W'].to_numpy()
F775Err=df['dF775W'].to_numpy()
F160Err=df['dF160W'].to_numpy()


plotResids(x275, F275Resid[idsort], F275Err[idsort], upperLimit, lowerLimit, 'tab:purple')
plotResids(x336, F336Resid[idsort], F336Err[idsort], upperLimit, lowerLimit, 'b')
plotResids(x475, F475Resid[idsort], F475Err[idsort], upperLimit, lowerLimit, 'g')
plotResids(x625, F625Resid[idsort], F625Err[idsort], upperLimit, lowerLimit, 'tab:orange')
plotResids(x775, F775Resid[idsort], F775Err[idsort], upperLimit, lowerLimit, 'r')
plotResids(x160, F160Resid[idsort], F160Err[idsort], upperLimit, lowerLimit, 'tab:brown')

plt.plot([np.min(x275[idsort]), np.max(x160[idsort])],[0,0],ls='--',lw=0.5, color='k')
plt.xticks([])
plt.ylim(lowerLimit, upperLimit)
plt.ylabel('HST Obs - Synth (mag)', fontsize=16)

plt.text(np.mean(x275[idsort]), upperLimit*0.90, 'F275W', style='italic', fontsize =
         16, color='tab:purple', horizontalalignment='center')
plt.text(np.mean(x336[idsort]), upperLimit*0.90, 'F336W', style='italic', fontsize =
         16, color='b', horizontalalignment='center')
plt.text(np.mean(x475[idsort]), upperLimit*0.90, 'F475W', style='italic', fontsize =
         16, color='g', horizontalalignment='center')
plt.text(np.mean(x625[idsort]), upperLimit*0.90, 'F625W', style='italic', fontsize =
         16, color='tab:orange', horizontalalignment='center')
plt.text(np.mean(x775[idsort]), upperLimit*0.90, 'F775W', style='italic', fontsize =
         16, color='r', horizontalalignment='center')
plt.text(np.mean(x160[idsort])-1, upperLimit*0.90, 'F160W', style='italic', fontsize =
         16, color='tab:brown', horizontalalignment='center')

m, s = wghtedStats(F275Resid[idsort], F275Err[idsort])
plt.text(np.mean(x275[idsort]), lowerLimit*0.95, 'Mean %.3f\nRMS %.3f' % (m, s), fontsize=10, horizontalalignment='center', color='tab:purple')

m, s = wghtedStats(F336Resid[idsort], F336Err[idsort])
plt.text(np.mean(x336[idsort]), lowerLimit*0.95, 'Mean %.3f\nRMS %.3f' % (m, s), fontsize=10, horizontalalignment='center', color='b')

m, s = wghtedStats(F475Resid[idsort], F475Err[idsort])
plt.text(np.mean(x475[idsort]), lowerLimit*0.95, 'Mean %.3f\nRMS %.3f' % (m, s), fontsize=10, horizontalalignment='center', color='g')

m, s = wghtedStats(F625Resid[idsort], F625Err[idsort])
plt.text(np.mean(x625[idsort]), lowerLimit*0.95, 'Mean %.3f\nRMS %.3f' % (m, s), fontsize=10, horizontalalignment='center', color='tab:orange')

m, s = wghtedStats(F775Resid[idsort], F775Err[idsort])
plt.text(np.mean(x775[idsort]), lowerLimit*0.95, 'Mean %.3f\nRMS %.3f'
         % (m, s), fontsize=10, horizontalalignment='center', color='r') 

m, s = wghtedStats(F160Resid[idsort], F160Err[idsort])
plt.text(np.mean(x160[idsort]), lowerLimit*0.95, 'Mean %.3f\nRMS %.3f' % (m, s), fontsize=10, horizontalalignment='center', color='tab:brown')

plt.tight_layout()
plt.show()
plt.savefig('HST_Resids_F275.pdf')
