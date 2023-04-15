# IPython log file

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
    idgood = np.where(mask==False)
    residsGood = resids[idgood]
    xResidsGood = xResids[idgood]
    errsGood = errs[idgood]
    upZap = None
    lowZap = None
#    upZap, lowZap = setLims(residsGood, upperLimit, lowerLimit, upperLimit*0.8, lowerLimit*0.8)
    plt.errorbar(xResidsGood, residsGood, yerr=errsGood, uplims=lowZap, lolims=upZap, marker='o', ls='',elinewidth=0.5, capsize=0.5, capthick=0.5, color=clr, mec=clr, markersize=3)



upperLimit = 0.04
lowerLimit = -0.08
plt.ion()
#def DECAM_ResidPlot(upperLimit, lowerLimit):
#plt.style.use('plot_style.txt')
plt.figure(figsize=(7, 4.5))
plt.rc('font',family='serif')
plt.rc('font',serif='Times New Roman')
plt.rc('ytick', labelsize=15)
df = pd.read_table('DECAM_Legacy.dat', sep=' & ')

G = df['g'].to_numpy()
gResid=-1.0*(df['Sg'].to_numpy() - df['g'].to_numpy())
rResid=-1.0*(df['Sr'].to_numpy() - df['r'].to_numpy())
zResid=-1.0*(df['Sz'].to_numpy() - df['z'].to_numpy())


gErr=0.001*df['gerr'].to_numpy()
rErr=0.001*df['rerr'].to_numpy()
zErr=0.001*df['zerr'].to_numpy()

mask = np.isnan(gResid)
idgood = np.where(mask==False)

G = G[idgood]
gResid = gResid[idgood]
rResid = rResid[idgood]
zResid = zResid[idgood]

gErr = gErr[idgood]
rErr = rErr[idgood]
zErr = zErr[idgood]

spanG = len(G)
xG = np.arange(spanG)
xR = xG + spanG*1.25
xZ = xR + spanG*1.25

idsort = np.argsort(G)


plotResids(xG, gResid[idsort], gErr[idsort], upperLimit, lowerLimit, 'm')
plotResids(xR, rResid[idsort], rErr[idsort], upperLimit, lowerLimit, 'b')
plotResids(xZ, zResid[idsort], zErr[idsort], upperLimit, lowerLimit, 'r')

plt.plot([np.min(xG), np.max(xZ)],[0,0],ls='--',lw=0.5, color='k')
plt.xticks([])
plt.ylim(lowerLimit, upperLimit)
plt.ylabel('DECaLS Obs - Synth (mag)', fontsize=16)

plt.text(np.mean(xG), upperLimit*0.88, 'g',style='italic', fontsize =
         16, color='m', horizontalalignment='center')
plt.text(np.mean(xR), upperLimit*0.88, 'r', style='italic', fontsize =
         16, color='b', horizontalalignment='center')
plt.text(np.mean(xZ), upperLimit*0.88, 'z', style='italic', fontsize =
         16, color='r', horizontalalignment='center')

m, s = wghtedStats(gResid[idsort], gErr[idsort])
plt.plot([np.min(xG), np.max(xG)],[m, m],ls='--',lw=0.5, color='m')
t = plt.text(np.mean(xG), lowerLimit*0.85, 'Mean %.3f\nRMS %.3f' % (m, s),fontsize=10, horizontalalignment='center', color='m')
t.set_bbox(dict(facecolor='white', alpha=0.8, edgecolor='white'))

m, s = wghtedStats(rResid[idsort], rErr[idsort])
plt.plot([np.min(xR), np.max(xR)],[m, m],ls='--',lw=0.5, color='b')
t = plt.text(np.mean(xR), upperLimit*0.35, 'Mean %.3f\nRMS %.3f' % (m, s),fontsize=10, horizontalalignment='center', color='b')
t.set_bbox(dict(facecolor='white', alpha=0.8, edgecolor='white'))

m, s = wghtedStats(zResid[idsort], zErr[idsort])
plt.plot([np.min(xZ), np.max(xZ)],[m, m],ls='--',lw=0.5, color='r')
t = plt.text(np.mean(xZ), upperLimit*0.35, 'Mean %.3f\nRMS %.3f' % (m, s),fontsize=10, horizontalalignment='center', color='r')
t.set_bbox(dict(facecolor='white', alpha=0.8, edgecolor='white'))

plt.savefig('DECAM_Resids_Legacy.pdf')
