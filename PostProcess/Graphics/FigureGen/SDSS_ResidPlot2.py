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
#    upZap, lowZap = setLims(residsGood, upperLimit, lowerLimit, upperLimit*0.8, lowerLimit*0.8)
    upZap = None
    lowZap = None
    plt.errorbar(xResidsGood, residsGood, yerr=errsGood, uplims=lowZap, lolims=upZap, marker='o', ls='',elinewidth=0.5, capsize=0.5, capthick=0.5, color=clr, mec=clr, markersize=3)


upperLimit = 0.1
lowerLimit = -0.125
plt.ion()
#def DECAM_ResidPlot(upperLimit, lowerLimit):
#plt.style.use('plot_style.txt')
plt.figure(figsize=(7, 4.5))
plt.rc('font',family='serif')
plt.rc('font',serif='Times New Roman')
plt.rc('ytick', labelsize=15)    

df = pd.read_table('SDSS_SynthObjMatchAll_mapped_reorder_s.dat', sep='\s+')

G = df['SDSSgmean'].to_numpy()
uResid=-1.0*(df['SDSSumean'].to_numpy() - df['SDSScat_umean'].to_numpy())
gResid=-1.0*(df['SDSSgmean'].to_numpy() - df['SDSScat_gmean'].to_numpy())
rResid=-1.0*(df['SDSSrmean'].to_numpy() - df['SDSScat_rmean'].to_numpy())
iResid=-1.0*(df['SDSSimean'].to_numpy() - df['SDSScat_imean'].to_numpy())
zResid=-1.0*(df['SDSSzmean'].to_numpy() - df['SDSScat_zmean'].to_numpy())

uErr=df['SDSScat_usigma'].to_numpy()
gErr=df['SDSScat_gsigma'].to_numpy()
rErr=df['SDSScat_rsigma'].to_numpy()
iErr=df['SDSScat_isigma'].to_numpy()
zErr=df['SDSScat_zsigma'].to_numpy()

mask = np.isnan(uResid)
idgood = np.where(mask==False)

G = G[idgood]
uResid = uResid[idgood]
gResid = gResid[idgood]
rResid = rResid[idgood]
iResid = iResid[idgood]
zResid = zResid[idgood]

uErr = gErr[idgood]
gErr = gErr[idgood]
rErr = rErr[idgood]
iErr = iErr[idgood]
zErr = zErr[idgood]

spanG = len(G)
xU = np.arange(spanG)
xG = xU + spanG*1.25
xR = xG + spanG*1.25
xI = xR + spanG*1.25
xZ = xI + spanG*1.25
idsort = np.argsort(G)


plotResids(xU, uResid[idsort], uErr[idsort], upperLimit, lowerLimit, 'm')
plotResids(xG, gResid[idsort], gErr[idsort], upperLimit, lowerLimit, 'b')
plotResids(xR, rResid[idsort], rErr[idsort], upperLimit, lowerLimit, 'g')
plotResids(xI, iResid[idsort], iErr[idsort], upperLimit, lowerLimit, 'r')
plotResids(xZ, zResid[idsort], zErr[idsort], upperLimit, lowerLimit, 'tab:brown')

plt.plot([np.min(xU), np.max(xZ)],[0,0],ls='--',lw=0.5, color='k')
plt.xticks([])
plt.ylim(lowerLimit, upperLimit)
plt.ylabel('SDSS Obs - Synth (mag)', fontsize=16)

plt.text(np.mean(xU), upperLimit*0.90, 'u', style='italic', fontsize =
         16, color='m', horizontalalignment='center')
plt.text(np.mean(xG), upperLimit*0.90, 'g', style='italic', fontsize =
         16, color='b', horizontalalignment='center')
plt.text(np.mean(xR), upperLimit*0.90, 'r', style='italic', fontsize =
         16, color='g', horizontalalignment='center')
plt.text(np.mean(xI), upperLimit*0.90, 'i', style='italic', fontsize =
         16, color='r', horizontalalignment='center')
plt.text(np.mean(xZ), upperLimit*0.90, 'z', style='italic', fontsize =
         16, color='tab:brown', horizontalalignment='center')
#t.set_bbox(dict(facecolor='white', alpha=0.8, edgecolor='white'))


m, s = wghtedStats(uResid[idsort], uErr[idsort])
plt.plot([np.min(xU), np.max(xU)],[m, m],ls='--',lw=0.5, color='m')
t = plt.text(np.mean(xU), lowerLimit*0.85, 'Mean %.3f\nRMS %.3f' % (m, s), fontsize=10, horizontalalignment='center', color='m')
t.set_bbox(dict(facecolor='white', alpha=0.8, edgecolor='white'))

m, s = wghtedStats(gResid[idsort], gErr[idsort])
plt.plot([np.min(xG), np.max(xG)],[m, m],ls='--',lw=0.5, color='b')
t = plt.text(np.mean(xG), lowerLimit*0.85, 'Mean %.3f\nRMS %.3f' % (m, s), fontsize=10, horizontalalignment='center', color='b')
t.set_bbox(dict(facecolor='white', alpha=0.8, edgecolor='white'))

m, s = wghtedStats(rResid[idsort], rErr[idsort])
plt.plot([np.min(xR), np.max(xR)],[m, m],ls='--',lw=0.5, color='g')
t = plt.text(np.mean(xR), lowerLimit*0.85, 'Mean %.3f\nRMS %.3f' % (m, s), fontsize=10, horizontalalignment='center', color='g')
t.set_bbox(dict(facecolor='white', alpha=0.8, edgecolor='white'))

m, s = wghtedStats(iResid[idsort], iErr[idsort])
plt.plot([np.min(xI), np.max(xI)],[m, m],ls='--',lw=0.5, color='r')
t = plt.text(np.mean(xI)-0.35, lowerLimit*0.85, 'Mean %.3f\nRMS %.3f' % (m, s), fontsize=10, horizontalalignment='center', color='r')
t.set_bbox(dict(facecolor='white', alpha=0.8, edgecolor='white'))

m, s = wghtedStats(zResid[idsort], zErr[idsort])
plt.plot([np.min(xZ), np.max(xZ)],[m, m],ls='--',lw=0.5, color='tab:brown')
t = plt.text(np.mean(xZ), lowerLimit*0.85, 'Mean %.3f\nRMS %.3f' % (m, s), fontsize=10, horizontalalignment='center', color='tab:brown')
t.set_bbox(dict(facecolor='white', alpha=0.8, edgecolor='white'))
    


plt.tight_layout()
plt.show()
plt.savefig('SDSS_Resids3.pdf')
