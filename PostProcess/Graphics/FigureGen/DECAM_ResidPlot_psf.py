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

upperLimit = 0.05
lowerLimit = -0.05
plt.ion()
#def DECAM_ResidPlot(upperLimit, lowerLimit):
#plt.style.use('plot_style.txt')
plt.figure(figsize=(7, 4.5))
plt.rc('font',family='serif')
plt.rc('font',serif='Times New Roman')
plt.rc('ytick', labelsize=15)
df = pd.read_table('DECAM_psf_synth_obs_match.dat', sep='\s+')

G = df['DECAM_gmean'].to_numpy()
gResid=-1.0*(df['DECAM_gmean'].to_numpy() - df['wavg_mag_psf_g'].to_numpy())
rResid=-1.0*(df['DECAM_rmean'].to_numpy() - df['wavg_mag_psf_r'].to_numpy())
iResid=-1.0*(df['DECAM_imean'].to_numpy() - df['wavg_mag_psf_i'].to_numpy())
zResid=-1.0*(df['DECAM_zmean'].to_numpy() - df['wavg_mag_psf_z'].to_numpy())
yResid=-1.0*(df['DECAM_Ymean'].to_numpy() - df['wavg_mag_psf_y'].to_numpy())


gErr=df['wavg_magerr_psf_g'].to_numpy()
rErr=df['wavg_magerr_psf_r'].to_numpy()
iErr=df['wavg_magerr_psf_i'].to_numpy()
zErr=df['wavg_magerr_psf_z'].to_numpy()
yErr=df['wavg_magerr_psf_y'].to_numpy()

mask = np.isnan(gResid)
idgood = np.where(mask==False)

G = G[idgood]
gResid = gResid[idgood]
rResid = rResid[idgood]
iResid = iResid[idgood]
zResid = zResid[idgood]
yResid = yResid[idgood]

gErr = gErr[idgood]
rErr = rErr[idgood]
iErr = iErr[idgood]
zErr = zErr[idgood]
yErr = yErr[idgood]

spanG = len(G)
xG = np.arange(spanG)
xR = xG + spanG*1.25
xI = xR + spanG*1.25
xZ = xI + spanG*1.25
xY = xZ + spanG*1.25

idsort = np.argsort(G)


plotResids(xG, gResid[idsort], gErr[idsort], upperLimit, lowerLimit, 'm')
plotResids(xR, rResid[idsort], rErr[idsort], upperLimit, lowerLimit, 'b')
plotResids(xI, iResid[idsort], iErr[idsort], upperLimit, lowerLimit, 'g')
plotResids(xZ, zResid[idsort], zErr[idsort], upperLimit, lowerLimit, 'r')
plotResids(xY, yResid[idsort], yErr[idsort], upperLimit, lowerLimit, 'tab:brown')

plt.plot([np.min(xG), np.max(xY)],[0,0],ls='--',lw=0.5, color='k')
plt.xticks([])
plt.ylim(lowerLimit, upperLimit)
plt.ylabel('DES Obs - Synth (mag)', fontsize=16)
plt.text(np.mean(xG), upperLimit*0.90, 'g', style='italic', fontsize =
         16, color='m', horizontalalignment='center')
plt.text(np.mean(xR), upperLimit*0.90, 'r', style='italic', fontsize =
         16, color='b', horizontalalignment='center')
plt.text(np.mean(xI), upperLimit*0.90, 'i', style='italic', fontsize =
         16, color='g', horizontalalignment='center')
plt.text(np.mean(xZ), upperLimit*0.90, 'z', style='italic', fontsize =
         16, color='r', horizontalalignment='center')
plt.text(np.mean(xY), upperLimit*0.90, 'y', style='italic', fontsize =
         16, color='tab:brown', horizontalalignment='center')

m, s = wghtedStats(gResid[idsort], gErr[idsort])
plt.plot([np.min(xG), np.max(xG)],[m, m],ls='--',lw=0.5, color='m')
t = plt.text(np.mean(xG), lowerLimit*0.85, 'Mean %.3f\nRMS %.3f' %
             (m, s), fontsize=10, horizontalalignment='center', color='m')
t.set_bbox(dict(facecolor='white', alpha=0.8, edgecolor='white'))

m, s = wghtedStats(rResid[idsort], rErr[idsort])
plt.plot([np.min(xR), np.max(xR)],[m, m],ls='--',lw=0.5, color='b')
t = plt.text(np.mean(xR), lowerLimit*0.85, 'Mean %.3f\nRMS %.3f' % (m, s), fontsize=10, horizontalalignment='center', color='b')
t.set_bbox(dict(facecolor='white', alpha=0.8, edgecolor='white'))

m, s = wghtedStats(iResid[idsort], iErr[idsort])
plt.plot([np.min(xI), np.max(xI)],[m, m],ls='--',lw=0.5, color='g')
t = plt.text(np.mean(xI), lowerLimit*0.85, 'Mean %.3f\nRMS %.3f' % (m, s), fontsize=10, horizontalalignment='center', color='g')
t.set_bbox(dict(facecolor='white', alpha=0.8, edgecolor='white'))

m, s = wghtedStats(zResid[idsort], zErr[idsort])
plt.plot([np.min(xZ), np.max(xZ)],[m, m],ls='--',lw=0.5, color='r')
t = plt.text(np.mean(xZ), lowerLimit*0.85, 'Mean %.3f\nRMS %.3f' %
             (m, s), fontsize=10, horizontalalignment='center', color='r')
t.set_bbox(dict(facecolor='white', alpha=0.8, edgecolor='white'))

m, s = wghtedStats(yResid[idsort], yErr[idsort])
plt.plot([np.min(xY), np.max(xY)],[m, m],ls='--',lw=0.5, color='tab:brown')
t = plt.text(np.mean(xY), lowerLimit*0.85, 'Mean %.3f\nRMS %.3f' % (m, s), fontsize=10, horizontalalignment='center', color='tab:brown')
t.set_bbox(dict(facecolor='white', alpha=0.8, edgecolor='white'))

plt.tight_layout()
plt.show()
plt.savefig('DECAM_Resids_psf2.pdf')
