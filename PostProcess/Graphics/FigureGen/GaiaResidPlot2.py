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

upperLimit = 0.045
lowerLimit = -0.04
plt.ion()
#def DECAM_ResidPlot(upperLimit, lowerLimit):
#plt.style.use('plot_style.txt')
plt.figure(figsize=(7, 4.5))
plt.rc('font',family='serif')
plt.rc('font',serif='Times New Roman')
plt.rc('ytick', labelsize=15)    

df = pd.read_table('GaiaSynthObsMatchTrimmedMappedReorderNoCalspec.dat', sep='\s+')

G = df['G'].to_numpy()
gResid=-1.0*(df['GVeg'].to_numpy() - df['G'].to_numpy())
RpResid=-1.0*(df['GRpVeg'].to_numpy() - df['RP'].to_numpy())
BpResid=-1.0*(df['GBpVeg'].to_numpy() - df['BP'].to_numpy())
gErr=df['Gerr'].to_numpy()
RPerr=df['RPerr'].to_numpy()
BPerr=df['BPerr'].to_numpy()

mask = np.isnan(G)
idgood = np.where(mask==False)

G = G[idgood]
gResid = gResid[idgood]
RpResid = RpResid[idgood]
BpResid = BpResid[idgood]
gErr = gErr[idgood]
RPerr = RPerr[idgood]
BPerr = BPerr[idgood]

spanG = len(G)
xG = np.arange(spanG)
xRP = xG + spanG*1.25
xBP = xRP + spanG*1.25

idsort = np.argsort(G)

upZap = None
lowZap = None
#    upZap, lowZap = setLims(gResid[idsort], upperLimit, lowerLimit, upperLimit*0.8, lowerLimit*0.8)
plt.errorbar(xG, gResid[idsort], yerr=gErr[idsort], uplims=lowZap, lolims=upZap, marker='o', ls='', elinewidth=0.5, capsize=0.5, capthick=0.5, color='g', mec='g', markersize=3)

#    upZap, lowZap = setLims(RpResid, upperLimit, lowerLimit, upperLimit*0.8, lowerLimit*0.8)
plt.errorbar(xRP, RpResid[idsort], yerr=RPerr[idsort], uplims=lowZap, lolims=upZap, marker='o', ls='',elinewidth=0.5, capsize=0.5, capthick=0.5, color='r', mec='r', markersize=3)

#    upZap, lowZap = setLims(BpResid, upperLimit, lowerLimit, upperLimit*0.8, lowerLimit*0.8)
plt.errorbar(xBP, BpResid[idsort], yerr=BPerr[idsort],  uplims=lowZap, lolims=upZap, marker='o', ls='',elinewidth=0.5 ,capsize=0.5, capthick=0.5, color='b', mec='b', markersize=3)

plt.plot([np.min(xG), np.max(xBP)],[0,0],ls='--',lw=0.5, color='k')
plt.xticks([])
plt.ylim(lowerLimit, upperLimit)
plt.ylabel('Gaia Obs - Synth (mag)', fontsize=16)
plt.text(np.mean(xG), upperLimit*0.88, 'G', style='italic', fontsize =
         16, color='g', horizontalalignment='center')
t= plt.text(np.mean(xRP), upperLimit*0.88, 'RP', style='italic', fontsize =
         16, color='r', horizontalalignment='center')
t.set_bbox(dict(facecolor='white', alpha=0.8, edgecolor='white'))
t = plt.text(np.mean(xBP), upperLimit*0.88, 'BP', style='italic', fontsize =
         16, color='b', horizontalalignment='center')
t.set_bbox(dict(facecolor='white', alpha=0.8, edgecolor='white'))


m, s = wghtedStats(gResid, gErr)
plt.plot([np.min(xG), np.max(xG)],[m, m],ls='--',lw=0.5, color='g')
t = plt.text(np.mean(xG), lowerLimit*0.70, 'Mean %.3f\nRMS %.3f' % (m, s), fontsize=10, horizontalalignment='center', color='g')
t.set_bbox(dict(facecolor='white', alpha=0.8, edgecolor='white'))

m, s = wghtedStats(RpResid, RPerr)
print(RpResid)
print(RPerr)
print(m, s)
plt.plot([np.min(xRP), np.max(xRP)],[m, m],ls='--',lw=0.5, color='r')
t = plt.text(np.mean(xRP), lowerLimit*0.70, 'Mean %.3f\nRMS %.3f' % (m, s), fontsize=10, horizontalalignment='center', color='r')
t.set_bbox(dict(facecolor='white', alpha=0.8, edgecolor='white'))

m, s = wghtedStats(BpResid, BPerr)
plt.plot([np.min(xBP), np.max(xBP)],[m, m],ls='--',lw=0.5, color='b')
t = plt.text(np.mean(xBP), lowerLimit*0.70, 'Mean %.3f\nRMS %.3f' % (m, s), fontsize=10, horizontalalignment='center', color='b')
t.set_bbox(dict(facecolor='white', alpha=0.8, edgecolor='white'))


plt.tight_layout()
plt.show()
plt.savefig('GaiaResids3.pdf')
