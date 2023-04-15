import matplotlib.pyplot as plt
import astropy.units as u
import pysynphot as S
import numpy as np
import re
import json
import os
import sys
import matplotlib.ticker
import matplotlib.gridspec as gs
import matplotlib.pyplot as plt
from matplotlib import rc
from mpl_toolkits.axes_grid1 import make_axes_locatable

rc('text', usetex=True)
rc('axes', unicode_minus=True)
#rc('text.latex', preamble=r'\usepackage{apjfonts}')
rc('font',**{'family':'serif','serif':['Times New Roman']})
plt.rcParams['ytick.labelsize']='large'

filter_names = ['F275W', 'F336W', 'F475W', 'F625W', 'F775W','F160W',]


def calcEfflam(SED, bp):
    source = S.ArraySpectrum(SED[:,0], SED[:,1], waveunits='angstrom', fluxunits='flam')
    obs = S.Observation(source, bp)
    efflam = obs.efflam()
    return efflam

def loadHSTphotometryFromTex(texFileName):
    photTable = np.loadtxt(texFileName, dtype=str, comments=('\\', '%'), delimiter='&')
    photDict = {}

    (nrows, ncols) = photTable.shape

    photRe=re.compile('\s*(\d+\.\d+)\(\d+\)\s*')
    
    for n in range(nrows):
        bandMags = np.zeros((ncols-1))
        objName = photTable[n,0].strip()
        for i in range(1, ncols):
            match = photRe.match(photTable[n, i])
            phot = float(match.group(1))
            bandMags[i-1] = phot
        photDict[objName] = bandMags

    return photDict

def convertSEDtoAB(SED):
    flm = (u.erg/u.second/u.cm/u.cm/u.angstrom)
    fluxJy = flm.to(u.Jy, SED[:,1], equivalencies=u.spectral_density(np.array(SED[:,0])*u.angstrom))
    fluxAB = -2.5*np.log10(fluxJy) + 2.5*np.log10(3631)
    SED_AB = np.zeros_like(SED)
    SED_AB[:,0] = SED[:,0]
    SED_AB[:,1] = fluxAB
    return SED_AB

def loadFilters(filterDefFileName):

    f = open(filterDefFileName, 'r')
    filterDict = json.load(f)
    f.close()

    filterFiles = {}
    filterBandpass = {}
    filterZP = {}
    filterNames = list(filterDict['filterList'].keys())
    nFilters = len(filterNames)
    for filter in filterNames:
        filterFiles[filter] = filterDict['filterList'][filter]['filterFile']
        filterBandpass[filter] = S.FileBandpass(filterFiles[filter])

    return filterBandpass

def processSedDirectory(directoryName):

    namePat = re.compile('(.*)_SED.dat')
    files = os.listdir(directoryName)

    sedDict = {}
    
    for f in files:
        match = namePat.match(f)
        if match:
            objName = match.group(1)
            SEDfile = os.path.join(directoryName, f)
            SED = np.loadtxt(SEDfile)
            sedDict[objName] = convertSEDtoAB(SED)

    return sedDict

def calcSedShifts(sedDict, magSpacing):

    nSed = len(sedDict)
    refFlux = np.zeros((nSed))
    shifts = np.zeros((nSed))
    names = np.zeros((nSed), dtype='<U12')

    shiftDict = {}

    for i, obj in enumerate(sedDict.keys()):
        refFlux[i] = sedDict[obj][-1][1]
        names[i] = obj

    fluxMin = np.amin(refFlux)
    fluxMax = np.amax(refFlux)
    refShifted = np.arange(fluxMin, fluxMin + nSed*magSpacing, magSpacing)

    namesOrdered = np.zeros((nSed), dtype='<U12')
    refFluxOrdered = np.zeros_like(refFlux)
    
    idx = np.argsort(refFlux)
    for i in idx:
        refFluxOrdered[i] = refFlux[idx[i]]
        namesOrdered[i] = names[idx[i]]

    for i, obj in enumerate(namesOrdered):
        shift = refFluxOrdered[i] - refShifted[i]
        shiftDict[obj] = shift
        print(obj, i, refFluxOrdered[i], refShifted[i], shift)

    return namesOrdered, refFluxOrdered, shiftDict
        
        
def plotHST(sedDict, filterDefFileName, HSTphotTexFileName, shiftDict):

    filterDict = loadFilters(filterDefFileName)
    photDict = loadHSTphotometryFromTex(HSTphotTexFileName)

    for obj in sedDict.keys():
        phot = photDict[obj]
        for i, bp in enumerate(filterDict.keys()):
            efflam = calcEfflam(sedDict[obj], filterDict[bp])
            photBp = photDict[obj][i]
            plt.plot(efflam, photBp-shiftDict[obj], '*', color='r')
    

def plotSedDict(sedDict, shiftDict):


    ymax = 0
    for obj in sedDict.keys():
        SED = sedDict[obj]
        shift = shiftDict[obj]
        print(obj, shift)
        plt.semilogx(SED[:,0], SED[:,1] - shift, color='k')
        plt.text(SED[-1,0], SED[-1,1] - shift, obj, rotation = -30.0, va='top', fontsize='x-small')
        ymax = max(ymax, SED[-1,1] - shift)

    plt.ylim(ymax,9)
    

def calcTweakedShift(objNamesOrdered, shiftDict, tweakDict=None):

    if tweakDict is None:
        return shiftDict

    newShiftDict = shiftDict.copy()

    nObj = len(shiftDict)

    for it in tweakDict.items():
        j = nObj
        for i, name in enumerate(objNamesOrdered):
            if name==it[0]:
                j=i
                print('matched ', it[0], it[1], j)
        for k in range(j, nObj):
            # find shift for that object in shiftDict
            newShiftDict[objNamesOrdered[k]] += it[1]
            print('tweaking ', k, it[1])

    for name in objNamesOrdered:
        print(name, shiftDict[name], newShiftDict[name])

    return newShiftDict
    
                         
def initTweak():
    tweakDict ={'WDFS1930-52': -1.0,
                'WDFS1557+55': -1.7,
                'WDFS0458-56': -1.0,
                'WDFS1110-17': -0.3,
                'WDFS0727+32': -1.3,
                'WDFS0639-57': -1.0,
                'WDFS0122-30': -1.0,
                'WDFS1111+39': -0.5,
                'WDFS0815+07': -0.7}

    return tweakDict
    
def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)
#def makeFigure2():

sedDirName = 'SED'
HSTphot = 'HSTphot.tex'
HSTfilters = 'ExtFiltersHST.json'

tweakDict = initTweak()

sedDict = processSedDirectory(sedDirName)

namesOrdered, refFluxOrdered, shiftDict = calcSedShifts(sedDict, 0.5)

tshiftDict = calcTweakedShift(namesOrdered, shiftDict, tweakDict)
colors='violet,deepskyblue,forestgreen,tomato,orange,chocolate'
colors = colors.split(',')

plt.ion()
fig = plt.figure(figsize=(8.5, 9.5))
ax = fig.add_subplot(1,1,1)



#ymax = 0
for obj in sedDict.keys():
    SED = sedDict[obj]
    shift = tshiftDict[obj]
    #print(obj, shift)
    ax.semilogx(SED[:,0], SED[:,1] - shift, color='k')
    speclab = '{} $({:+.1f})$'.format(obj, shift)
    ax.text(SED[-1,0], SED[-1,1] - shift, speclab, rotation = -20.0,
            va='top', fontsize='x-small')
    if obj == 'WDFS0815+07':
        line = ax.semilogx(SED[:,0], SED[:,1] - shift, color='k')
        legend = plt.legend([line[0],],\
                            [r'$\mathbf{F}_{\nu} \left( T_{\mathrm{eff}}, \log g, A_V, R_V=3.1, \mu \right)$',],\
                            frameon=False,\
                            loc='lower left',\
                            fontsize='large' )
        ax.add_artist(legend)
    #ymax = max(ymax, SED[-1,1] - shift)

#plt.ylim(ymax,9)
    
filterDict = loadFilters(HSTfilters)
photDict = loadHSTphotometryFromTex(HSTphot)

for obj in sedDict.keys():
    phot = photDict[obj]
    for i, bp in enumerate(filterDict.keys()):
        efflam = calcEfflam(sedDict[obj], filterDict[bp])
        photBp = photDict[obj][i]
        if obj == 'WDFS0815+07':
            plt.plot(efflam, photBp-tshiftDict[obj], 'o', ms = 6,
                     alpha = 0.7,
                 color=colors[i],label=r'$\textit{%s}$'%filter_names[i])
        elif obj =='GD153' or obj == 'GD71' or obj == 'G191B2B':
            plt.plot(efflam, photBp-tshiftDict[obj], '*', ms = 8,
                     color=colors[i], alpha = 0.7)
        else:
            plt.plot(efflam, photBp-tshiftDict[obj], 'o', ms = 6,
                     color=colors[i], alpha = 0.7)


    
xticks = [1000, 2000, 3000, 5000, 7000, 9000, 11000, 13000, 15000, 18000, 21000, 25000]
ax.set_xticks(xticks)
ax.set_xticklabels([str(x) for x in xticks], rotation=45, fontsize='large')
ax.set_xlim(1200, 57000)
plt.minorticks_off()
ax.set_xlabel(r'Wavelength ($\lambda$, \AA)', fontsize='xx-large')


# setup y-axis
yticks = np.arange(9, 40.5, 1)
ax.set_yticks(yticks)
ax.set_yticklabels(['%d'%y for y in yticks], fontsize='large')
ax.set_ylabel(r'AB (mag)', fontsize='xx-large')
ymin, ymax = ax.get_ylim()
ax.set_ylim(ymax+0.75, ymin)
ax2 = ax.twinx()
#plotSedDict(sedDict, tshiftDict)

#plotHST(sedDict, HSTfilters, HSTphot, tshiftDict)
# set twin scale 
Jy = lambda AB: 10.**(-0.4*(AB-2.5*np.log10(3631)))
ymin, ymax = ax.get_ylim()
ax2.plot([],[])
ax2.set_yscale('log')

subs = []
for exp in np.arange(-13, 0, 1.):
    this_pow = []
    for x in [10., 8., 6., 4., 2.]:
        this_pow.append(x*(10**exp))
    subs += this_pow[::-1]
subs = sorted(np.unique(subs))
ax2.set_yticks(subs)
ax2.set_ylabel(r"$\mathbf{F}_{\nu}$ (Jy)", fontsize='xx-large')
ax2.set_ylim((Jy(ymin),Jy(ymax)))

ax3 = ax.twinx()
ax3.spines["right"].set_position(("axes", 1.15))
make_patch_spines_invisible(ax3)
ax3.spines["right"].set_visible(True)
cgs = lambda JY: 10.**(-23)*JY
ymin, ymax = ax2.get_ylim()
#ax3.plot([],[])
ax3.set_yscale('log')
subs = []
for exp in np.arange(-35, -22, 1.):
    this_pow = []
    for x in [10., 8., 6., 4., 2.]:
        this_pow.append(x*(10**exp))
    subs += this_pow[::-1]
subs = sorted(np.unique(subs))
ax3.set_yticks(subs)
ax3.set_ylim((cgs(ymin),cgs(ymax)))
ax3.set_ylabel(r"$\mathbf{F}_{\nu}$ (ergs cm$^{-2}$ s$^{-1}$ Hz$^{-1}$)", fontsize='xx-large')
ax3.patch.set_alpha(0.)


ax.legend(frameon=False, fontsize='large', numpoints=1)
ax.set_xlim(1200, 60000)
fig.tight_layout(rect=[0,0,1,1])
plt.savefig('WD_SED_sequence_all.pdf')
