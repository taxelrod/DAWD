# Calculate AB mag in an external filter from a SED file produced by dumpSED

import os
import pysynphot
import json
import re
import numpy as np

# 1. read in SED file, get DM
# 2. read in filter file - use code from calcExternalMags.py
# 3. use pysynphot to calculate mag

def decodeSEDfile(SEDfileName):
    f = open(SEDfileName, 'r')
    hdr = f.readline()
    f.close()

    Objre = re.compile('^#.*obj:\s+(\w+)')
    Objmatch = Objre.match(hdr)
    ObjId = Objmatch.group(1)

    DMre = re.compile('^#.*DM:\s+(\d+\.\d+)')
    DMmatch = DMre.match(hdr)
    DM = float(DMmatch.group(1))

    SED = np.loadtxt(SEDfileName)
    return SED, ObjId

def loadFilters(filterDict):

    wlArray = np.arange(1000., 20000., 0.1) # for calculating filter ZP
    
    filterFiles = {}
    filterBandpass = {}
    filterZP = {}
    filterNames = list(filterDict['filterList'].keys())
    nFilters = len(filterNames)
    for filter in filterNames:
        filterFiles[filter] = filterDict['filterList'][filter]['filterFile']
        filterBandpass[filter] = pysynphot.FileBandpass(filterFiles[filter])
        try:
            filterZP[filter] = -2.5*np.log10(float(filterDict['filterList'][filter]['zpABflux']))
        except KeyError:
            filterFluxZP, filterZP[filter] = calcFilterABZP(filterBandpass[filter], wlArray)
            print('Calculating ZP from bandpass for ', filter, filterFluxZP, filterZP[filter])
            
        print(filter, filterFiles[filter], filterZP[filter])


    return filterNames, filterBandpass, filterZP


def calcFilterABZP(filterBandpass, wlArray):
    spectrumAB = pysynphot.ArraySpectrum(wlArray, 2.9979e-5*3631/wlArray**2, waveunits='angstrom', fluxunits='flam')
    obs = pysynphot.Observation(spectrumAB, filterBandpass, binset=wlArray)
    fluxZP = obs.effstim('flam')
    synmag = -2.5*np.log10(obs.effstim('flam'))
    return fluxZP, synmag

def calcExtSynMag(filterBandpass, filterZP, SED):
    spectrum = pysynphot.ArraySpectrum(SED[:,0], SED[:,1], waveunits='angstrom', fluxunits='flam')
    obs = pysynphot.Observation(spectrum, filterBandpass, binset=SED[:,0])
    synmag = -2.5*np.log10(obs.effstim('flam')) - filterZP
    return synmag

def setupPhotEnv():
    pass

def main(filterDefFileName, SEDdirectory, outputFileName):
    setupPhotEnv()

    f = open(filterDefFileName, 'r')
    filterDict = json.load(f)
    f.close()

    filterNames, filterBandpass, filterZP = loadFilters(filterDict)
    nFilters = len(filterNames)
    synmag = np.zeros(nFilters)

    outputFile = open(outputFileName, 'w')
    
    # print hdr line
    print('# objId', end=' ', file=outputFile)
    for filter in filterNames:
        print(filter, end=' ', file=outputFile)
    print('', file=outputFile)

    SEDpathRe = re.compile('.*Sed.dat')

    fileList = os.listdir(SEDdirectory)
    for file in fileList:
        if SEDpathRe.match(file):
            filePath = os.path.join(SEDdirectory, file)
            SED, ObjId = decodeSEDfile(filePath)
            print(ObjId, end=' ', file=outputFile)
            for i, filter in enumerate(filterNames):
                synmag[i] = calcExtSynMag(filterBandpass[filter], filterZP[filter], SED)
                print(synmag[i], end= ' ', file=outputFile)

            print('', file=outputFile)


    outputFile.close()

                

