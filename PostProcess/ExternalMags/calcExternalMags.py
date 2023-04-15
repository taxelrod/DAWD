#!/usr/bin/env python
"""
Load chain data from completed run.  Calculate mags and uncertainties for external filters.

"""

import json
import importlib
import ioWD
import os
import sys
import numpy as np
from astropy.table import Table
import WDmodel
import h5py

# objectPhotometry encapsulates all photometric data and fit results for a WD

        
class objectPhotometry(object):
    def __init__(self, objName=None, paramDict=None):
        if objName is None:
            return
        
        self.objName = objName
        self.paramDict = paramDict

        self.CRNLbandName = paramDict['CRNL']['bandName']

        objList = paramDict['objList']
        if objList is None:
            print('objectPhotometry init: objlist not in paramDict')
        elif objName not in objList:
            print('objectPhotometry init: ', objName, 'not in objlist')
        else:
            self.objParams = objList[objName]
            self.photFileName = self.objParams['photFile']

        self.grid_file = paramDict['grid_file']
        self.grid_name = None
        self.model = WDmodel.WDmodel(self.grid_file, self.grid_name)

        self.synMags = {}

    def loadHSTPhotometry(self):
        self.phot = ioWD.get_phot_for_obj(self.objName, self.photFileName)
        # ignore bands we don't want
        try:
            ignore_pb = self.paramDict['ignore_band']
        except KeyError:
            ignore_pb = None

        if ignore_pb is not None:
            pbnames = self.phot.pb
            pbnames = list(set(pbnames) - set(ignore_pb))
            # filter the photometry recarray to use only the passbands we want
            useind = [x for x, pb in enumerate(self.phot.pb) if pb in pbnames]
            useind = np.array(useind)
            self.phot = self.phot.take(useind)

            # set the pbnames from the trimmed photometry recarray to preserve order
            pbnames = list(self.phot.pb)


        # initialize self.pb
        self.nBands = len(self.phot.pb)
        self.pb = passband.get_pbmodel(self.phot.pb,self.model,  None)

        # set iCRNL to flag band that CRNL is applied to
        pbNames = np.array(list(self.pb.keys()))  # pbNames is an odict, hence this ugliness
        self.iCRNL = np.where(pbNames==self.CRNLbandName)[0][0]  

    def calcSed(self, teff, logg, Av):
        self.modelSed = self.model._get_model(teff, logg)
        self.modelSed = self.model.reddening(self.model._wave, modelSed, Av)
        return self.model._wave, modelSed

    def calcHSTsynMags(self, teff, logg, Av, deltaZp, CRNL):
        self.teff = teff
        self.logg = logg
        self.Av = Av
        self.photCRNL = np.copy(self.phot['mag'])
        self.photCRNL[self.iCRNL] += CRNL[1]*(self.photCRNL[self.iCRNL] - CRNL[0])
        self.modelSed = self.model._get_model(teff, logg)
        self.modelSed = self.model.reddening(self.model._wave, self.modelSed, Av)
        sedPack = np.rec.array([self.model._wave, self.modelSed], dtype=[('wave', '<f8'), ('flux', '<f8')])
                        # needed due to GN interface inconsistency
        self.HSTsynMags = passband.get_model_synmags(sedPack, self.pb) # recarray with dtype=[('pb', '<U5'), ('mag', '<f8')])
        self.HSTsynMags['mag'] += deltaZp
        self.optDM = np.sum((self.photCRNL-self.HSTsynMags['mag'])/self.phot['mag_err']**2)/np.sum(1./self.phot['mag_err']**2)
        self.HSTsynMags['mag'] += self.optDM

    def calcExtSynMag(self, filter, filterBandpass, filterZP):
#        dbug = open('debug.out', 'w')
        self.spectrum = pysynphot.ArraySpectrum(self.model._wave, self.modelSed, waveunits='angstrom', fluxunits='flam')
        self.obs = pysynphot.Observation(self.spectrum, filterBandpass, binset=self.model._wave)

#        obsWave = self.obs.wave
#        obsFlux = self.obs.flux
#        binWave = self.obs.binwave
#        binFlux = self.obs.binflux
#        print(filter, filterZP, self.optDM, file=dbug)
#        for i in range(len(obsWave)):
#            print(obsWave[i], obsFlux[i], binWave[i], binFlux[i], filterBandpass.sample(obsWave[i]), file=dbug)
#        dbug.close()
#        raise ValueError

        synmag = -2.5*np.log10(self.obs.effstim('flam')) - filterZP + self.optDM
        self.synMags[filter] = synmag
#        print('calculating mags for:', self.objName, self.teff, self.logg, self.optDM, filter, synmag) 

class objectCollectionPhotometryHST(object):

    def __init__(self, paramDict):

        self.paramDict = paramDict
        self.objNames = list(paramDict['objList'])  # keeps it pickleable
        self.nObj = len(self.objNames)
        self.nObjParams = 3  # increase this if additional per-object variables need to be added, eg Rv
        self.objPhot = {}
        self.objSlice = {}
 
        print(self.objNames)
        for (i, objName) in enumerate(self.objNames):
            iLo = i*self.nObjParams
            iHi = iLo + self.nObjParams
            self.objSlice[objName] = np.s_[iLo:iHi]
            self.objPhot[objName] = objectPhotometry(objName, paramDict)
            self.objPhot[objName].loadHSTPhotometry()
            if i==0:
                self.nBands = len(self.objPhot[objName].pb)
            else:
                checkNbands = len(self.objPhot[objName].pb)
                assert checkNbands == self.nBands


        self.ZpSlice = np.s_[iHi:iHi+self.nBands-1]  # last element of deltaZp is not explicitly carried because deltaZp sume to 0
        self.CRNLSlice = np.s_[-1:] # CRNL1
        self.nParams = self.nObj*self.nObjParams + self.nBands - 1 + 1 # yes, I know!
        self.CRNL0fixed = paramDict['CRNL']['fixed0']

    def calcHSTsynMags(self, theta):
        # unpack theta into self.nObj arrays of 3, to be interpreted by each objPhot + an array of length self.nBands, which
        # becomes deltaZp + an array of length 2, which is CRNL

        self.theta = theta

        # expand deltaZp by two elements.  The last element, deltaZPF160W is forced to always be zero.
        # the second to last element enforces sum(deltaZp)=0
        
        deltaZp = np.resize(theta[self.ZpSlice], (self.nBands)) # extend by one element
        deltaZp[-1] = -np.sum(theta[self.ZpSlice]) # F160W

        CRNL1 = theta[self.CRNLSlice]  # one element array
        CRNL = np.array([self.CRNL0fixed, CRNL1], dtype = object)
        
        for objName in self.objNames:
            objSlice = self.objSlice[objName]
            obj = self.objPhot[objName]
            (teff, logg, Av) = theta[objSlice]
            obj.calcHSTsynMags(teff, logg, Av, deltaZp, CRNL) # also sets obj.modelSED and self.optDM
     

class objectCollectionPhotometryExt(objectCollectionPhotometryHST):


    def loadFilters(self, filterDict):
        # for each filter set, read in the passbands, then interpolate the m onto model grid
        self.filterFiles = {}
        self.filterBandpass = {}
        self.filterZP = {}
        self.filterMags = {}
        self.filterNames = list(filterDict['filterList'].keys())
        self.nFilters = len(self.filterNames)
        for filter in self.filterNames:
            self.filterFiles[filter] = filterDict['filterList'][filter]['filterFile']
            self.filterZP[filter] = -2.5*np.log10(float(filterDict['filterList'][filter]['zpABflux']))
            try:
                self.filterBandpass[filter] = pysynphot.FileBandpass(self.filterFiles[filter])
                print(filter, self.filterFiles[filter])
            except:
                print(filter, self.filterFiles[filter], 'not found')

    def setIndices(self):
        self.lenPerObject = 4 + self.nFilters  # teff, logg, Av, DM + mag for each filter
        self.offsetTeff = 0
        self.offsetLogg = 1
        self.offsetAv = 2
        self.offsetDM = 3
        # set slice for each object
        self.objSynSlice = {}
        for i, objName in enumerate(self.objNames):
            iLo = i*self.lenPerObject
            iHi = iLo + self.lenPerObject
            self.objSynSlice[objName] = np.s_[iLo:iHi]


    def calcExtMags(self):
        # calculates mags in filter set for theta already set by call to objectCollectionPhotometryHST.calcHSTSynMags()
        #
        self.setIndices()
        for objName in self.objNames:
            for filter in self.filterNames:
                self.objPhot[objName].calcExtSynMag(filter, self.filterBandpass[filter], self.filterZP[filter])

    def addMagData(self, dataDest):
        # iterate through objects.  Fill each object's slice with Teff, eetc and synmags
        temp = np.zeros(self.lenPerObject)
        
        for objName in self.objNames:
            temp[self.offsetTeff] = self.objPhot[objName].teff
            temp[self.offsetLogg] = self.objPhot[objName].logg
            temp[self.offsetAv] = self.objPhot[objName].Av
            temp[self.offsetDM] = self.objPhot[objName].optDM
            for n in range(self.nFilters):
                temp[self.offsetDM + n + 1] = self.objPhot[objName].synMags[self.filterNames[n]]
            dataDest[self.objSynSlice[objName]] = temp.copy()
            

        
# setupPhotEnv sets the environment variable PYSYN_CDBS prior to importing bandpass

def setupPhotEnv(pbPath):
    global passband, pysynphot
    
    if pbPath is None:
#        pbPath = '/home/tsa/Dropbox/WD/PyWD/WDmodel/WDdata/photometry/synphot/'
# for CE1035
        pbPath = '/home/tim/Dropbox/WD/PyWD/WDmodel/WDdata/photometry/synphot/'
# for windows
#        pbPath = 'C:\\Users\\TandR\\Dropbox\\WD\\PyWD\\WDmodel\\WDdata\\photometry\\synphot\\'

    os.environ['PYSYN_CDBS'] = pbPath

    passband = importlib.import_module('passband')
    pysynphot = importlib.import_module('pysynphot')



def dumpSed(objectName, paramFileName, sedFileName, teff, logg, Av, DM):

    setupPhotEnv(None)

    f = open(paramFileName, 'r')
    paramDict = json.load(f)
    f.close()

    objPhot = objectPhotometry(objectName, paramDict)

    wave, sed = objPhot.calcSed(teff, logg, Av)

    norm = 10**(-0.4*DM)

    sedOut = open(sedFileName, 'w')
    for (i, waveLen) in enumerate(wave):
        print(waveLen, norm*sed[i], file=sedOut)

    sedOut.close()
    
    
def calcExternalMags(paramDict, objCollectionExt, hdfChainFileName):
    
    # open the hdf file and get the chain points
    f = h5py.File(hdfChainFileName,'r')
    runData = f['chain']
    runPosition = runData['position']
    (nPts, nObjnParams) = runPosition.shape

    try:
        indexStart = paramDict['indexStart']
    except:
        indexStart = 0

    try:
        indexFinish = paramDict['indexFinish']
    except:
        indexFinish = nPts

    # allocate ouput array for external mag data
    
    nObjects = objCollectionExt.nObj
    nFilters = objCollectionExt.nFilters
    lenPerObject = 4 + nFilters
                     
    extMagData = np.zeros((indexFinish - indexStart, lenPerObject*nObjects))
    for i, n in enumerate(np.arange(indexStart, indexFinish)):
        theta = runPosition[n, :]  # parameter vector for all objects together
        objCollectionExt.calcHSTsynMags(theta)    # calculates model mags in HST bands and the DM for each object
        objCollectionExt.calcExtMags() # calculate model mags in external filters
        objCollectionExt.addMagData(extMagData[i, :])

    f.close()
    return extMagData

def writeHDF5(objCollectionExt, extMagData, hdfOutputFileName):
    f = h5py.File(hdfOutputFileName, 'a')

    grpTables = f.create_group('Tables')

    # create object index table
    nObjects = objCollectionExt.nObj
    nFilters = objCollectionExt.nFilters
    
    objTable = np.empty((nObjects), dtype=[('objName',h5py.string_dtype()), ('indexLo', 'i4'), ('indexHi', 'i4')])
    for (i, objName) in enumerate(objCollectionExt.objNames):
        objTable[i]['objName'] = objName
        objTable[i]['indexLo'] = objCollectionExt.objSynSlice[objName].start
        objTable[i]['indexHi'] = objCollectionExt.objSynSlice[objName].stop

    print(objTable)
    grpTables.create_dataset('objectTable', data=objTable)

    # create variable index table
    varTable = np.empty((4 + nFilters), dtype=[('varName',h5py.string_dtype()), ('index', 'i4')])
    varTable[0] = ('Teff', objCollectionExt.offsetTeff)
    varTable[1] = ('logg', objCollectionExt.offsetLogg)
    varTable[2] = ('Av', objCollectionExt.offsetAv)
    varTable[3] = ('DM', objCollectionExt.offsetDM)
    for (i, filterName) in enumerate(objCollectionExt.filterNames):
        varTable[i+4] = (filterName, i+4)

    print(varTable)
    grpTables.create_dataset('varTable', data=varTable)

    grpData = f.create_group('Data')
    grpData.create_dataset('ChainData', data=extMagData)
    
    f.close()
    
    # create variable index table
                          

def main(paramFileName, hdfChainFileName, filterFileName, hdfOutputFileName, pbPath = None):

    setupPhotEnv(pbPath)  # note that this imports passband

    f = open(paramFileName, 'r')
    paramDict = json.load(f)
    f.close()

    objCollectionExt = objectCollectionPhotometryExt(paramDict)

    # get filters for which to calculate mags; filterFile should be json

    f = open(filterFileName, 'r')
    filterDict = json.load(f)
    f.close()

    objCollectionExt.loadFilters(filterDict)

    # calculate mags and uncertainties in filters


    extMags = calcExternalMags(paramDict, objCollectionExt, hdfChainFileName)

    writeHDF5(objCollectionExt, extMags, hdfOutputFileName)

    return extMags

    # output extMags as hdf file

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
