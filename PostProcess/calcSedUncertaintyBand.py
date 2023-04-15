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

class objectCollectionPhotometryHST(object):

    def __init__(self, paramDict):

        self.paramDict = paramDict
        self.objNames = list(paramDict['objList'])  # keeps it pickleable
        self.nObj = len(self.objNames)
        self.nObjParams = 3  # increase this if additional per-object variables need to be added, eg Rv
        self.nFilters = 6 # KLUGE
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


    def calcHSTsynMags(self, theta, calcOnly=None):
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
            if calcOnly is not None and objName != calcOnly :
                continue
            objSlice = self.objSlice[objName]
            obj = self.objPhot[objName]
            (teff, logg, Av) = theta[objSlice]
            obj.calcHSTsynMags(teff, logg, Av, deltaZp, CRNL) # also sets obj.modelSED and self.optDM
     

# setupPhotEnv sets the environment variable PYSYN_CDBS prior to importing bandpass

def setupPhotEnv(pbPath):
    global passband, pysynphot
    
    if pbPath is None:
#        pbPath = '/home/tsa/Dropbox/WD/PyWD/WDmodel/WDdata/photometry/synphot/'
# for CE1035
        pbPath = '/home/tim/Dropbox/WD/PyWD/WDmodel/WDdata/photometry/synphot/'
# for windows
        pbPath = 'C:\\Users\\TandR\\Dropbox\\WD\\PyWD\\WDmodel\\WDdata\\photometry\\synphot\\'

    os.environ['PYSYN_CDBS'] = pbPath

    passband = importlib.import_module('passband')
    pysynphot = importlib.import_module('pysynphot')


def calcSedUncertaintyBand(paramDict, objCollection, hdfChainFileName, outFileName):
    
    # open the hdf file and get the chain points
    # NOTE that the paramDict and the hdf file must match!
    
    f = h5py.File(hdfChainFileName,'r')
    runData = f['chain']
    runPosition = runData['position']
    (nPts, nObjnParams) = runPosition.shape
    print(nPts, ' in chain')

    try:
        indexStart = paramDict['indexStart']
    except:
        indexStart = 0

    try:
        indexFinish = paramDict['indexFinish']
    except:
        indexFinish = nPts

    # get object name and reference spectrum for the uncertainty calc

    calcObjectName = paramDict['calcObjectName']
    refSpectrumFile = paramDict['refSpectrum']

    refSpectrum = np.loadtxt(refSpectrumFile)

    # need to check its wl is same as that for the calcObject
                     

    for name in objCollection.objNames:
        obj = objCollection.objPhot[name]
        obj.init = False
        
    for i, n in enumerate(np.arange(indexStart, indexFinish)):
        theta = runPosition[n, :]  # parameter vector for all objects together
        objCollection.calcHSTsynMags(theta, calcOnly=calcObjectName) # calculates seds for each object
 
        for name in objCollection.objNames:
            if name == calcObjectName:
                obj = objCollection.objPhot[name]
                if not obj.init:
                    obj.wl = obj.model._wave
                    nWL = len(obj.wl)
                    obj.sedRef = refSpectrum[:,1]
                    obj.sedDeltaSq = np.zeros((nWL))
                    obj.init = True

                objFlux = obj.modelSed * 10**(-0.4*obj.optDM)
                obj.sedDeltaSq += (objFlux - obj.sedRef)**2
 
    nChainPts = float(indexFinish - indexStart)
    print('num pts used:', nChainPts)
    
    for name in objCollection.objNames:
        if name == calcObjectName:
            obj = objCollection.objPhot[name]
            obj.sedSigma = np.sqrt(obj.sedDeltaSq)/nChainPts
            dataStruct = np.zeros((nWL, 3))
            dataStruct[:,0] = obj.wl
            dataStruct[:,1] = obj.sedRef
            dataStruct[:,2] = obj.sedSigma
            np.savetxt(outFileName, dataStruct)
    
    f.close()


def main(paramFileName, hdfChainFileName, outDirectoryName, pbPath = None):

    setupPhotEnv(pbPath)  # note that this imports passband

    f = open(paramFileName, 'r')
    paramDict = json.load(f)
    f.close()

    objCollection = objectCollectionPhotometryHST(paramDict)

    calcSedUncertaintyBand(paramDict, objCollection, hdfChainFileName, outDirectoryName)

    return objCollection

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
