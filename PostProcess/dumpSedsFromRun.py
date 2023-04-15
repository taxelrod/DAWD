import numpy as np
import photMCMC as pmc

def dumpSedsFromRun(runFileName, paramFileName, outFileNameSuffix):
    runDtype=np.dtype([('obj','<U18'), ('teff', float), ('logg', float), ('Av', float), ('DM', float)])
    colList = (0, 1, 2, 3, 4)
    runData = np.loadtxt(runFileName, dtype=runDtype, usecols=colList, skiprows=1)
    for objData in runData:
        obj = objData['obj']
        teff = objData['teff']
        logg = objData['logg']
        Av = objData['Av']
        DM = objData['DM']
        sedFileName = obj + outFileNameSuffix
        pmc.dumpSed(obj, paramFileName, sedFileName, teff, logg, Av, DM)
    return runData
