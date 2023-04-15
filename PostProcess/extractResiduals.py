import numpy as np
import matplotlib.pyplot as plt

def ingestMapFile(mapFileName):
    mapDtype=np.dtype([('obj','<U12'), ('r275', float), ('r336', float), ('r475', float), ('r625', float), ('r775', float), ('r160', float), ('sig275', float), ('sig336', float), ('sig475', float), ('sig625', float), ('sig775', float), ('sig160', float)])
    colList = (0, 4, 5, 6, 7, 8, 9, 22, 23, 24, 25, 26, 27)

    mapData = np.loadtxt(mapFileName, dtype=mapDtype, usecols=colList, skiprows=1)

    wl = [0.275, 0.336, 0.475, 0.625, 0.775, 1.60]

    objResids = {}

    nObj = len(mapData)
    
    resids = np.zeros((nObj,6))
    residSigmas = np.zeros((nObj,6))

    resids[:, 0] = mapData[:]['r275']
    resids[:, 1] = mapData[:]['r336']
    resids[:, 2] = mapData[:]['r475']
    resids[:, 3] = mapData[:]['r625']
    resids[:, 4] = mapData[:]['r775']
    resids[:, 5] = mapData[:]['r160']

    residSigmas[:, 0] = mapData[:]['sig275']
    residSigmas[:, 1] = mapData[:]['sig336']
    residSigmas[:, 2] = mapData[:]['sig475']
    residSigmas[:, 3] = mapData[:]['sig625']
    residSigmas[:, 4] = mapData[:]['sig775']
    residSigmas[:, 5] = mapData[:]['sig160']

    for i in range(nObj):
        objResids[mapData[i]['obj']] = (resids[i, :], residSigmas[i, :])
            
    return wl, objResids, resids, residSigmas

def plotObjResiduals(objResids, objName):
    plt.errorbar(range(6),objResids[objName][0], yerr=objResids[objName][1], fmt='x')
