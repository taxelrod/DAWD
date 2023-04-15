"""
Note:  before running need "from photMCMC import objectPhotometry, objectCollectionPhotometry"
"""
import h5py
import pickle
import corner
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from numpy import quantile

def makeQuantiles(hdfFileName, objPickleName, quantiles, nSkip=None, quantOutFileName=None):

    if quantOutFileName is None:
        quantOutFileName = 'testQuant.dat'

    quantOutFile = open(quantOutFileName, 'w')
    
    f=h5py.File(hdfFileName,'r')
    runData = f['chain']

    fpkl=open(objPickleName, 'rb')
    objRun = pickle.load(fpkl)
    fpkl.close()

    objNames = objRun.objNames

    bandNames = objRun.objPhot[objNames[0]].pb.keys()
    bands = {}
    for i, bandName in enumerate(bandNames):
        bands[bandName] = i

    nBands = len(bands)
    
    if objNames is None:
        objNames = runData['objnames']
        
    runPosition = runData['position']
    (nPts, nObjnParams) = runPosition.shape

    #
    # make slices
    #
    if nSkip is None:
        quantSlice = np.s_[:]
    else:
        quantSlice = np.s_[-nSkip:]

    objParams = {'Teff':0, 'logg':1, 'Av':2}
    nObjParams = len(objParams)
    
    for i,obj in enumerate(objNames):
        objIdx = i*nObjParams
        objSlice = np.s_[quantSlice,objIdx:objIdx+nObjParams]
        objChain = runPosition[objSlice]
        objQuants = quantile(objChain, quantiles, axis=0)
        print(obj, file = quantOutFile)
        print(objQuants, file = quantOutFile)
        

    #return runPosition[objSlice['gd71']]


def makePlots(hdfFileName, objPickleName, nPlot=None, outFileName=None, mapOutFileName=None, mapFileName=None, singleObjName=None,
              dumpMapThetaName=None, dumpMedianThetaName=None, noCorners=False):

    if outFileName is None:
        outFileName = 'test.pdf'
    
    f=h5py.File(hdfFileName,'r')
    runData = f['chain']

    fpkl=open(objPickleName, 'rb')
    objRun = pickle.load(fpkl)
    fpkl.close()

    if mapFileName is not None:
        nameMap = np.loadtxt(mapFileName, dtype=str)
        numNames = len(nameMap)
        nameDict = {}
        for i in range(numNames):
            oldName = nameMap[i,0]
            newName = nameMap[i,1]
            nameDict[oldName]=newName
    else:
        nameDict = None
        

    objNames = objRun.objNames

    bandNames = objRun.objPhot[objNames[0]].pb.keys()
    bands = {}
    for i, bandName in enumerate(bandNames):
        bands[bandName] = i

    nBands = len(bands)
    
    if objNames is None:
        objNames = runData['objnames']
        
    runPosition = runData['position']
    runLnProb = runData['lnprob']
    runMagerr = runData['magerr']
    (nPts, nObjnBands) = runMagerr.shape

    #
    # identify MAP point
    #
    mapIdx = np.argmax(runLnProb)
    mapTheta = runPosition[mapIdx, :] # should be nObj*nParams in length ??  + nBands - 1 + 2
    mapMagerr = runMagerr[mapIdx, :]
    print(mapIdx, runLnProb[mapIdx], mapTheta, mapMagerr)

    if dumpMapThetaName is not None:
        np.savetxt(dumpMapThetaName, mapTheta)
    #
    # make slices
    #
    if nPlot is None:
        plotSlice = np.s_[:]
    else:
        plotSlice = np.s_[-nPlot:]

    objParams = {'Teff':0, 'logg':1, 'Av':2}
    nObjParams = len(objParams)

    assert(nObjnBands%nBands == 0)
    nObj = nObjnBands/nBands

    bandSlice = {}
    objSlice = {}
    objMagSlice = {}
    for band in bands.keys():
        indx = bands[band]
        bandSlice[band] = np.s_[indx::nBands]
    
    for i,obj in enumerate(objNames):
        objIdx = i*nObjParams
        bandIdx = i*nBands
        objSlice[obj] = np.s_[plotSlice,objIdx:objIdx+nObjParams]
        objMagSlice[obj] = np.s_[plotSlice, bandIdx:bandIdx+nBands]

    pdfOut = PdfPages(outFileName)

    if mapOutFileName is not None:
        fMap = open(mapOutFileName, 'w')
    else:
        fMap = None

    #
    # make object corner plots
    #
    mapPt = {}
    bandSigmas = {}
    bandMedians = {}
    paramSigmas = {}
    paramMedians = {}

    bandLabels = [r'$\delta$275', r'$\delta$336', r'$\delta$475', r'$\delta$625', r'$\delta$775', r'$\delta$160']
    
    for obj in objNames:
        # if there's a name mapping, use it to get the figure title
        if nameDict is not None:
            titleName = nameDict[obj]
        else:
            titleName = obj
            
        doCorner = False

        if singleObjName is not None and titleName == singleObjName:
            doCorner = True

        if singleObjName is None:
            doCorner = True

        if noCorners:
            doCorner = False
        
        cornerData = np.hstack((runPosition[objSlice[obj]], runMagerr[objMagSlice[obj]]))
        mapPts = np.zeros((nObjParams+nBands))
        mapPts[0:nObjParams] = mapTheta[objSlice[obj][1]]
        mapPts[nObjParams:] = mapMagerr[objMagSlice[obj][1]]
        mapPt[obj] = mapPts
        bandSigmas[obj] = np.std(runMagerr[objMagSlice[obj]], axis=0)
        bandMedians[obj] = np.median(runMagerr[objMagSlice[obj]], axis=0)
        paramSigmas[obj] = np.std(runPosition[objSlice[obj]], axis=0)
        paramMedians[obj] = np.median(runPosition[objSlice[obj]], axis=0)
        if doCorner == True:
            fig=corner.corner(cornerData,labels=list(objParams)+bandLabels, show_titles=True, label_kwargs={'fontsize': 13}, title_kwargs={'fontsize': 13, 'rotation':45.0, 'ha':'center'}, title_fmt='.1e' )

            # resize the figure so that 45 deg text doesn't get cut off
            scale = 0.9
            axes = fig.axes
            for ax in axes:
                bbox = ax.get_position()
                (x0, y0, x1, y1) = bbox.bounds
                ax.set_position(pos=[x0*scale, y0*scale, x1*scale, y1*scale])

            fig.text(0.7,0.95,titleName,fontsize=16)

            plt.savefig(pdfOut, format='pdf')
            plt.close(fig)
    
    #
    # make band delta zp  + CRNL corner plots
    #
    nBandsM1 = nBands - 1
    cornerData = runPosition[plotSlice, -nBandsM1-1:]
    #
    # kluge alert!
    #
    zpF160W = -np.sum( runPosition[plotSlice, -nBandsM1-1:-1], axis=1 )
#    print('zpF160W:', runPosition[plotSlice, -nBandsM1-1:-1].shape, cornerData.shape, np.mean(zpF160W))

    (nRow, nCol) = cornerData.shape
    
    cornerDatax = np.hstack((cornerData[:,-nBandsM1-1:-1], np.reshape(zpF160W, (nRow, 1)), np.reshape(cornerData[:,-1], (nRow, 1))))
    
    zpLabels = []
    for (i,bandName) in enumerate(bandNames):
        if i == nBandsM1:
            break
        zpLabels.append('zp' + bandName)

    zpLabels.append('zpF160W')
    zpLabels.append('CRNL1')

    zpLabels = [r'$\Delta$275', r'$\Delta$336', r'$\Delta$475', r'$\Delta$625', r'$\Delta$775', r'$\Delta$160', r'$\alpha_{CRNL}$']
    
    if not noCorners:
        fig=corner.corner(cornerDatax,labels=zpLabels, show_titles=True, label_kwargs={'fontsize': 13}, title_kwargs={'fontsize': 13, 'rotation':45.0, 'ha':'center'}, title_fmt='.1e'  )

        # resize the figure so that 45 deg text doesn't get cut off
        scale = 0.9
        axes = fig.axes
        for ax in axes:
            bbox = ax.get_position()
            (x0, y0, x1, y1) = bbox.bounds
            ax.set_position(pos=[x0*scale, y0*scale, x1*scale, y1*scale])

        plt.savefig(pdfOut, format='pdf')
        plt.close(fig)

    pdfOut.close()

# ----------------------------------------
# Form a theta from the product of the median object params + median zp and CRNL
#
    if dumpMedianThetaName is not None:
        medianTheta = np.zeros_like(mapTheta)
        for obj in objNames:
            medianTheta[objSlice[obj][1]] = paramMedians[obj]

        zp_CRNL_medians = np.median(runPosition[plotSlice, -nBandsM1-1:], axis=0)
        medianTheta[-nBandsM1-1:] = zp_CRNL_medians

        np.savetxt(dumpMedianThetaName, medianTheta)
# ----------------------------------------
    print('# quantity median std')
    for (i, label) in enumerate(zpLabels):
        print(label, np.median(cornerDatax[:, i]), np.std(cornerDatax[:, i]))

    if fMap is not None:
        fileHdrLine = '# obj teff_map logg_map Av_map '
        for bandName in bandNames:
            fileHdrLine += 'r' + bandName + '_map '
        fileHdrLine += 'teff_med logg_med Av_med '
        for bandName in bandNames:
            fileHdrLine += 'r' + bandName + '_med '
        fileHdrLine += 'teff_sigma logg_sigma Av_sigma '
        for bandName in bandNames:
            fileHdrLine += 'r' + bandName + '_sigma '
        print(fileHdrLine, file=fMap)
        
        for obj in objNames:
            print(obj, end=' ', file=fMap)
            for n, v in enumerate(mapPt[obj]):
                print(v, end=' ', file=fMap)
            for n, v in enumerate(paramMedians[obj]):
                print(v, end=' ', file=fMap)
            for n, v in enumerate(bandMedians[obj]):
                print(v, end=' ', file=fMap)
            for n, v in enumerate(paramSigmas[obj]):
                print(v, end=' ', file=fMap)
            for n, v in enumerate(bandSigmas[obj]):
                print(v, end=' ', file=fMap)
            print(' ', file=fMap)
            
        print('# quantity median std', file=fMap)
        for (i, label) in enumerate(zpLabels):
            print(label, np.median(cornerDatax[:, i]), np.std(cornerDatax[:, i]), file=fMap)
            
        fMap.close()
    f.close()

