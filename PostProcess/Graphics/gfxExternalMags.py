import sys
import h5py
import corner
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from numpy import quantile

def makePlots(hdfFileName, nskipPlot=None, outFileName=None, statOutFileName=None, noPlot=None):

    if outFileName is None:
        outFileName = 'test.pdf'

    if statOutFileName is not None:
        fStat = open(statOutFileName, 'w')
    else:
        fStat = sys.stdout
        
    f=h5py.File(hdfFileName,'r')

    tableData = f['Tables']
    objTable = tableData['objectTable']
    varTable = tableData['varTable']

    magData = f['Data']['ChainData']
    nPts = magData.shape[0]

    if noPlot is None:
        pdfOut = PdfPages(outFileName)
        plotLabels = varTable['varName']
        filterNames = plotLabels[4:]

    print('# obj', end=' ', file=fStat)
    for filterVar in varTable[4:]:
#        print(filterVar[0])
        print(filterVar[0] + 'mean', filterVar[0] + 'med', filterVar[0] + 'sigma', end=' ', file=fStat)
    print('', file=fStat)
    
    nObj = objTable.shape[0]
    for n in range(nObj):
        objSlice = np.s_[objTable[n][1]:objTable[n][2]]
        objName = objTable[n][0]
        objData = magData[:, objSlice]
        if noPlot is None:
            fig=corner.corner(objData,labels=plotLabels, show_titles=True)
            fig.text(0.7, 0.95, objName, fontsize=16)
            pdfOut.savefig()
#            plt.savefig(pdfOut, format='pdf')
            plt.close(fig)
        # calculate mag stats for filters
        print(objName, end=' ', file=fStat)
        objStats = {}
        for filterVar in varTable[4:]:
            filterMags = objData[:, filterVar['index']]
            objStats[filterVar['varName']] = (np.mean(filterMags), np.median(filterMags), np.std(filterMags))
            print(np.mean(filterMags), np.median(filterMags), np.std(filterMags), end=' ', file=fStat)
        print('', file=fStat)

#        print(objName, objStats, file=fStat)


    if noPlot is None:
        pdfOut.close()
    f.close()
    fStat.close()
