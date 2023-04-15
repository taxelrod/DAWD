import numpy as np

import pandas as pd
import matplotlib.pyplot as plt

def CRNL_Uncorrect(zpFileName, hstMagsFileName, outputFileName):
    
    dfzp = pd.read_table(zpFileName, sep='\s+')
    CRNL1 = dfzp[dfzp['quantity']=='CRNL1']['median'].to_numpy()[0]
    zpF275W = dfzp[dfzp['quantity']=='zpF275W']['median'].to_numpy()[0]
    zpF336W = dfzp[dfzp['quantity']=='zpF336W']['median'].to_numpy()[0]
    zpF475W = dfzp[dfzp['quantity']=='zpF475W']['median'].to_numpy()[0]
    zpF625W = dfzp[dfzp['quantity']=='zpF625W']['median'].to_numpy()[0]
    zpF775W = dfzp[dfzp['quantity']=='zpF775W']['median'].to_numpy()[0]
    zpF160W = dfzp[dfzp['quantity']=='zpF160W']['median'].to_numpy()[0]

    zps = np.array([zpF275W, zpF336W, zpF475W, zpF625W, zpF775W, zpF160W])

    print(CRNL1, zps)

    dfhst = pd.read_table(hstMagsFileName, sep='\s+')

    objIds = dfhst['obj'].to_numpy()

    of = open(outputFileName, 'w')

    print('# obj F275W uF275W F336W uF336W F475W uF475W F625W uF625W F775W uF775W F160W uF160W', file=of)
    
    for objId in objIds:
        hstMags = getHSTmags(dfhst, objId)
        umags = uncorrectHSTmags(hstMags, zps, CRNL1)
        print(objId, hstMags, umags)
        print(objId, end=' ', file=of)
        for i in range(len(hstMags)):
            print(hstMags[i], end=' ', file=of)
            print(umags[i], end=' ', file=of)
        print('', file=of)

    of.close()
    
def CRNL_Recorrect(zpFileName, hstMagsFileName, CRNL1):
    pass

def getHSTmags(df, objId):
    print(objId)
    dfRow = df[df['obj']==objId]
    F275W = dfRow['F275W'].to_numpy()[0]
    F336W = dfRow['F336W'].to_numpy()[0]
    F475W = dfRow['F475W'].to_numpy()[0]
    F625W = dfRow['F625W'].to_numpy()[0]
    F775W = dfRow['F775W'].to_numpy()[0]
    F160W = dfRow['F160W'].to_numpy()[0]

    
    return [F275W, F336W, F475W, F625W, F775W, F160W]

# These are the AB mags that would be observed with the filter zeropoints
# corrected, and CRNL corrected.  These should be compared to the synthetic
# mags from the SED models

def uncorrectHSTmags(mags, zps, CRNL1):
    CRNL0 = 15.0

    umags = mags.copy()
    umags[5] = mags[5]*(1 + CRNL1) - CRNL1*CRNL0
    umags -= zps

    return umags

def correctionHSTmags(mags, zps, CRNL1):
    CRNL0 = 15.0

    correction = np.zeros_like(mags)
    correction[5] = CRNL1*(mags[5] - CRNL0)
    correction -= zps

    return umagscorrection
