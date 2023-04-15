import numpy as np
import re
import os

def nameConverter(mapFileName, inFileName, outFileName, nameCol=1):
    nameMap = np.loadtxt(mapFileName, dtype=str)
    numNames = len(nameMap)
    nameDict = {}
    for i in range(numNames):
        oldName = nameMap[i,0]
        newName = nameMap[i,1]
        nameDict[oldName]=newName

    print(nameDict)

    inTable = np.loadtxt(inFileName, dtype=str)

    for line in inTable:
        name = line[nameCol-1]
        try:
            newName = nameDict[name]
        except KeyError:
            print('name ', name, ' not in map. Unchanged')
            newName = name

        line[nameCol-1] = newName

    np.savetxt(outFileName, inTable, fmt='%s')

def nameConvertDirectory(directoryName, mapFileName):
    nameMap = np.loadtxt(mapFileName, dtype=str)
    numNames = len(nameMap)
    nameDict = {}
    for i in range(numNames):
        oldName = nameMap[i,0]
        newName = nameMap[i,1]
        nameDict[oldName]=newName

    namePat = re.compile('(.*)_MedianSed_vac.dat')
    files = os.listdir(directoryName)
    for f in files:
        match = namePat.match(f)
        if match:
            newName = nameDict[match.group(1)] + '_SED.dat'
            print(f, newName)
            os.rename(os.path.join(directoryName, f), os.path.join(directoryName, newName))
    
    

    
    
