import numpy as np
import re

def tableMerge(table1Name, table2Name, mapping, tableOutName, delim1=None, delim2=None):
    # load table1 and 2 as string arrays (accept any format)
    # decode mapping
    # make new table with nrows = table1, ncols = len(mapping)
    table1 = np.loadtxt(table1Name, dtype=str, delimiter=delim1)
    table2 = np.loadtxt(table2Name, dtype=str, delimiter=delim2)

    nrows = table1.shape[0]
    ncols1 = table1.shape[1]
    ncols2 = table2.shape[1]
    
    ncols = len(mapping)
    tableOut = np.zeros((nrows, ncols), dtype=table1.dtype)

    for n in range(nrows):
        if table1[n,0].strip() != table2[n,0].strip():
            print('mismatch: ', table1[n,0], table2[n,0])
        for j in range(ncols1):
            table1[n, j] = cleanNumString(table1[n, j])
        for j in range(ncols2):
            table2[n, j] = cleanNumString(table2[n, j])
            
        
    for i,map in enumerate(mapping):
        tableNum, columnNum = decodeMap(map)
        if tableNum == 1:
            tableOut[:,i] = table1[:, columnNum]
        elif tableNum == 2:
            tableOut[:,i] = table2[:, columnNum]
        else:
            print('oops', tableNum, columnNum)

    np.savetxt(tableOutName, tableOut, fmt='%s', delimiter=' & ')
                        
def decodeMap(mapStr):
    mapRe = re.compile('(\d+)\.(\d+)')
    mapMatch = mapRe.match(mapStr)
    tableNum = mapMatch.group(1)
    columnNum = mapMatch.group(2)
    return int(tableNum), int(columnNum)

def cleanNumString(string):
    fpRe = re.compile('(\d\d\.\d\d\d\d)\d+')
    fpMatch = fpRe.match(string)
    if fpMatch is not None:
        fp = float(fpMatch.group(1))
        return str(np.round(fp, decimals=3))
    else:
        return string
