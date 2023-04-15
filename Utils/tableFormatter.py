import numpy as np

def tableFormatter(inFileName, outFileName, startCol=None, endCol=None):

    inTable = np.loadtxt(inFileName, dtype=str, comments=None)

    fOut = open(outFileName, 'w')

    (nRows, nCol) = inTable.shape

    if startCol==None:
        colRange = range(0, nCol)
        startCol = 0
        endCol = nCol
    elif endCol==None:
        colRange = range(startCol, nCol)
        endCol = nCol
    else:
        colRange = range(startCol, endCol)

    nCol = endCol - startCol
    print(colRange)
    print(nCol)

    # print out latex table header

    print('\\begin{table}', file=fOut)
    print('\\tiny', file=fOut)
    print('\\begin{tabular}{*{%d}{c}}' % nCol, file=fOut)
    columnHeadings = inTable[0][colRange].copy()
    tableHeader=''
    for col in colRange:
        if col==endCol-1:
           tableHeader += '%d \\\\' % endCol
        else:
            tableHeader += '%d & ' % (col+1)
    print(tableHeader, file=fOut)
    
    for line in inTable[1:nRows]:
        for col in colRange:
            field = line[col]
            try:
                fieldPrint = float(field)
                if col < endCol-1:
                    print('%8.3f & ' % fieldPrint, end='', file=fOut)
                else:
                    print('%8.3f ' % fieldPrint, end='', file=fOut)
            except ValueError:
                if col < endCol-1:
                    if field != '\"\"':
                        print('%s & '%field, end='', file=fOut)
                    else:
                        print('& ', end='', file=fOut)
                else:
                    if field != '\"\"':
                        print('%s '%field, end='', file=fOut)
                    else:
                        print(' ', end='', file=fOut)
                    

        print('\\\\', file=fOut)

    # print out latex table footer

    print('\\end{tabular}', file=fOut)
    print('\\end{table}', file=fOut)

    # print out column key table

    print('\\begin{table}', file=fOut)
    print('\\tiny', file=fOut)
    print('\\begin{tabular}{c c}', file=fOut)
    for col in colRange:
        print('%d & %s \\\\ ' % (col+1, columnHeadings[col-startCol]), file=fOut)

    print('\\end{tabular}', file=fOut)
    print('\\end{table}', file=fOut)
   
    
    

    
    fOut.close()

def tableFormatterNoSplit(inFileName, outFileName, sortFirst=True):

    inTable = np.loadtxt(inFileName, dtype=str, comments=None)

    if sortFirst:
        (nRow, nCol) = inTable.shape
        idx = np.argsort(inTable[:,0])
        sortTable = np.zeros_like(inTable)
        for i in range(nRow):
            sortTable[i,:] = inTable[idx[i],:]
        inTable = sortTable

    fOut = open(outFileName, 'w')

    (nRows, nCol) = inTable.shape

    startCol = 0
    endCol = nCol
    colRange = range(startCol, endCol)

    # print out latex table header

    print('\\begin{table*}[h]', file=fOut)
    print('\\scriptsize', file=fOut)
    print('\\begin{centering}', file=fOut)
    print('\\begin{tabular}{*{%d}{c}}' % nCol, file=fOut)

    for line in inTable[0:nRows]:
        for col in colRange:
            field = line[col]
            try:
                fieldPrint = float(field)
                if col < endCol-1:
                    print('%8.3f & ' % fieldPrint, end='', file=fOut)
                else:
                    print('%8.3f ' % fieldPrint, end='', file=fOut)
            except ValueError:
                if col < endCol-1:
                    if field != '\"\"':
                        print('%s & '%field, end='', file=fOut)
                    else:
                        print('& ', end='', file=fOut)
                else:
                    if field != '\"\"':
                        print('%s '%field, end='', file=fOut)
                    else:
                        print(' ', end='', file=fOut)
                    

        print('\\\\', file=fOut)

    # print out latex table footer

    print('\\end{tabular}', file=fOut)
    print('\\end{centering}', file=fOut)
    print('\\end{table*}', file=fOut)

    fOut.close()

    
    
    
    
