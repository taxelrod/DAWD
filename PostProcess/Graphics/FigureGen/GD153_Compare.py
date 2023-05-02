# Run from .../WDmodel/WDmodel

import plotSEDsigma as pss
import numpy as np
import matplotlib.pyplot as plt

def makePlot():
    gd153_calspec011=np.loadtxt('RunAll/CRNL_7_paper/MedianSeds/CALSPECcompare/GD153_Calspec_011_Interp.dat')
    gd153_calspec010=np.loadtxt('RunAll/CRNL_7_paper/MedianSeds/CALSPECcompare/GD153_Calspec_010_Interp.dat')
    gd153_calspec011[:,1] = pss.convertSEDtoAB(gd153_calspec011[:,0], gd153_calspec011[:,1])
    gd153_calspec010[:,1] = pss.convertSEDtoAB(gd153_calspec010[:,0], gd153_calspec010[:,1])
    gd153_median = np.loadtxt('RunAll/CRNL_7_paper/DataRelease/SED/GD153_SED.dat')
    gd153_median[:,1] = pss.convertSEDtoAB(gd153_median[:,0], gd153_median[:,1])
    plt.plot(gd153_median[:,0], gd153_calspec011[:,1]-gd153_median[:,1], 'r')
    plt.plot(gd153_median[:,0], gd153_calspec010[:,1]-gd153_median[:,1], 'b')
    plt.ylim(-0.02,0.02)
    plt.xlim(2000, 16000)
    plt.ylabel('AB Mag')
    plt.xlabel('Wavelength (Angstrom)')
    #title('GD153 CALSPEC_010(b) and -011(r) - Median Sed')
    plt.show()
