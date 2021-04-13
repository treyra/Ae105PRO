import OrbitInterpolation
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

def main():
    """
    Test of orbit interpolation software
    """

    #Import generated data from matlab simulation
    tTruth = np.genfromtxt('C:\\Users\\trey_\\Documents\\Research\\DARTS\\Code\\t.csv',delimiter=',')
    xTruth = np.genfromtxt('C:\\Users\\trey_\\Documents\\Research\\DARTS\\Code\\x.csv',delimiter=',')
    paramsTruth = np.genfromtxt('C:\\Users\\trey_\\Documents\\Research\\DARTS\\Code\\params.csv',delimiter=',')

    #Subsample the data for fitting by taking every second (truth has millisecond data)

    tSub = tTruth[0::2000]
    xSub = xTruth[0::2000]

    paramsFit = OrbitInterpolation.orbitInterp(tSub,xSub)
    print("Fit Params: " + str(paramsFit))
    print("True Params: " + str(paramsTruth))

    #Interpolate
    xFit = OrbitInterpolation.keplerainOrbitModel(paramsFit,tTruth)

    #Compare with a spline fit
    xPosSplineFun = interp1d(tSub,xSub[:,0], kind='cubic')
    
    
    #Plot and compare
    plt.figure()

    plt.plot(tTruth,xTruth[:,0])
    plt.plot(tTruth,xFit[:,0])
    plt.plot(tTruth,xPosSplineFun(tTruth))


    plt.xlabel("Time (s)")
    plt.ylabel("X coordinate (km)")
    plt.title("Comparing fit and truth")
    plt.show()

    plt.figure()

    plt.plot(tTruth,xFit[:,0]-xTruth[:,0])
    plt.plot(tTruth,xPosSplineFun(tTruth)-xTruth[:,0])
    plt.show()





if __name__ == "__main__":
    main()