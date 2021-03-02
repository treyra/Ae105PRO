
import SplineOrbitInterpolation
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

#Demonstration that spline fitting of simulated orbit data works better than fitting orbit parameters 
#or even modeling as a Keplarian orbit with the same orbit parameters used to generate the simulated 
#orbit (the simulation includes effects up to J2 and the sun and moon perturbations)


def main():
    #Load sample orbit data (pseudo randomly generated)
    tTruth = np.genfromtxt('C:\\Users\\trey_\\Documents\\Research\\DARTS\\Code\\t.csv',delimiter=',')
    xTruth = np.genfromtxt('C:\\Users\\trey_\\Documents\\Research\\DARTS\\Code\\x.csv',delimiter=',')
    paramsTruth = np.genfromtxt('C:\\Users\\trey_\\Documents\\Research\\DARTS\\Code\\params.csv',delimiter=',')

     #Subsample the data for fitting by taking every second (truth has millisecond data)
    
    tSub = tTruth[0::1000]
    xSub = xTruth[0::1000]

    #Compare with a spline fit
    xSpline = SplineOrbitInterpolation.splineOrbitInterplation(tSub,xSub,tTruth)

    print("Univariate Splines")
    #Return the actual splines for differentiation
    posSplines = SplineOrbitInterpolation.splineOrbitInterplationPosSpline(tSub,xSub)

    #Store interpolated values
    xSpline2 = np.zeros((len(tTruth),3))
    vSpline2 = np.zeros((len(tTruth),3))
    aSpline2 = np.zeros((len(tTruth),3))
    #Evaluate spline and derivatives
    for i in range(3):
        xSpline2[:,i] = posSplines[i](tTruth)
        vSpline2[:,i] = posSplines[i](tTruth,1)
        aSpline2[:,i] = posSplines[i](tTruth,2)

    print("Univariate Splines Complete")

    #Compare with numerical differentiation (Time step is 1 second)
    xSubVCenterDiff = xSub.copy()
    for i in range(3):
        xSubVCenterDiff[:,i+3] = np.gradient(xSub[:,i], 1)
    

    #Spline fit (reusing interpolation method)
    vSpline = SplineOrbitInterpolation.splineOrbitInterplation(tSub,xSubVCenterDiff,tTruth)
    vSpline = vSpline[:,3:]

    #Numerical derivative on the spline fit
    vNumerical = np.zeros((len(xSpline),3))
    aNumerical = np.zeros((len(xSpline),3))
    for i in range(3):
        #millisecond data rate
        vNumerical[:,i] = np.gradient(xSpline[:,i], .001)
        aNumerical[:,i] = np.gradient(vNumerical[:,i], .001)
        
    

    #Compute position and velocity errors as a function of time
    xSplineError = xSpline - xTruth
    vSplineError = vSpline - xTruth[:,3:]
    vNumericalError = vNumerical - xTruth[:,3:]
    xSpline2Error = xSpline2 - xTruth[:,:3]
    vSpline2Error = vSpline - xTruth[:,3:]

    
    #Plot and compare
    
    plt.figure()

    plt.plot(tTruth,1000*np.linalg.norm(xSpline2Error,axis=1),label="Spline")
    

    plt.xlabel("Time (s)")
    plt.ylabel("Position Error (m)")
    plt.title("Position Error in Orbit Univariate Spline Interpolation (Simulation Performed Every Second)")
    plt.legend()
    plt.draw()

    plt.figure()

    plt.plot(tTruth,1000*np.linalg.norm(vSpline2Error,axis=1),label="Spline")
    

    plt.xlabel("Time (s)")
    plt.ylabel("Velocity Error (m/s)")
    plt.title("Velocity Error in Orbit Univariate Spline Interpolation (Simulation Performed Every Second)")
    plt.legend()
    plt.draw()
    
    
    
    ################## OLD
    plt.figure()

    plt.plot(tTruth,1000*np.linalg.norm(xSplineError[:,:3],axis=1),label="Spline")
    

    plt.xlabel("Time (s)")
    plt.ylabel("Position Error (m)")
    plt.title("Position Error in Orbit Spline Interpolation (Simulation Performed Every Second)")
    plt.legend()
    plt.draw()

    plt.figure()

    plt.plot(tTruth,1000*np.linalg.norm(xSplineError[:,3:],axis=1),label="Spline")

    plt.xlabel("Time (s)")
    plt.ylabel("Velocity Error (m/s)")
    plt.title("Velocity Error in Orbit Spline Interpolation (Simulation Performed Every Second)")
    plt.legend()
    plt.draw()

    plt.figure()

    plt.plot(tTruth,1000*np.linalg.norm(vSplineError,axis=1),label="Spline")

    plt.xlabel("Time (s)")
    plt.ylabel("Velocity Error (m/s)")
    plt.title("Velocity Error in Orbit Spline Interpolation on Numerical Derivative (Simulation Performed Every Second)")
    plt.legend()
    plt.draw()

    plt.figure()

    plt.plot(tTruth,1000*np.linalg.norm(vNumericalError,axis=1),label="Spline")

    plt.xlabel("Time (s)")
    plt.ylabel("Velocity Error (m/s)")
    plt.title("Velocity Error in Orbit Spline Interpolation on Numerical Derivative of Spline Position (Simulation Performed Every Second)")
    plt.legend()
    plt.draw()



    plt.figure()

    plt.plot(tTruth,1000*np.linalg.norm(vNumericalError,axis=1),label="Spline on Numerical Derivative of Spline")
    plt.plot(tTruth,1000*np.linalg.norm(vSplineError,axis=1),label="Spline on Numerical Derivative")
    plt.plot(tTruth,1000*np.linalg.norm(xSplineError[:,3:],axis=1),label="Spline")

    plt.xlabel("Time (s)")
    plt.ylabel("Velocity Error (m/s)")
    plt.title("Errors in interpolation methods")
    plt.legend()
    plt.draw()

    plt.figure()

    plt.plot(tTruth,1000*np.linalg.norm(vNumericalError,axis=1),label="Spline on Numerical Derivative of Spline")
    plt.plot(tTruth,1000*np.linalg.norm(xSplineError[:,3:],axis=1),label="Spline")

    plt.xlabel("Time (s)")
    plt.ylabel("Velocity Error (m/s)")
    plt.title("Errors in interpolation methods")
    plt.legend()
    plt.draw()


    plt.figure()

    plt.plot(tTruth,1000*np.linalg.norm(aNumerical,axis=1),label="Spline on Numerical Derivative of Spline")
    plt.plot(tTruth,1000*np.linalg.norm(aSpline2,axis=1),label="Spline")
    plt.plot(tTruth,1000*np.linalg.norm(aNumerical-aSpline2,axis=1),label="Diff")

    plt.xlabel("Time (s)")
    plt.ylabel("Acceleration (m/s^2)")
    plt.title("Differences in interpolation methods")
    plt.legend()
    plt.draw()

    plt.show()








if __name__ == "__main__":
    main()