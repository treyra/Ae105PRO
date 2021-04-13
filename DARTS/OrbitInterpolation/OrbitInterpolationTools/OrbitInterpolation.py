import numpy as np

from scipy.optimize import least_squares


def orbitInterp(t,x):
    """
    Orbit interpolation function, fits keplarian orbital parameters
    (a,e,i,Omega,omega,nu0) to interpolate between two or more measurement 
    points along an orbit. This is intended for short time steps, so ignores
    all higher order corrections such as J2, and assumes the points given are
    equally accuarte (weighting may come). Points given must be time, and
    state, where the state is a 6 dimensional vector with positon followed by
    velocity
    
    """

    #Since Keplarian orbit paramteters can't be fit directly by a linear
    #regression model, we want to use matlabs built in option to perform a
    #non linear least squares fit
    
    #Some useful constants:
    earthRadius = 6378.136
    
    #We'll start with a cirular orbit with slight inclination (to resolve
    #ambiguities in Omega and omega if i = 0) in LEO
    #(a,e,i,Omega,omega,nu0)
    params0 = np.array([earthRadius + 310 , .001, 0*np.pi/180, 0, 0, 0])

    fitResult = least_squares(fittingResiduals, params0, args=(t, x),xtol = 10**-16)
    print(fitResult)

    return fitResult.x



def fittingResiduals(params,t,x):

    return np.linalg.norm(keplerainOrbitModel(params,t) - x,axis = 1) 


    
def keplerainOrbitModel(params,t):
    """
    Convert from orbital elements to ECI postition and velocity 
    This assumes an elliptical orbit, but could be extended to other types
    Vectorized for the fitting function
    """
    
    

    a = params[0]
    e = params[1]
    inclin = params[2]
    Omega = params[3]
    w = params[4]
    nu0 = params[5] #In Radians
    
    #Assuming Earth
    muEarth = 398600.432896939164493230
    
    #Need to compute the true anomoly nu by solving the Kepler time problem
    #get M0 (Mean anomoly) of nu0
    M0 = KeplerTime(e,nu0)
    
    #Compute mean motion n= 2pi/T
    n = np.sqrt(muEarth/a**3)
    
    #Compute Mean anomoly now
    
    M = M0 + n*t

    #Solve
    nu = KeplerAngle(e,M)
    
    #Compute the postion and veloctiy vectors in perifocal coordinates,
    #assuming an eliptical orbit
    #Magnitudes of the vectors (Elliptical orbit)  
    rmag = a * (1-e**2)/(1+e*np.cos(nu))
    vmag = (muEarth*(2/rmag - 1/a))**(1/2)* 1000 # Want this in m/s, for better weighting of the fit

    rperi = rmag*np.array([np.cos(nu),np.sin(nu),np.zeros(len(nu))])
    vperi = vmag*np.array([-np.sin(nu),np.cos(nu),np.zeros(len(nu))])
    
    #Rotate to ECI
    (r,v) = perifocalToGeoCentric(rperi,vperi,Omega,inclin,w)
    #Stack as a state vector
    x = np.vstack((r,v))
    #Curve fitting we actually want the vector transposed
    x = x.transpose()
    return x

def perifocalToGeoCentric(Rin,Vin,Omega,inclin,w):
    """
    Convert from perifocal frame to ECI
    """
    #Create the rotation matricies
    rotationMatrix = np.matmul(np.matmul(EulerZ(Omega), EulerX(inclin)), EulerZ(w))
    
    R = np.matmul(rotationMatrix , Rin)
    V = np.matmul(rotationMatrix , Vin)
    return (R,V)



def EulerX(theta):
    """
    EulerX Matrix
    """
    Rx = np.array([[1,0,0],[0,np.cos(theta),-np.sin(theta)],[0,np.sin(theta),np.cos(theta)]])
    return Rx


def EulerZ(theta):
    """    
    EulerZ Matrix
    """
    Rz = np.array([[np.cos(theta),-np.sin(theta),0],[np.sin(theta),np.cos(theta),0],[0,0,1]])
    return Rz


def KeplerTime(e,nu):
    """
    Given nu, find M 
    """
    E = 2 * np.tan(np.sqrt((1-e)/(1+e))*np.tan(nu/2))
    M = E - e * np.sign(E)
    return M



#def KeplerAngle(e,M):
#    """
#    Given M, find nu: (Using newton's method) 
#    """
#    #Pick initial guess for E using Curtis's Guidlines    
#    E = M - e/2;
#    print(M)
#    if( M < np.pi):
#            E = M + e/2;
#    #Set Eold negative to garuntee we enter loop
#    Eold = -1;
#    #Set the convergence bound requested
#    while(np.abs(E-Eold) > 10**-15):
#        Eold = E;
#        E = E - (E - e * np.sin(E) - M)/(1 - e * np.cos(E));
#    nu = 2 * np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E/2))  ;
#    return nu

#Trying a vectorized version of the kepler angle

def KeplerAngle(e,M):
    """
    Given M, find nu: (Using newton's method). Operates on an array of M's
    """
    #Pick initial guess for E using Curtis's Guidlines    
    E = M - e/2;
    nu = np.zeros(len(M))
    for i,meanAnomoly in enumerate(M):

        if( meanAnomoly < np.pi):
            E[i] = meanAnomoly + e/2
        #Set Eold negative to garuntee we enter loop
        Eold = -1
        #Set the convergence bound requested
        j = 0
        while(np.abs(E[i]-Eold) > 10**-15):
            Eold = E[i]
            E[i] = E[i]- (E[i] - e * np.sin(E[i]) - meanAnomoly)/(1 - e * np.cos(E[i]))
            j += 1
            if (j>1000):
                break #Catches a possible inability to converge to desired tolerances
        nu[i] = 2 * np.arctan(np.sqrt((1+e)/(1-e))*np.tan(E[i]/2))

    return nu