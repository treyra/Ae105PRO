import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import InterpolatedUnivariateSpline

def splineOrbitInterplation(t_meas,x_meas,t_interp):
    """
    Given a set of orbit state estimates, simply interpolates between them using a cubic interpolation method. 
    This method takes advantage of the fact that orbits are known to be smooth functions to give an accurate estimate.
    :param t_meas: Time of estimates
    :param x_meas: State vector (position and velocity in a nx6 dimensional vector with each state having format [rx,ry,rz,vx,vy,vz])
     estimates to  interpolate between
    :param t_interp: Times to provide interpolation at
    :return x_interp: Interpolated state vector at each t_interp
    """

    #Make the output array
    x_interp = np.zeros((len(t_interp),6))

    #Make a spline interpretation for each axis
    for i in range(6):

        interpFun = interp1d(t_meas,x_meas[:,i], kind='cubic')
        x_interp[:,i] =  interpFun(t_interp)
    #Return interpolated states
    return x_interp


def splineOrbitInterplationPos(t_meas,x_meas,t_interp):
    """
    Given a set of orbit state estimates, simply interpolates between them using a cubic interpolation method. 
    This method takes advantage of the fact that orbits are known to be smooth functions to give an accurate estimate.
    :param t_meas: Time of estimates
    :param x_meas: State vector (position in a nx3 dimensional vector with each state having format [rx,ry,rz])
     estimates to  interpolate between
    :param t_interp: Times to provide interpolation at
    :return x_interp: Interpolated state vector at each t_interp
    """

    #Make the output array
    x_interp = np.zeros((len(t_interp),3))

    #Make a spline interpretation for each axis
    for i in range(3):

        interpFun = interp1d(t_meas,x_meas[:,i], kind='cubic')
        x_interp[:,i] =  interpFun(t_interp)
    #Return interpolated states
    return x_interp


def splineOrbitInterplationPosMulti(t_meas,x_meas,t_interp):
    """
    Given a set of orbit state estimates for NumSC of space craft, simply interpolates between them using a cubic interpolation method. 
    This method takes advantage of the fact that orbits are known to be smooth functions to give an accurate estimate.
    :param t_meas: Time of estimates
    :param x_meas: State vector (position in a nxNumSCx3 dimensional vector with each state having format [rx,ry,rz])
     estimates to  interpolate between
    :param t_interp: Times to provide interpolation at
    :return x_interp: Interpolated state vector at each t_interp
    """

    #Make the output array
    x_interp = np.zeros((len(t_interp),len(x_meas[0,:,1]),3))
    for i in range(len(x_meas[0,:,1])):
        x_interp[:,i,:] = splineOrbitInterplationPos(t_meas,x_meas[:,i,:],t_interp)

    




def splineOrbitInterplationPosSpline(t_meas,x_meas):
    """
    Given a set of orbit state estimates, simply interpolates between them using a cubic interpolation method. 
    This method takes advantage of the fact that orbits are known to be smooth functions to give an accurate estimate.
    Note that the velocity data in the input is not used, but is accepted to provide compatibility with the state vector produced by simulations

    :param t_meas : Time of estimates
    :param x_meas : Position vector (position (and velocity) in a nx3(+3) dimensional vector with each state having format [rx,ry,rz(,vx,vy,vz)]). Velocity components unused
     estimates to  interpolate between

    :return posSplines : list, shape(numDeputies,3) Spline objects in each axis (x,y,z) that can interpolate the position vector at each t_interp
    """

    #Output list
    posSplines = []

    #Make a spline interpretation for each axis
    for i in range(3):

        posSplines.append(InterpolatedUnivariateSpline(t_meas,x_meas[:,i],k=3))
      
    #Return splines
    return posSplines






