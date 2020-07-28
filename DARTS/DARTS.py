import sys 
import os
sys.path.append(os.path.abspath("C:\\Users\\trey_\\Documents\\GitHub\\PROpy"))

import matplotlib.pyplot as plt
import scipy
import pro_lib
import numpy as np
import math

import scipy.spatial.transform.Rotation


def main():
    #computeScienceMerit(0,0)
    #return
    
    # assigning parameters 
    """
    NoRev = input[0] # number of orbits

    # orbital parameters

    altitude = input[1]
    ecc = input[2]
    INC = input[3]
    Om = input[4]
    om = input[5]
    f = input[6]


    deputy_num = input[7] 

    fpo = input[8] # time steps per orbit
    """

    print(math.factorial(3))
    return

    mu = 398600.4418  #gravitational constant
    r_e = 6378.1363   # Earth Radius 
    J2 = 1082.64*10**(-6) #J2 Constant
    # info
    input_info = np.zeros([9]) 
    input_info[0] = 1  #no of orbits
    input_info[1] = 500 #altitude
    input_info[2] = 0.0 #ecc
    input_info[3] = 60 #INC
    input_info[4] = 0.0 #Om
    input_info[5] = 60 #om
    input_info[6] = 20 #f
    input_info[7] = 3 #num_deputies
    input_info[8] = 100 #output per orbit

    # intial position on deputies
    
    num_deputy = int(input_info[7])

    init_deputy = np.zeros([num_deputy,3])


    for i in range(num_deputy):
        #UNITS ??
        init_deputy[i,0] = 50*np.sin(2*np.pi*i/num_deputy)
        init_deputy[i,1] = 50*np.cos(2*np.pi*i/num_deputy) 
        init_deputy[i,2] = 10 + (-1)**(i)

    dyn = main2(input_info,init_deputy,mu,r_e,J2)

    print("Output") 
    print(np.shape(dyn))
    print(dyn)


def costFunction(t,stateVector):
    """
    Return a cost value for the given orbit trajectory. Intended to compare 
    for a single orbit. TODO: Lifetime better?
    
    Parameters
    ----------
    t : array
        Time at which each state is provided
    stateVector : array, shape (len(t), len(state)) TODO: DETERMINE IF WE HAVE ANY USE FOR ABS STATES
        State throughout the orbit at each time step. State should be 
        (x,y,z,vx,vy,vz) stacked for each space craft
        First s/c is the chief and should be in ECI coordinates
        Each other s/c is a deputy and should be in RELATIVE
        coordinates TODO: Error check?

    Returns
    -------
    J : double
        Numerical cost of the orbit. A lower cost means a more 
        efficient/desirable orbit.
    
    """

    #Cost function we are trying to evaluate
    #J = 


    #Compute the cost to set up the orbit, and maintain it
    orbitCost = computeOrbitCost(t,stateVector)

    #Compute the science merit that offsets the cost
    scienceMerit = computeScienceMerit(t,stateVector)

 


def computeScienceMerit(t,stateVector,lookAngle=0):
    """
    Return a science merit score for the given orbit trajectory. 
    Intended to be over a single orbit 

    Assumptions: Earth is modeled as a perfect sphere for this 
    approximation with radius 6378.1363 km
    
    Parameters
    ----------
    t : array
        Time at which each state is provided
    stateVector : array, shape (len(t), len(state))
        State throughout the orbit at each time step. State should be 
        (x,y,z,vx,vy,vz) stacked for each space craft
    lookAngle : double (default 0)
        Angle the target direction makes with the nadir direction, 
        positive defined as a right handed rotation about the along
        track direction, 0 is nadir.
    Returns
    -------
    J : double
        Numerical cost of the orbit. A lower cost means a more 
        efficient/desirable orbit.
    
    """

    #Compute the number of space craft in the formation
    numSC = int(len(stateVector[0]/6))
    
    #Define targets/ s/c position to see them

    #Load the elevation data in from https://webmap.ornl.gov/ogc/dataset.jsp?ds_id=10023 at .25 degree resolution
    #If finer resolution is desired the hard coding will need to be changed
    elevationData = np.genfromtxt('elevationData.txt',delimiter=' ',skip_header=5) #Units, meters
    #bottom left corner corresponds to -180 W -90 S in this projection, so upper right is 179.75E 89.75N

    #Compute the ground track of the center of the swath of the target location over the orbit

#Old method, unsure if we need more than chief ground track
#    #Time,s/c,ground track
#    groundTracks = np.zeros((len(t),numSC,2))
#    
#    #Default ground looking approach, using the s/c ground track
#    if lookAngle == 0:
#        #Compute ground track over orbit, including the effect of earth rotating
#        for i in range(numSC):
#            groundTracks = groundTrackComputation(stateVector[:,0+i*6:3+i*6])
#    #Non zero look angle (TODO: is the look angle constant?)
#    #TODO: Likely redundant with angle=0 computation
#    else:
#        #Compute ground track over orbit, including the effect of earth rotating
#        for i in range(numSC):
#            groundTracks = groundAngleTrackComputation(stateVector[:,0+i*6:6+i*6],lookAngle)
    
    #Time,ground track
    groundTracks = np.zeros((len(t),2))
    #Distances to the target
    r0s = np.zeros(len(t))
    for i in range(len(groundTracks)):
        (groundTracks[i],r0s[i]) = groundAngleTrackComputation(stateVector[i,0:6],lookAngle)

    
    #When over target, compute baseline, ambiguity
    baselines = np.zeros(len(t))
    seperations = np.zeros(len(t))
    vegH = np.zeros(len(t))

    #Use the current target ground position to determine if we are over a target
    #TODO: Incorporate the look angle functionality (add offset to ground track computation? or to the positioning function?)
    for i,position in enumerate(groundTracks):
        #Find lat/lon of the chief
        lat = position[0]
        lon = position[1]

        #Compare with the vegetation data to see if we are over a target
        #Will compute where to look in elevationData matrix
        #Bound by extremes of array to avoid falling of the end of the array
        row = int(np.min(719,np.max(0,np.round((-lat+89.75)*4)))) 
        col = int(np.min(1439,np.max(0,np.round((lon+180)*4))))
        vegH[i] = elevationData[row,col]
        #Check if we're over a target!
        if vegH[i]  > 0:
            #Compute the baseline
            baselines[i] = computeBaseline(stateVector[i])
            #Compute the ambiguity
            seperations[i] = computeMaxSeperation(stateVector[i])

    #Now loop through and see how often we violate our constraints:
    #   Resolution > vegH/5
    #   ambiguity > vegH
    #
    #First need some data
    #Resolution due to baseline:
    #delta_n = lambda * r0 / (p_delta*Baseline)
    #
    #Nearest ambiguity location
    #h_n = k lambda * r0 / (p_a * mu)
    #
    #Where the variables are:
    #lambda: wavelength, using NISAR S-band as a reference (9 cm)
    #r0 = H/cos(theta), distance to target where H is the orbit altitude
    #p_delta,p_a = 2 for SAR
    #k = 1,2,3 (kth ambiguity, we'll take k = 1)
    
    lam = .09 #Units in meters
    resolutions = np.zeros(len(t))
    ambiguities = np.zeros(len(t))
    numViolateRes = 0
    numViolateAmb = 0
    for i in range(len(t)):
        resolutions[i] = lam * r0s[i] / (2 * baselines[i])
        ambiguities[i] = lam * r0s[i] / (2 * seperations[i])
        #Note if violate baseline/ambiguity constraint, and how often

        if resolutions[i] < vegH[i]/5:
            numViolateRes +=1
        if ambiguities[i] < vegH[i]:
            numViolateAmb +=1



    #Score
    return np.sum(resolutions)

    
    #TODO, account for starting position of the Earth
    #TODO: Verify that the origin lines up with 0 W 0 N
def groundTrackComputation(r,t0):
    """
    Return a latitude/longitude trajectory along the Earth given the
    position vector along the orbit and the starting time difference from
    0:00 UT
    
    Parameters
    ----------
    r : array, shape (number of positions, 3)
        Position vector of the space craft at each step, as (x,yz) in 
        the ECI frame
    t0 : double
        Time since 0:00 UT 

    Returns
    -------
    groundTrack : array, shape (len(r), 2)
        The latitude and longitude of the trajectory along the Earth's 
        surface
    
    """

    #Convert into lon/lat using the fact that we are on a sphere
    longitude = np.arctan2(r[:,2],r[:,1])
    latitude = np.arctan2(r[:,3],np.sqrt(r[:,2]**2+r[:,1]**2))

    #And account for the Earth's rotation by subtracting off Earth's rotation
    longitude = longitude - t0* EARTH.rotD *pi/180;
    
    #stack up the output
    groundTrack = np.vstack(latitude,longitude)
    return groundTrack.transpose()


    #TODO, account for starting position of the Earth
    #TODO: Verify that the origin lines up with 0 W 0 N
def groundAngleTrackComputation(x,t0,lookAngle):
    """
    Return a latitude/longitude trajectory along the Earth of the look target
    given the position vector along the orbit and the starting time difference 
    from 0:00 UT and the look angle
    
    Parameters
    ----------
    x : array, shape (6)
        Position vector of the space craft at the time step, as 
        (x,y,z,vx,vy,vz) in the ECI frame
    t0 : double
        Time since 0:00 UT 
    lookAngle : double
        Angle the target direction makes with the nadir direction, 
        positive defined as a right handed rotation about the along
        track direction, 0 is nadir. (Degrees)

    Returns
    -------
    groundTrack : array, shape (len(r), 2)
        The latitude and longitude of the look trajectory along the 
        Earth's surface
    """

    #Compute the radial and along track unit vectors
    r = x[0:3]
    rhat = r/np.linalg.norm(r) #Also xhat, as this is the x direction in LVLH

    v = x[3:6]
    vhat = v/np.linalg.norm(v)

    yhat = vhat - np.linalg.dot(rhat,vhat)*vhat #y direction is along velocity perpendicular to radial direction
    yhat = yhat/np.linalg.norm(yhat)

    #Compute the view unit vector by rotating the -rhat vector around the yhat vector by the look angle
    rotate = Rotation.from_rotvec(np.radians(lookAngle) * yhat)
    lookVector = rotate.apply(-rhat) #UNUSED

    #Compute the angle we need to rotate rhat by to get look target using law of sines
    #(we know the look angle and its opposite side (radius of the Earth)
    rotationAngle =  180 - ( + np.linalg.norm(r)* np.abs(lookAngle)/(6378.1363 * 1000) + lookAngle)
    
    #Rotate rhat the opposite way we rotated the lookVector
    rotate2 = Rotation.from_rotvec(-np.sign(lookAngle)*np.radians(rotationAngle) * yhat)
    targetVector = rotate2.apply(rhat)

    #Convert into lon/lat using the fact that we are on a sphere
    longitude = np.arctan2(targetVector[:,2],targetVector[:,1])
    latitude = np.arctan2(targetVector[:,3],np.sqrt(targetVector[:,2]**2+targetVector[:,1]**2))

    #And account for the Earth's rotation by subtracting off Earth's rotation
    longitude = longitude - t* EARTH.rotD *pi/180;

    #Compute r0, distance from chief to the target, also using the law of sines
    r0 = rotationAngle*(6378.1363 * 1000)/lookAngle
    
    #stack up the output
    groundTrack = np.vstack(latitude,longitude)
    return (groundTrack.transpose(),r0)

#TODO: Implement projection onto look angle
def computeBaseline(stateVector, lookAngle=0):
    """
    Return the baseline for the formation at the current time of the orbit.
    
    Parameters
    ----------
    stateVector : array
        State Vector of the swarm at the time to compute the baseline for.
        State should be (x,y,z,vx,vy,vz) stacked for each space craft. 
        x,y,z directions defined by the LVLH coordinate system, where
        x is radially out from the Earth to the s/c, y is along track and
        z is the orbit angular momentum.
    lookAngle : double (default is 0)
        Angle between the nadir and target directions in degrees. Alternatively 
        the angle between the -x direction and the target direction. Note that
        the target is always in the cross track plane, or the plane normal
        to the orbit trajectory. Assumed a positive look angle is a positive
        rotation about the along track vector y.

    Returns
    -------
    Ln : double
        Baseline for the formation, the spatial extent of the spacecraft
        projected normal to the look angle
    """

    #Pull out the positions of each space craft
    positions = np.zeros((int(len(stateVector)/6),3))
    for i in range(len(positions)):
        positions[i] = stateVector[i*6:i*6+3]
    
    #Compute the look angle unit vector
    lookVector = np.array([-np.cos(np.radians(lookAngle)),0,np.sin(np.radians(lookAngle))])


    #Compute candidate baselines (since we don't know a priori which two s/c are farthest apart)
    #TODO: Can this be more efficient? Currently O(n^2), but small numbers (There are n chose 2 combos)
    #bestPair = np.zeros(2)
    bestLn = 0
    #Loop over combinations
    for i in range(len(position)-1):
        for j in range(len(position) - 1 - i):
            candiadateBaseline = positions[i] - position[j+i+1]
            #Project onto the look angle and subtract this off to get component perpendicular to the look angle
            candiadateLn = candiadateBaseline - np.linalg.dot(candiadateBaseline,lookVector)*candiadateBaseline/np.linalg.norm(candiadateBaseline)
            if candiadateLn > bestLn:
                bestLn = candiadateLn
    return bestLn

def computeMaxSeperation(stateVector, lookAngle=0):
    """
    Return the maximum cross track separation to estimate the ambiguity.
    
    Parameters
    ----------
    stateVector : array
        State Vector of the swarm at the time to compute the baseline for.
        State should be (x,y,z,vx,vy,vz) stacked for each space craft. 
        x,y,z directions defined by the LVLH coordinate system, where
        x is radially out from the Earth to the s/c, y is along track and
        z is the orbit angular momentum.
    lookAngle : double (default is 0)
        Angle between the nadir and target directions in degrees. Alternatively 
        the angle between the -x direction and the target direction. Note that
        the target is always in the cross track plane, or the plane normal
        to the orbit trajectory. Assumed a positive look angle is a positive
        rotation about the along track vector y.

    Returns
    -------
    mu : double
        Maximum gap between any two space craft perpendicular to the look angle
        in the cross track plane
    """
    
    #Pull out the positions of each space craft
    positions = np.zeros((int(len(stateVector)/6),3))
    for i in range(len(positions)):
        positions[i] = stateVector[i*6:i*6+3]
    
    #Compute the look angle unit vector
    lookVector = np.array([-np.cos(np.radians(lookAngle)),0,np.sin(np.radians(lookAngle))])
    
    #Project onto the look angle and subtract this off to get component perpendicular to the look angle
    #Want coordinates in the plane centered at the origin of the LVLH system, perpendicular to the look
    #angle. Will then remove the along track dimension and get a 1D arrangement and sort them
    projectedPositions = np.zeros((int(len(stateVector)/6),2))
    for i in range(len(positions)):
        positions[i] =  positions[i] - np.linalg.dot(positions[i],lookVector)*positions[i]/np.linalg.norm(positions[i])
        #Remove the along track (y) component
        projectedPositions[i] = np.array([positions[i,0],positions[i,2]])

    #Now we know they are all along a straight line, so we can sort by our x axis!
    sortIndex = np.argsort(projectedPositions[:,0])
    sortedPositions = projectedPositions[sortIndex]

    #Now loop through looking for max separation
    mu = 0
    for i in range(len(sortedPostions)-1):
        sep = np.linalg.norm(sortedPositions[i+1]-sortedPositions[i])
        if sep > mu:
            mu = sep

    return mu





def main2(input_info,initial_xyz,mu,r_e,J2):

    initial_condition_type = "nonlinear_correction_linearized_j2_invariant"

    # compute initial conditions for the chief and deputy
    ys = pro_lib.initial_conditions_deputy(initial_condition_type, input_info, initial_xyz, mu,r_e,J2)

    # # assigning parameters 

    NoRev = input_info[0] # number of orbits

    # # orbital parameters

    altitude = input_info[1]
    ecc = input_info[2]
    INC = input_info[3]
    Om = input_info[4]
    om = input_info[5]
    f = input_info[6]


    deputy_num = int(input_info[7]) 

    fpo = input_info[8] # time steps per orbit

    # # Energy Matched 

    # # J2 disturbance 

    # # No drag

    # ##----Constants--------------

    # mu = 398600.4418  #gravitational constant
    # r_e = 6378.1363   # Earth Radius 
    # J2 = 1082.64*10**(-6) #J2 Constant

    # k_J2 = (3/2)*J2*mu*(r_e**2)

    # Orbital Elements
    a = r_e + altitude           # semimajor axis [km]
    inc = INC*np.pi/180             # inclination [rad]


    # Xu Wang Parameters
    h = np.sqrt(a*(1 - ecc**2)*mu)           # angular momentum [km**2/s]

    # # Number of Satellites
    # sat_num = int(deputy_num + 1) 
    

     # effective period of the orbit
    a_bar = a*(1 + 3*J2*r_e**2*(1-ecc**2)**0.5/(4*h**4/mu**2)*(3*np.cos(inc)**2 - 1))**(-2/3)  
    period = 2*np.pi/np.sqrt(mu)*a_bar**(3/2)  
 

    #  # simulation time
    time = np.linspace(0,NoRev*period,int(period/(fpo)))  
    # print(time)
    # orbit_num = time/period               # time vector with units of orbits instead of seconds
                                          # orbit = period of non J2 orbit
  

    #
    # param for simulation

    param = np.array([mu, r_e, J2,deputy_num])
    
    # run the dynamics

    # T = len(time)
    
    # total_dyn = np.zeros([sat_num*6,T]) 
    # total_dyn[:,0] = ys

    # dt = time[1]
    # print(dt)
    
    sol  = scipy.integrate.odeint(pro_lib.dyn_chief_deputies,ys,time,args=(param[0],param[1],param[2],param[3]))
    print(sol[:,1])

    #for tt in range(T-1):
    #    total_dyn[:,tt+1] = pro_lib.EulerInt(pro_lib.Dyn_Chief_Deputies(total_dyn[:,tt],param),dt,total_dyn[:,tt])

    # result is table with rows = time
    # columns = orbital parameters (chief + non-chief)

    # fig = plt.figure()
    # ax = plt.axes(projection='3d')
    #------------------------------------------------------------------------------------------#
    #----------------------------------------------------- Chief Solution------_---------------#
    #------------------------------------------------------------------------------------------#
    # # chief
    # r = sol[:,0] # (geocentric distance)
    # v_x = sol[:,1] # (radial velocity)
    # h = sol[:,2] # (angular momentum)
    # Omega = sol[:,3] # (right ascension of the ascending node)
    # i = sol[:,4] #(orbit inclination)
    # theta = sol[:,5] # (argument of latitude)

    # # Convert to cartesian coordinates
    # # see http://farside.ph.utexas.edu/teaching/celestial/Celestialhtml/node34.html
    # X = r*(np.cos(Omega)*np.cos(theta) - np.sin(Omega)*np.sin(theta)*np.cos(i))
    # Y = r*(np.sin(Omega)*np.cos(theta) + np.cos(Omega)*np.sin(theta)*np.cos(i))
    # Z = r*(np.sin(i)*np.sin(theta))

    # # ax.plot3D(X,Y,Z)

    # # deputies are in relative coordinates
    # for k in range(1, 4):
    #     # ax.plot3D(X + sol[:,k*6+0],Y + sol[:,k*6+1],Z + sol[:,k*6+2])
    #     ax.plot3D(sol[:,k*6+0],sol[:,k*6+1],sol[:,k*6+2])
    # plt.show()

    return sol


if __name__ == "__main__":
    main()

   