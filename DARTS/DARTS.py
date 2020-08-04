import sys 
import os
sys.path.append(os.path.abspath("C:\\Users\\trey_\\Documents\\GitHub\\PROpy"))

import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.spatial.transform import Rotation
import pro_lib
import numpy as np
import math
#For visualization
import pandas as pd
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import matplotlib.animation as animation



def main():
    #Run on nominal NISAR mission
    #Sun Sync periodic, N = 173, D = 12
    #alt = 747, i = 98.4, e = 0 (assumed), other params unknown
    #Swath width: 240 km
    
    #Assume a nominal number of 6 s/c, initialized in a line .25 km apart cross track, .05 km apart along track
    
    
    # assigning parameters 
    
    NoRev = 173
    # orbital parameters

    altitude = 747
    ecc = 0
    inc = 98.4
    Om = 0
    om = 0
    f = 0

    num_deputy = 5


    mu = 398600.4418  #gravitational constant
    r_e = 6378.1363   # Earth Radius 
    J2 = 1082.64*10**(-6) #J2 Constant

    #Initial positions of DEPUTIES (Chief at 0,0,0)
    init_deputy = np.array([[-.5,-.1,0],
                            [-.25,-.05,0],
                            [.25,.05,0],
                            [.5,.1,0],
                            [.75,.15,0]])

    print(init_deputy)

    # compute initial conditions for the chief and deputy
    ys = pro_lib.initial_conditions_deputy("nonlinear_correction_linearized_j2_invariant", 
                                           [NoRev,altitude,ecc,inc,Om,om,f,5], init_deputy, mu,r_e,J2)
    print(ys)
    #Compute orbit duration
    '''
    fpo = input_info[8] # time steps per orbit

    # # Energy Matched 
    # # J2 disturbance 
    # # No drag

    # k_J2 = (3/2)*J2*mu*(r_e**2)

    # Orbital Elements
    a = r_e + altitude           # semimajor axis [km] (Assumed circular orbit)
    inc = inc*np.pi/180             # inclination [rad]


    # Xu Wang Parameters
    h = np.sqrt(a*(1 - ecc**2)*mu)           # angular momentum [km**2/s]


     # effective period of the orbit
    a_bar = a*(1 + 3*J2*r_e**2*(1-ecc**2)**0.5/(4*h**4/mu**2)*(3*np.cos(inc)**2 - 1))**(-2/3)  
    period = 2*np.pi/np.sqrt(mu)*a_bar**(3/2)  
 

    #  # simulation time
    time = np.linspace(0,NoRev*period,int(period/(fpo)))  
    # print(time)
    # orbit_num = time/period               # time vector with units of orbits instead of seconds
                                          # orbit = period of non J2 orbit
    '''

    #We know the orbit is periodic over 12 days, so just compute over those 12 days, every 10 seconds
    time = np.arange(0,12*24*3600,60)
    print("start")
    
 
    
    orbitState  = odeint(pro_lib.dyn_chief_deputies,ys,time,args=(mu,r_e,J2,num_deputy))
    print(np.shape(orbitState))
    print(orbitState)
    print(orbitState[0])
    #Plot the computed dynamics
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for i in range(num_deputy):
    
        ax.plot(orbitState[:,6*(i+1)],orbitState[:,6*(i+1)+1],orbitState[:,6*(i+1)+2])
    plt.draw()
    
    #Try feeding into our cost function
    print(costFunction(time,orbitState,30))
  
    
    #fig = go.Figure()
    ## Add traces, one for each slider step
    #for state in orbitState:
    #    xs = state[6::6]
    #    ys = state[7::6]
    #    zs = state[8::6]
    #    fig.add_trace(
    #        go.Scatter3d(
    #            visible=False,
    #            x=xs,y=ys,z=zs),
    #            range_x=[-1,1], range_y=[-1,1], range_z=[-1,1])
    #
    ## Make 0th trace visible
    #fig.data[0].visible = True
    #
    ## Create and add slider
    #steps = []
    #for i in range(len(fig.data)):
    #    step = dict(
    #        method="update",
    #        args=[{"visible": [False] * len(fig.data)}],
    #    )
    #    step["args"][0]["visible"][i] = True  # Toggle i'th trace to "visible"
    #    steps.append(step)
    #
    #sliders = [dict(
    #    active=0,
    #    currentvalue={"prefix": "Frequency: "},
    #    #pad={"t": 50},
    #    steps=steps
    #)]
    #
    #fig.update_layout(
    #    sliders=sliders
    #)
    #
    #
    #fig.show()


    #orbitData = np.array([orbitState[:,6::6],orbitState[:,7::6],orbitState[:,8::6]])
    #Writer = animation.writers['ffmpeg']
    #writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)
    #
    #fig = plt.figure(figsize=(10,10))
    #ax = fig.add_subplot(111, projection='3d')
    #ax.set_xlim(-1, 1)
    #ax.set_ylim(-1, 1)
    #ax.set_zlim(-1, 1)
    #print(len(time))
    #ani = animation.FuncAnimation(fig, animate, frames=len(time), fargs=(orbitData,ax))
    #
    #ani.save("demo2.mp4", writer=writer)


def animate(i,orbitData,ax):
    print(i)
    ax.clear()
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_zlim(-1, 1)
    data = orbitData[:,int(i):int(i+1),:] #select data range
    for i in range(len(data[0,0,:])):
        ax.plot(data[0,:,i],data[1,:,i],data[2,:,i],"o")

def costFunction(t,stateVector,lookAngle):
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
    scienceMerit = computeScienceMerit(t,stateVector,lookAngle)

    return (orbitCost - scienceMerit)

 
def computeOrbitCost(t,stateVector):
    return 0

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
        (groundTracks[i],r0s[i]) = groundAngleTrackComputation(stateVector[i,0:6],t[i],lookAngle)

    
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
        row = int(np.min((719,np.max((0,np.round((-lat+89.75)*4))))))
        col = int(np.min((1439,np.max((0,np.round((lon+180)*4))))))
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
        print(i)
        print(baselines[i])
        print (seperations[i])
        resolutions[i] = lam * r0s[i] / (2 * baselines[i])
        ambiguities[i] = lam * r0s[i] / (2 * seperations[i])
        #Note if violate baseline/ambiguity constraint, and how often

        if resolutions[i] < vegH[i]/5:
            numViolateRes +=1
        if ambiguities[i] < vegH[i]:
            numViolateAmb +=1


    print(resolutions)

    print("Violations:")
    print(numViolateRes)
    print(numViolateAmb)
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

    yhat = vhat - np.dot(rhat,vhat)*vhat #y direction is along velocity perpendicular to radial direction
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
    longitude = np.arctan2(targetVector[1],targetVector[0])
    latitude = np.arctan2(targetVector[2],np.sqrt(targetVector[1]**2+targetVector[0]**2))

    #And account for the Earth's rotation by subtracting off Earth's rotation
    longitude = longitude - t0* 4.178074346064814**(-3) *np.pi/180;

    #Compute r0, distance from chief to the target, also using the law of sines
    r0 = rotationAngle*(6378.1363 * 1000)/lookAngle
    
    #stack up the output
    groundTrack = np.vstack((latitude,longitude))
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
    #First space craft is the chief, which has position 0,0,0 by definition (the state vector data is the orbit elements, which we don't need)
    positions[0] = np.zeros(3)
    
    #Compute the look angle unit vector
    lookVector = np.array([-np.cos(np.radians(lookAngle)),0,np.sin(np.radians(lookAngle))])


    #Compute candidate baselines (since we don't know a priori which two s/c are farthest apart)
    #TODO: Can this be more efficient? Currently O(n^2), but small numbers (There are n chose 2 combos)
    #bestPair = np.zeros(2)
    bestLn = 0
    #Loop over combinations
    print("Positions:")
    print(positions)
    print("candiadateBaseline")
    for i in range(len(positions)-1):
        for j in range(len(positions) - 1 - i):
            candiadateBaseline = positions[i] - positions[j+i+1]
            print(candiadateBaseline)
            #Project onto the look angle and subtract this off to get component perpendicular to the look angle
            candiadateLn = np.linalg.norm(candiadateBaseline - np.dot(candiadateBaseline,lookVector)*candiadateBaseline/np.linalg.norm(candiadateBaseline))
            if candiadateLn > bestLn:
                bestLn = candiadateLn
    print("best")
    print(bestLn)
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
        positions[i] =  positions[i] - np.dot(positions[i],lookVector)*positions[i]/np.linalg.norm(positions[i])
        #Remove the along track (y) component
        projectedPositions[i] = np.array([positions[i,0],positions[i,2]])

    #Now we know they are all along a straight line, so we can sort by our x axis!
    sortIndex = np.argsort(projectedPositions[:,0])
    sortedPositions = projectedPositions[sortIndex]

    #Now loop through looking for max separation
    mu = 0
    for i in range(len(sortedPositions)-1):
        sep = np.linalg.norm(sortedPositions[i+1]-sortedPositions[i])
        if sep > mu:
            mu = sep

    return mu




if __name__ == "__main__":
    main()

   