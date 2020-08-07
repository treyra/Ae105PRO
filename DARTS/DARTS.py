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
from mayavi import mlab

#TODO: Make Earth radius / J2/ etc class constants?

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
    #init_deputy = np.array([[-.5,-.1,0],
    #                        [-.25,-.05,0],
    #                        [.25,.05,0],
    #                        [.5,.1,0],
    #                        [.75,.15,0]])

    init_deputy = np.array([[0,-.5,-.5,],
                            [0,-.25,-.25],
                            [0,.25,.25],
                            [0,.5,.5],
                            [0,.75,.75]])


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
    #Plot the computed dynamics
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title("6 space craft formation in NISAR J2 dynamic orbit, LVLH frame")
    ax.set_xlabel("x, radial out from Earth (km)")
    ax.set_ylabel("y, along track (km)")
    ax.set_zlabel("z, cross track (km)")
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_zlim(-1, 1)
    ax.azim = -100
    ax.elev = 43
    for i in range(num_deputy):
    
        ax.plot(orbitState[:,6*(i+1)],orbitState[:,6*(i+1)+1],orbitState[:,6*(i+1)+2])
    
    ax.plot([0],[0],[0],"ko")
    plt.draw()
    
    #Try feeding into our cost function
    print("Cost:")
    print(costFunction(time,orbitState,30))
    plt.show()
    
    print('ax.azim {}'.format(ax.azim))
    print('ax.elev {}'.format(ax.elev))
    #Save the user selected "best" veiw for animation
    azimuth = ax.azim
    elevation = ax.elev
    
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
    #ax.azim = azimuth
    #ax.elev = elevation
    #print(len(time))
    #ani = animation.FuncAnimation(fig, animate, frames=len(time), fargs=(orbitData,ax))
    #
    #ani.save("demo4.mp4", writer=writer)

#Methods for performing a genetic algorithm on various swarm setups
def optimize():
    for iterations in range(10):
        pass


def mutate(state,numOffspring,stdDeviation=.05):
    """
    Generates new random initial swarm configurations
    given a state to mutate from. The new states will
    be arranged in Gaussian fashion around the given 
    state with the specified standard deviation.

    Parameters
    ----------
    state : array, shape(3*numDeputies)
        initial deputy spatial configurations of the
        swarm. Should be (x,y,z) of each deputy in order
    numOffspring : int
        number of new positions to generate about this one
    stdDeviation : double (default = .05)
        standard Deviation of the normal distribution used
        to generate new arrangements

    Returns
    ---------
    offspring : array shape(numOffspring,3)
        Mutated swarm initial positions
    """

    #Create the output array
    offspring = np.zeros((numOffspring,3))

    #





def animate(i,orbitData,ax):
    print(i)
    ax.clear()
    ax.set_title("6 space craft formation in NISAR J2 dynamic orbit, LVLH frame")
    ax.set_xlabel("x, radial out from Earth (km)")
    ax.set_ylabel("y, along track (km)")
    ax.set_zlabel("z, cross track (km)")
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_zlim(-1, 1)
    data = orbitData[:,int(i):int(i+1),:] #select data range
    for i in range(len(data[0,0,:])):
        ax.plot(data[0,:,i],data[1,:,i],data[2,:,i],"o")
    
    ax.plot([0],[0],[0],"ko")

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
    print(orbitCost)

    #Compute the science merit that offsets the cost
    scienceMerit = computeScienceMerit(t,stateVector,lookAngle)

    return (orbitCost - scienceMerit)

 
def computeOrbitCost(t,stateVector):
    """
    Return a cost for the formation based on the initial delta-V
    needed to assume the formation. This is assumed to be the initial
    relative velocity of each deputy, assuming the initial drift to the
    starting positions is negligible 
    
    Parameters
    ----------
    t : array
        Time at which each state is provided (unused, left for future costs
        that use t)
    stateVector : array, shape (len(t), len(state))
        State throughout the orbit at each time step. State should be 
        (x,y,z,vx,vy,vz) stacked for each space craft. First 6 elements
        should be the chief orbit in Xu Wang Parameters
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
    numSC = int(len(stateVector[0])/6)

    #Compute delta-V to initialize orbit
    deltaV = 0
    for i in range(numSC - 1):
        initialV = stateVector[0,3+6*(i+1):6+6*(i+1)] * 1000 #Want in m/s
        print(initialV)
        deltaV += np.linalg.norm(initialV)

    #Exponentiate to give cost
    return np.exp(deltaV)

#TODO: Consider consolidating to clean up
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
        (x,y,z,vx,vy,vz) stacked for each space craft. First 6 elements
        should be the chief orbit in Xu Wang Parameters
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
    numSC = int(len(stateVector[0])/6) #UNUSED
    
    #Define targets/ s/c position to see them

    #Load the elevation data in from https://webmap.ornl.gov/ogc/dataset.jsp?ds_id=10023 at .25 degree resolution
    #If finer resolution is desired the hard coding will need to be changed
    elevationData = np.genfromtxt('elevationData.txt',delimiter=' ',skip_header=5) #Units, meters
    #bottom left corner corresponds to -180 W -90 S in this projection, so upper right is 179.75E 89.75N

    #Compute the ground track of the center of the swath of the target location over the orbit
    #Time,ground track
    groundTracks = np.zeros((len(t),2))
    #Distances to the target
    r0s = np.zeros(len(t))
    for i in range(len(groundTracks)):
        #translating back to state vecotr from the orbital elements
        chiefState = stateVector[i,0] * np.dot(pro_lib.rotation_matrix_lvlh_to_eci(stateVector[i,3],stateVector[i,5],stateVector[i,4]),[1,0,0])
        #Compute along track unit vector
        yhat = np.dot(pro_lib.rotation_matrix_lvlh_to_eci(stateVector[i,3],stateVector[i,5],stateVector[i,4]),[0,1,0])
        (groundTracks[i],r0s[i]) = groundAngleTrackComputation(chiefState,yhat,t[i],lookAngle)

    
    #Function for visualizing the ground tracks
    visualize(np.radians(groundTracks[:,0]),np.radians(groundTracks[:,1]),elevationData)

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
        ##Check if we're over a target!
        #if vegH[i]  > 0:
        #    #Compute the baseline
        #    baselines[i] = computeBaseline(stateVector[i])
        #    print(baselines)
        #    #Compute the ambiguity
        #    seperations[i] = computeMaxSeperation(stateVector[i])
        #    print(seperations)
        
        #Compute the baseline
        baselines[i] = computeBaseline(stateVector[i],lookAngle)
        #Compute the ambiguity
        seperations[i] = computeMaxSeperation(stateVector[i],lookAngle)
    
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
    #Currently giving score as # of times better than constraint (1/5 of veg height)
    score = 0
    for i in range(len(t)):
        resolutions[i] = lam * r0s[i] / (2 * baselines[i]*1000) #Convert to meters
        ambiguities[i] = lam * r0s[i] / (2 * seperations[i])
        #Check if over a  target!
        #Note if violate baseline/ambiguity constraint, and how often
        if vegH[i]  > 0:
            if resolutions[i] > vegH[i]/5:
                numViolateRes +=1
            if ambiguities[i] < vegH[i]:
                numViolateAmb +=1
            score += (vegH[i]/resolutions[i])/5
    print(len(vegH))
    print(len(np.where(vegH > 0)[0]))
    print(resolutions)

    print("Violations:")
    print(numViolateRes)
    print(numViolateAmb)
    return score

    
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
def groundAngleTrackComputation(x,yhat,t0,lookAngle):
    """
    Return a latitude/longitude trajectory along the Earth of the look target
    given the position vector along the orbit and the starting time difference 
    from 0:00 UT and the look angle
    
    Parameters
    ----------
    x : array, shape (3)
        Position vector of the space craft at the time step, as 
        (x,y,z) in the ECI frame
    yhat : array, shape (3)
        Along track unit vector in the ECI frame
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
        Earth's surface in degrees
    r0 : double
        The distance in km from the space craft to the target position
        on the surface of the Earth.
    """

    #Convert look angle to rads
    rlookAngle = np.radians(lookAngle)

    #Compute the radial and unit vector
    rhat = x/np.linalg.norm(x) #Also xhat, as this is the x direction in LVLH



    #Compute the view unit vector by rotating the -rhat vector around the yhat vector by the look angle
    rotate = Rotation.from_rotvec(rlookAngle * yhat)
    lookVector = rotate.apply(-rhat) #UNUSED

    #Compute the angle we need to rotate rhat by to get look target using law of sines
    #(we know the look angle and its opposite side (radius of the Earth)
    wideAng = np.arcsin(np.linalg.norm(x)* np.sin(np.abs(rlookAngle))/(6378.1363))
    rotationAngle =  np.pi - (wideAng + rlookAngle)
    #Compute r0, distance from chief to the target, also using the law of sines
    r0 = np.sin(rotationAngle)*(6378.1363)/np.sin(rlookAngle)

    #Check if physical or we got the wrong quadrant
    if r0 > np.linalg.norm(x):
        #wideAng was in wrong quadrant so should have been wideAng = np.pi - wideAng
        rotationAngle = wideAng - rlookAngle
        r0 = np.sin(rotationAngle)*(6378.1363)/np.sin(rlookAngle)
    
    #Rotate rhat the opposite way we rotated the lookVector
    rotate2 = Rotation.from_rotvec(-np.sign(lookAngle)*rotationAngle * yhat)
    targetVector = rotate2.apply(rhat)

    #Convert into lon/lat using the fact that we are on a sphere
    longitude = np.arctan2(targetVector[1],targetVector[0])
    latitude = np.degrees(np.arctan2(targetVector[2],np.sqrt(targetVector[1]**2+targetVector[0]**2)))

    #And account for the Earth's rotation by subtracting off Earth's rotation
    longitude = np.degrees(longitude - t0* 4.178074346064814**(-3) *np.pi/180);
    longitude = np.mod(longitude,360)              
    
    #stack up the output
    groundTrack = np.vstack((latitude,longitude))
    return (groundTrack.transpose(),r0)

def visualize(latitude,longitude,elevationData):
    """
    Plots a 3D visualization of the target ground track

    Parameters
    ----------
    latitude : array shape (len(t))
        list of latitudes at each time to visualize
    longitude : array shape (len(t))
        list of longitudes at each time to visualize
    """

    earthRadius = 6378.1363

    ##Plot Ground track
    #fig = plt.figure()
    #ax = fig.gca(projection='3d')
    #
    #
    #groundTrack = np.array([np.cos(longitude)*np.cos(latitude),np.sin(longitude)*np.cos(latitude),np.sin(latitude)])*earthRadius
    #print(np.shape(groundTrack))
    #print(groundTrack[:,0])
    #ax.plot(groundTrack[0],groundTrack[1],groundTrack[2])
    #
    ##Plot a sphere
    #
    #u, v = np.mgrid[0:2*np.pi:100j, 0:np.pi:50j]
    #x = np.cos(u)*np.sin(v)*earthRadius*.9
    #y = np.sin(u)*np.sin(v)*earthRadius*.9
    #z = np.cos(v)*earthRadius*.9
    #
    #
    #surf = ax.plot_surface(x,y,z,cmap=plt.cm.Spectral)
    #surf.set_alpha(.5)

    #Plot the elevation data



    #plt.show()
    

    #Plot Ground track
    #fig = plt.figure()
    #ax = fig.gca(projection='3d')
    
    
    groundTrack = np.array([np.cos(longitude)*np.cos(latitude),np.sin(longitude)*np.cos(latitude),np.sin(latitude)])*earthRadius
    #print(np.shape(groundTrack))
    #print(groundTrack[:,0])
    #ax.plot(groundTrack[0],groundTrack[1],groundTrack[2])
    
    #Plot a sphere
    
    u, v = np.mgrid[0:2*np.pi:100j, 0:np.pi:50j]
    x = np.cos(u)*np.sin(v)*earthRadius
    y = np.sin(u)*np.sin(v)*earthRadius
    z = np.cos(v)*earthRadius
    
    
    #s = mlab.mesh(x, y, z)
    mlab.plot3d(groundTrack[0],groundTrack[1],groundTrack[2],line_width = 100)
    mlab.show()

    #Plot the elevation data
    pass



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
    for i in range(len(positions)-1):
        for j in range(len(positions) - 1 - i):
            candiadateBaseline = positions[i] - positions[j+i+1]
            #Project onto the look angle and subtract this off to get component perpendicular to the look angle
            candiadateLn = np.linalg.norm(candiadateBaseline - np.dot(candiadateBaseline,lookVector)*candiadateBaseline/np.linalg.norm(candiadateBaseline))
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
    
    #First space craft is the chief, which has position 0,0,0 by definition (the state vector data is the orbit elements, which we don't need)
    positions[0] = np.zeros(3)

    #Compute the look angle unit vector
    lookVector = np.array([-np.cos(np.radians(lookAngle)),0,np.sin(np.radians(lookAngle))])
    
    #Project onto the look angle and subtract this off to get component perpendicular to the look angle
    #Want coordinates in the plane centered at the origin of the LVLH system, perpendicular to the look
    #angle. Will then remove the along track dimension and get a 1D arrangement and sort them
    projectedPositions = np.zeros((int(len(stateVector)/6),2))
    #Only loop through deputies, as chief will be at [0,0,0] by definition
    for i in range(len(positions)-1):
        positions[i+1] =  positions[i+1] - np.dot(positions[i+1],lookVector)*positions[i+1]/np.linalg.norm(positions[i+1])
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

if __name__ == "__main__":
    main()

   