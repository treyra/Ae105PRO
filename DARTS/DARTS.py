import sys 
import os
#For finding PROpy if not in the base directory (which it is in the default repo, developed seperatly though)
#sys.path.append(os.path.abspath("C:\\Users\\trey_\\Documents\\GitHub\\PROpy"))

import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.spatial.transform import Rotation
import pro_lib
import numpy as np
#For parallel evaluation
from multiprocessing import Pool
#For data export
import netCDF4 as nc


#Class constants (Earth parameters) can be overwritten locally or when calling the utility methods
mu = 398600.432896939164493230  #gravitational constant
r_e = 6378.136  # Earth Radius 
J2 = 0.001082627 #J2 Constant




def main():

    #Demo cost computation and visualization
    #demo()

    #Perform optimization
    #optimize()


    #Exporting


    bestInit = np.array([[-0.11474911, -0.28277769, -0.41374357],
                              [ 0.11759513,  0.35926752, -0.29420415],
                              [-0.4380676,   0.13452622,  0.16840377],
                              [-0.80578882,  0.30974955,  0.41924403],
                              [ 0.18555127,  0.55242093,  0.95920626]])
    #Time we integrate each orbit over
    time = np.arange(0,12*24*3600,1)
     
    # orbital parameters, wrapped in a dictionary

    orbParams = {"time":time, #time steps over which to evaluate
                "NoRev":173, # revolutions considered, orbit altitude, orbit 
                "altitude":747, #orbit altitude
                "ecc":0, #orbit eccentricity
                "inc":98.4, #orbit inclination (deg)
                "Om":0, #orbit right ascension of ascending node (deg)
                "om":0, #orbit argument of periapsis (deg)
                "f":0, #orbit true anomaly (deg)
                "num_deputy":5, 
                "lookAngle":30, #(deg)
                "mu":mu,  #gravitational parameter of Earth
                "r_e":r_e,  # Earth Radius 
                "J2":J2} #J2 Constant of Earth
    exportOrbit(bestInit,orbParams)




def demo():
    """
    Function that demonstrates the scoring of a
    sample trajectory. The orbit used is based off the
    NISAR mission
    """


    #Run on nominal NISAR mission
    #Sun Sync periodic, N = 173, D = 12 (slight errors, probably due to higher order terms than J2 or minor corrections to the publicly availible data)
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


  

    #Initial positions of DEPUTIES (Chief at 0,0,0)
    
    #Original arbitrary seed
    #init_deputy = np.array([[0,-.5,-.5,],
    #                        [0,-.25,-.25],
    #                        [0,.25,.25],
    #                        [0,.5,.5],
    #                        [0,.75,.75]])

    init_deputy = np.array([[-0.11474911, -0.28277769, -0.41374357],
                              [ 0.11759513,  0.35926752, -0.29420415],
                              [-0.4380676,   0.13452622,  0.16840377],
                              [-0.80578882,  0.30974955,  0.41924403],
                              [ 0.18555127,  0.55242093,  0.95920626]])


    # compute initial conditions for the chief and deputy
    ys = pro_lib.initial_conditions_deputy("nonlinear_correction_linearized_j2_invariant", 
                                           [NoRev,altitude,ecc,inc,Om,om,f,5], init_deputy, mu,r_e,J2)
    

    #We know the orbit is periodic over 12 days, so just compute over those 12 days, every 10 seconds
    time = np.arange(0,12*24*3600,60)
    print("start")
    
 
    #Integrate the relative dynamics and chief orbital elements using pro_lib's dynamics function
    orbitState  = odeint(pro_lib.dyn_chief_deputies,ys,time,args=(mu,r_e,J2,num_deputy))


    #Plot the computed dynamics (set optional parameters to configure for animation)
    animationTools(orbitState, time)



    #Compute orbit cost

    #For displaying plots 3D visulaizations
    from mayavi import mlab
    cost = costFunction(time,orbitState,30,visualize=True)
    print("Cost:")
    print(cost)


    mlab.show()

def optimize():
    """
    Function that searches for an optimal formation for the weighting
    of science objectives versus set up costs. Utilizes a genetic algorithm
    to mutate the given swarm initial configuration and search for better solutions.
    As this is a random method, a sub optimal solution will be found, and there is
    no guarantee a better solution will be found than those provided.

    Currently all weightings are hard coded.
    
    """


    #Initialize the first guesses
    bestInit = np.array([[-0.11474911, -0.28277769, -0.41374357],
                              [ 0.11759513,  0.35926752, -0.29420415],
                              [-0.4380676,   0.13452622,  0.16840377],
                              [-0.80578882,  0.30974955,  0.41924403],
                              [ 0.18555127,  0.55242093,  0.95920626]])
    secondBestInit = np.array([[-8.45096242e-02, -3.19864051e-01, -5.20371558e-01],
                              [-4.31575222e-04,  3.79166091e-01, -2.77960141e-01],
                              [-2.93617169e-01,  7.21188588e-02,  1.57492712e-01],
                              [-7.38234180e-01,  3.49160003e-01,  4.40035707e-01],
                              [ 2.71468850e-01,  7.09862976e-01,  9.35197969e-01]])
    thirdBestInit = np.array([[ 0.04914105, -0.38146571, -0.51124651],
                         [ 0.0730455,   0.06972367, -0.50061731],
                         [-0.19072884,  0.09952473,  0.01128832],
                         [-0.51840071,  0.35818882,  0.47213096],
                         [ 0.27650221,  1.00772683,  0.73223388]])
    numDeputies = len(bestInit)

    #Time we integrate each orbit over
    time = np.arange(0,12*24*3600,60)
     
    # orbital parameters, wrapped in a dictionary

    orbParams = {"time":time, #time steps over which to evaluate
                "NoRev":173, # revolutions considered, orbit altitude, orbit 
                "altitude":747, #orbit altitude
                "ecc":0, #orbit eccentricity
                "inc":98.4, #orbit inclination (deg)
                "Om":0, #orbit right ascension of ascending node (deg)
                "om":0, #orbit argument of periapsis (deg)
                "f":0, #orbit true anomaly (deg)
                "num_deputy":5, 
                "lookAngle":30, #(deg)
                "mu":mu,  #gravitational parameter of Earth
                "r_e":r_e,  # Earth Radius 
                "J2":J2} #J2 Constant of Earth

    #Mutate, evaluate, and select
    for iterations in range(10):
        #Mutate to get states to evaluate:
        statesToEval = np.zeros((10,numDeputies,3))
        
        #First 5 states are bestInit and mutations, next 3 are second best, last 2 are thirdBest
        statesToEval[0] = bestInit
        statesToEval[5] = secondBestInit
        statesToEval[8] = thirdBestInit

        statesToEval[1:5] = mutate(bestInit,4)
        statesToEval[6:8] = mutate(secondBestInit,2)
        statesToEval[9] = mutate(thirdBestInit,1)

        #Pack up these and the orbit params in a tuple for multiprocessing
        params = [[statesToEval[0],orbParams]]
        for i in range(len(statesToEval)-1):
            params.append([statesToEval[i+1],orbParams])
        params = tuple(params)

        #Now loop through and score!
        #Will do this using multi-processing to parallelize
        #Using Pool as a context manager to ensure we close out properly 
        with Pool(processes=10) as pool:
            #Applies each state to the evaluate method and submits to the pool,
            #then waits for them to all return
            costs = pool.starmap(evalOrbit,params)
        print(costs)

        #Evaluate what has the lowest cost

        #Get indexes of lowest cost
        costs = np.array(costs)
        indexes = costs.argsort()[:3]

        #reassign
        bestInt = statesToEval[indexes[0]]
        secondBestInit = statesToEval[indexes[1]]
        thirdBestInit = statesToEval[indexes[2]]

    #Now print out the optimal formation and visualize!
    print("Best Initial Conditions:")
    print(bestInt)

    #Heavy module, import locally
    from mayavi import mlab
    

    #Compute orbit and cost again for display
    orbitState = computeOrbitDynamics(bestInit,orbParams) 
    cost =  costFunction(orbParams["time"],orbitState,30,visualize=True)

    print("Cost:")
    print(cost)

    #Show visualization
    #Blocks. Could multi thread to show the plot at the same time, or view one at a time as now.
    mlab.show()

    
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
    for i in range(orbParams["num_deputy"]):
    
        ax.plot(orbitState[:,6*(i+1)],orbitState[:,6*(i+1)+1],orbitState[:,6*(i+1)+2])
    
    ax.plot([0],[0],[0],"ko")
    plt.show()

def exportOrbit(initDeputy,orbParams):
    """
    A function to take a initial formation input and export it
    to a local .csv file at a 10kHz rate via interpolation. 
    Output is converted to ECI coordinates.

    Parameters
    ----------
    initDeputy : array, shape(num_deputy,3)
        Initial deputy positions to initialize the formation from
    orbParams : dict
        Dictionary of the orbital parameters used to compute
        the orbit. Required values are:
            time, NoRev, altitude, ecc, inc, Om, om, f, 
            num_deputy, lookAngle, mu, r_e, J2

            These are: time steps over which to evaluate, 
            revolutions considered, orbit altitude, orbit 
            eccentricity, orbit inclination (deg), orbit 
            right ascension of ascending node (deg), orbit
            argument of periapsis (deg), orbit true anomaly 
            (deg), number of deputies, look angle (deg), 
            earth gravitational parameter, radius of earth,
            earth J2
    """

    #Compute orbit
    orbitState = computeOrbitDynamics(initDeputy,orbParams)

    #TO DO: Account for frame rotation rate? Feel like we're missing terms here
    #Accelerations
    accelerations = np.zeros((len(orbitState),orbParams["num_deputy"]+1,3))

    #Attitude quaternions
    attitudes = np.zeros((len(orbitState),orbParams["num_deputy"]+1,4))

    #Reshape orbitState into a position and velocity array with same dimensions as acclerations and attitudeds
    positions = np.zeros((len(orbitState),orbParams["num_deputy"]+1,3))
    velocities = np.zeros((len(orbitState),orbParams["num_deputy"]+1,3))

    #Convert to ECI    
    for i in range(len(orbitState)):
        

        
        #Get chief position and velocity
        positions[i,0] = orbitState[i,0] * np.dot(pro_lib.rotation_matrix_lvlh_to_eci(orbitState[i,3],orbitState[i,5],orbitState[i,4]),np.array([1,0,0]))

        #Compute chief velocity from knowing the angular momentum (orbitStat[i,2]) and radial velocity (orbitStat[i,1])
        velocities[i,0] = np.dot(pro_lib.rotation_matrix_lvlh_to_eci(orbitState[i,3],orbitState[i,5],orbitState[i,4]),np.array([orbitState[i,1],orbitState[i,2]/orbitState[i,0],0]))

        #Compute acceleration vector
        orbitStateChange = pro_lib.dyn_chief_deputies(orbitState[i],orbParams["time"][i],orbParams["mu"],orbParams["r_e"],orbParams["J2"],orbParams["num_deputy"])
        
        #Need to compute the chief acceleration, take advantage of fact that we have vx dot (orbitStateChange[1]), r dot (orbitStateChange[0]) and h dot (orbitStateChange[2]) 
        accelerations[i,0] = np.dot(pro_lib.rotation_matrix_lvlh_to_eci(orbitState[i,3],orbitState[i,5],orbitState[i,4]),
                                         np.array([orbitStateChange[1],orbitStateChange[2]/orbitState[i,0]-orbitState[i,2]/(orbitState[i,0]**2)*orbitStateChange[0],0]))

        #Convert each deputy state
        for j in range(orbParams["num_deputy"]):
            #Rotate to ECI frame using the rotation matrix
            #Rotate position and add chief position
            positions[i,j+1] = np.dot(pro_lib.rotation_matrix_lvlh_to_eci(orbitState[i,3],orbitState[i,5],orbitState[i,4]),orbitState[i,6*(j+1):6*(j+1)+3]) + positions[i,0]
            #Rotate velocity and add chief velocity
            velocities[i,j+1] = np.dot(pro_lib.rotation_matrix_lvlh_to_eci(orbitState[i,3],orbitState[i,5],orbitState[i,4]),orbitState[i,6*(j+1)+3:6*(j+2)]) + velocities[i,0]
            #Rotate acceleration and add chief acceleration
            accelerations[i,j+1] = np.dot(pro_lib.rotation_matrix_lvlh_to_eci(orbitState[i,3],orbitState[i,5],orbitState[i,4]),orbitStateChange[6*(j+1)+3:6*(j+2)]) + accelerations[i,0]
        
        #Compute the quaternion associated with this orientation by using the built in rotation class to convert the pro_lib rotation matrix to a quaternion
        #Computes using the roation matrix from eci to lvlh

        orientation =  Rotation.from_matrix(pro_lib.rotation_matrix_eci_to_lvlh(orbitState[i,3],orbitState[i,5],orbitState[i,4]))
        attitudes[i] = orientation.as_quat()

    #Write to netCDF

    output = nc.Dataset('orbitOutput.nc','w', format='NETCDF4')

    time = output.createDimension('time', None)
    spaceCraft = output.createDimension('spaceCraft',orbParams["num_deputy"]+1)
    spatialDim = output.createDimension('spatialDim',3)
    quatDim = output.createDimension('quatDim',4)

    timeOut = output.createVariable("time","f4",("time",))
    positionOut = output.createVariable("position","f4",("time","spaceCraft","spatialDim",))
    VelocityOut = output.createVariable("velocity","f4",("time","spaceCraft","spatialDim",))
    accelOut = output.createVariable("acceleration","f4",("time","spaceCraft","spatialDim",))
    attOut = output.createVariable("attitude","f4",("time","spaceCraft","quatDim",))
    
    timeOut[:] = orbParams["time"]
    positionOut[:] = positions
    VelocityOut[:] = velocities
    accelOut[:] = accelerations
    attOut[:] = attitudes


def evalOrbit(state,orbParams):
    """
    Method to compute the cost of a proposed initial formation 
    position over the course of the orbit
    
    Parameters
    ----------
    state : array, shape(numDeputies,3)
        initial deputy spatial configurations of the
        swarm. Should be (x,y,z) of each deputy in order
    orbParams : dict
        Dictionary of the orbital parameters used to compute
        the orbit. Required values are:
            time, NoRev, altitude, ecc, inc, Om, om, f, 
            num_deputy, lookAngle, mu, r_e, J2

            These are: time steps over which to evaluate, 
            revolutions considered, orbit altitude, orbit 
            eccentricity, orbit inclination (deg), orbit 
            right ascension of ascending node (deg), orbit
            argument of periapsis (deg), orbit true anomaly 
            (deg), number of deputies, look angle (deg), 
            earth gravitational parameter, radius of earth,
            earth J2

    Returns
    ---------
    cost : double
        Cost of the orbit associated with this initial position
    """

    print("Entered eval")
    #Get orbit dyanmics
    orbitState = computeOrbitDynamics(state,orbParams)
    
    print("left eval")
    return costFunction(orbParams["time"],orbitState,orbParams["lookAngle"])

def computeOrbitDynamics(state,orbParams):
    """
    Method to compute the orbit dynamics of from an initial formation 
    position

    Parameters
    ----------
    state : array, shape(numDeputies,3)
        initial deputy spatial configurations of the
        swarm. Should be (x,y,z) of each deputy in order
    orbParams : dict
        Dictionary of the orbital parameters used to compute
        the orbit. Required values are:
            time, NoRev, altitude, ecc, inc, Om, om, f, 
            num_deputy, lookAngle, mu, r_e, J2

            These are: time steps over which to evaluate, 
            revolutions considered, orbit altitude, orbit 
            eccentricity, orbit inclination (deg), orbit 
            right ascension of ascending node (deg), orbit
            argument of periapsis (deg), orbit true anomaly 
            (deg), number of deputies, look angle (deg), 
            earth gravitational parameter, radius of earth,
            earth J2

    Returns
    ---------
    orbitState : array shape(len(time),6*(num_deputy+1))
        State vector of the orbit at each time specified.
        First 6 states are the chief orbital parameters.
        Each subsequent 6 states are a deputy's relative
        state in LVLH as (x,y,z,vx,vy,vz)
    """


    #compute initial conditions for the chief and deputy
    ys = pro_lib.initial_conditions_deputy("nonlinear_correction_linearized_j2_invariant",
                                            [orbParams["NoRev"],orbParams["altitude"],orbParams["ecc"],orbParams["inc"],orbParams["Om"],orbParams["om"],orbParams["f"],orbParams["num_deputy"]],
                                            state,orbParams["mu"],orbParams["r_e"],orbParams["J2"])
    
    #Integrate the relative dynamics and chief orbital elements using pro_lib's dynamics function
    orbitState  = odeint(pro_lib.dyn_chief_deputies,ys,orbParams["time"],args=(orbParams["mu"],orbParams["r_e"],orbParams["J2"],orbParams["num_deputy"]))
    return orbitState

def mutate(states,numOffspring,stdDeviation=.05):
    """
    Generates new random initial swarm configurations
    given a state to mutate from. The new states will
    be arranged in Gaussian fashion around the given 
    state with the specified standard deviation.

    Parameters
    ----------
    states : array, shape(numDeputies,3)
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
    offspring = np.zeros((numOffspring,len(states),3))

    #Now create random children 

    for i in range(len(offspring)):

        for j in range(len(states)):
            offspring[i,j,0] = np.random.normal(states[j,0],stdDeviation)
            offspring[i,j,1] = np.random.normal(states[j,1],stdDeviation)
            offspring[i,j,2] = np.random.normal(states[j,2],stdDeviation)

    return offspring

def animate(i,orbitData,ax,title):
    """
    Method to animate the ith frame of the orbit visualization
    by drawing the position of each space craft at this time step

    Parameters
    ----------
    i : int
        Frame to animate. This is the element of the orbit data 
        used to plot the space craft's positions
    orbitData : array, shape(3,len(time),num_deputy)
        The positions of each deputy spacecraft. The arrangement 
        is set to make slicing the original orbit data easier.
        The first dimension corresponds to x,y,z
        The second dimension corresponds to time step
        The third dimension corresponds to the deputy
    ax : matplotlib.pyplot axis
        The active axis we are using for animating with the 
        matplotlib.animation library 
    title : string
        The title for the figure, as we clear the axis each
        step to trace only the points, not the full line
    """

    #Print out progress
    print(i)
    #Plot the 
    ax.clear()
    ax.set_title(title)
    ax.set_xlabel("x, radial out from Earth (km)")
    ax.set_ylabel("y, along track (km)")
    ax.set_zlabel("z, cross track (km)")
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_zlim(-1, 1)
    data = orbitData[:,int(i):int(i+1),:] #select data range
    for j in range(len(data[0,0,:])):
        ax.plot(data[0,:,j],data[1,:,j],data[2,:,j],"o")
    
    ax.plot([0],[0],[0],"ko")

def costFunction(t,stateVector,lookAngle,visualize=False):
    """
    Return a cost value for the given orbit trajectory. This is the
    cost to set up the orbit minus the science merit of the orbit.
    The science merit is a per day score of the observational value
    for the formation. For best results, give a representative length
    of time (ie, time to repeat for a periodic orbit), and the score will
    be normalized per day using the last time in t.

    Parameters
    ----------
    t : array
        Time at which each state is provided, in seconds from the start of
        the computation
    stateVector : array, shape (len(t), len(state))
        State throughout the orbit at each time step. State should be 
        (x,y,z,vx,vy,vz) stacked for each space craft in LVLH frame.
        First s/c is the chief and should be instead the Xu Wang Parameters,
        (r,vx,h,Omega,inc,theta)
    lookAngle : double
        Angle at which the formation is looking towards its target, defined 
        as a right-handed rotation around the along track (+y) direction. 
        0 is the nadir (-x) direction
    visualize : Boolean (default=False)
        When true, a 3D rendering of the orbit is generated

    Returns
    -------
    J : double
        Numerical cost of the orbit, as setup cost minus science merit. 
        A lower cost means a more efficient/desirable orbit.
    
    """

    #Cost function we are trying to evaluate
    #J = e^(deltaV to initialize) - Sum[(targetHeight[i])/(5*resolution[i])]/numDays
    #Subject to constraints: (enforced by large cost penalty)
    #   resolution[i] < targetHeight[i]/5
    #   ambiguity[i] > targetHeight[i]

    #Compute the cost to set up the orbit, and maintain it
    orbitCost = computeOrbitCost(t,stateVector)
    #print(orbitCost)

    #Compute the science merit that offsets the cost
    scienceMerit = computeScienceMerit(t,stateVector,lookAngle,visualize)

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
        (x,y,z,vx,vy,vz) stacked for each deputy in LVLH. First 6 elements
        are assumed to be the chief orbit in Xu Wang Parameters and are
        ignored
    Returns
    -------
    J : double
        Numerical cost of the orbit to setup. A lower cost means a more 
        efficient/desirable orbit.
    
    """
    
    #Compute the number of space craft in the formation
    numSC = int(len(stateVector[0])/6)

    #Compute delta-V to initialize orbit
    deltaV = 0
    for i in range(numSC - 1):
        initialV = stateVector[0,3+6*(i+1):6+6*(i+1)] * 1000 #Want in m/s
        #print(initialV)
        deltaV += np.linalg.norm(initialV)

    #Exponentiate to give cost
    return np.exp(deltaV)

def computeScienceMerit(t,stateVector,lookAngle=0,visualizeTrajectory=False):
    """
    Return a science merit score for the given orbit trajectory. 
    Intended to be over a single orbit 

    Assumptions: Earth is modeled as a perfect sphere for this 
    approximation with radius 6378.1363 km as set at the module level
    
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
    visualize : Boolean (default=False)
        When true, a 3D rendering of the orbit is generated
    
    Returns
    -------
    score : double
        Numerical evaluation of the science merit of the orbit. A higher
        value indicates higher resolution compared to the target heights
    """
    
    #Define targets/ s/c position to see them

    #Load the elevation data in from https://webmap.ornl.gov/ogc/dataset.jsp?ds_id=10023 at .25 degree resolution
    #If finer resolution is desired the hard coding will need to be changed
    elevationData = np.genfromtxt('elevationData.txt',delimiter=' ',skip_header=5) #Units, meters
    #bottom left corner corresponds to -180 W -90 S in this data, so upper right is 179.75E 89.75N

    #Compute the ground track of the center of the swath of the target location over the orbit
    #Time,ground track
    targetGroundTracks = np.zeros((len(t),2))
    
    
    #Don't create these unless we're going to plot them
    if visualizeTrajectory:
        groundTracks = np.zeros((len(t),2))

    #Distances to the target
    r0s = np.zeros(len(t))
    for i in range(len(targetGroundTracks)):
        #translating back to state vector from the orbital elements
        chiefState = stateVector[i,0] * np.dot(pro_lib.rotation_matrix_lvlh_to_eci(stateVector[i,3],stateVector[i,5],stateVector[i,4]),[1,0,0])
        #Compute along track unit vector
        yhat = np.dot(pro_lib.rotation_matrix_lvlh_to_eci(stateVector[i,3],stateVector[i,5],stateVector[i,4]),[0,1,0])
        #Compute where we are looking
        (targetGroundTracks[i],r0s[i]) = groundAngleTrackComputation(chiefState,yhat,t[i],lookAngle)

        if visualizeTrajectory:

            #Compute ground track with no look angle for visualization as well
            groundTracks[i] = groundAngleTrackComputation(chiefState,yhat,t[i],0)[0]

    
    #When over target, compute baseline, ambiguity
    baselines = np.zeros(len(t))
    separations = np.zeros(len(t))
    vegH = np.zeros(len(t))

    #Use the current target ground position to determine if we are over a target
    for i,position in enumerate(targetGroundTracks):
        #Find lat/lon of the chief
        lat = position[0]
        lon = position[1]

        #Compare with the vegetation data to see if we are over a target
        #Will compute where to look in elevationData matrix
        #Bound by extremes of array to avoid falling of the end of the array
        row = int(np.min((719,np.max((0,np.round((-lat+89.75)*4))))))
        col = int(np.min((1439,np.max((0,np.round((lon+180)*4))))))
        vegH[i] = elevationData[row,col]
        
        #Compute the baseline
        baselines[i] = computeBaseline(stateVector[i],lookAngle)
        #Compute the ambiguity
        separations[i] = computeMaxSeperation(stateVector[i],lookAngle)
    
    #Function for visualizing the ground tracks
    if visualizeTrajectory:
        #visualize(np.radians(targetGroundTracks[:,0]),np.radians(targetGroundTracks[:,1]),np.radians(groundTracks[:,0]),np.radians(groundTracks[:,1]),elevationData,baselines,separations)
        #Want to not have the groundtracks, for less clutter
        visualize(np.radians(targetGroundTracks[:,0]),np.radians(targetGroundTracks[:,1]),elevationData=elevationData,baselines=baselines,separations=separations)


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
        ambiguities[i] = lam * r0s[i] / (2 * separations[i])
        #Check if over a  target!
        #Note if violate baseline/ambiguity constraint, and how often
        if vegH[i]  > 0:
            if resolutions[i] > vegH[i]/5:
                numViolateRes +=1
                score -= 50000 #Heavily penalize constraint violations 
            if ambiguities[i] < vegH[i]:
                numViolateAmb +=1
                score -= 50000 #Heavily penalize constraint violations
            score += (vegH[i]/resolutions[i])/5
    
    print(resolutions)   
    print(f"Time over targets: {len(np.where(vegH > 0)[0])/len(vegH)*100}%")

    print("Violations of Scientific Constraints")
    print(f"Resolution Violations (>1/5 target height): {numViolateRes}")
    print(f"Ambiguity Violations (< target height): {numViolateAmb}")
    return score


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


    #Compute the radial and unit vector
    rhat = x/np.linalg.norm(x) #Also xhat, as this is the x direction in LVLH

    #Check if lookAngle != 0
    if lookAngle == 0:
        targetVector = -rhat
        r0 = np.linalg.norm(x) - r_e
    else:

        #Convert look angle to rads
        rlookAngle = np.radians(lookAngle)

        



        #Compute the view unit vector by rotating the -rhat vector around the yhat vector by the look angle
        rotate = Rotation.from_rotvec(rlookAngle * yhat)
        lookVector = rotate.apply(-rhat) #UNUSED

        #Compute the angle we need to rotate rhat by to get look target using law of sines
        #(we know the look angle and its opposite side (radius of the Earth)
        wideAng = np.arcsin(np.linalg.norm(x)* np.sin(np.abs(rlookAngle))/(r_e))
        rotationAngle =  np.pi - (wideAng + rlookAngle)
        #Compute r0, distance from chief to the target, also using the law of sines
        r0 = np.sin(rotationAngle)*(r_e)/np.sin(rlookAngle)

        #Check if physical or we got the wrong quadrant
        if r0 > np.linalg.norm(x):
            #wideAng was in wrong quadrant so should have been wideAng = np.pi - wideAng
            rotationAngle = wideAng - rlookAngle
            r0 = np.sin(rotationAngle)*(r_e)/np.sin(rlookAngle)
        
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

def visualize(targetLatitude,targetLongitude,latitude=None,longitude=None,elevationData=None,baselines=None,separations=None):
    """
    Plots a 3D visualization of the target ground track,
    including optional visualizations of targets,
    baselines and separations

    Parameters
    ----------
    latitude : array shape (len(t))
        list of target latitudes at each time to 
        visualize, targets are where the formation
        is looking
    longitude : array shape (len(t))
        list of target longitudes at each time to visualize
    latitude : array shape (len(t)) (optional)
        list of ground track latitudes at each time to 
        visualize
    longitude : array shape (len(t)) (optional)
        list of ground track longitudes at each time to 
        visualize
    elevationData : array shape(len(t)) (optional)
        list of vegetation heights to superimpose
    baselines : array shape(len(t)) (optional)
        list of baselines at each time step to visualize
    separations : array shape(len(t)) (optional)
        list of separations at each time step to visualize
    """
    
    #Heavy module, import locally
    from mayavi import mlab
    

   
    #Compute target and s/c ground tracks
    targetGroundTrack = np.array([np.cos(targetLongitude)*np.cos(targetLatitude),np.sin(targetLongitude)*np.cos(targetLatitude),np.sin(targetLatitude)])#*earthRadius
    
    #Check if we are using the formation ground tracks, then plot them
    if latitude is not None:
        groundTrack  = np.array([np.cos(longitude)*np.cos(latitude),np.sin(longitude)*np.cos(latitude),np.sin(latitude)])#*.95#*earthRadius
    
    
    #Plot a sphere for the earth
    u, v = np.mgrid[0:2*np.pi:100j, 0:np.pi:50j]
    x = np.cos(u)*np.sin(v)#*earthRadius
    y = np.sin(u)*np.sin(v)#*earthRadius
    z = np.cos(v)#*earthRadius

    #Create a Mayavi mlab figure. This is a package that provides matlab-like plotting capabilites, and better
    #3D capabilities than matplot lib. Requires vtk to be installed for visualization, and may require a downgraded
    #version of vtk to run (9.0 instead of 8.2 caused issues, developed with vtk 8.1 and Mayavi 4.7.1
    mlab.figure()
    s = mlab.mesh(x, y, z,colormap="pink") #Pick a neutral color map for the planet
    
    #Check if we are plotting vegetation
    if elevationData is not None: 
        #Compute where the elevation data goes and plot it
        #First find where their is vegetation
        mask = np.where(elevationData >0)
        #Create flat arrays for lat long and veg height
        flatVegData = np.zeros((len(mask[0]),3))

        #Loop through and flatten
        #bottom left corner corresponds to -180 W -90 S in this data, so upper right is 179.75E 89.75N
        for i in range(len(mask[0])):
            #Compute lat
            flatVegData[i,0] = -.25 * mask[0][i] + 89.75
            #Compute lon
            flatVegData[i,1] = .25 * mask[1][i] - 180
            #Get vegH
            flatVegData[i,2] = elevationData[mask[0][i],mask[1][i]]
        
        
            
        #Compute the 3d coords
        vegLocations = np.array([np.cos(np.radians(flatVegData[:,1]))*np.cos(np.radians(flatVegData[:,0])),
                                 np.sin(np.radians(flatVegData[:,1]))*np.cos(np.radians(flatVegData[:,0])),
                                 np.sin(np.radians(flatVegData[:,0]))])
        
        mlab.points3d(vegLocations[0],vegLocations[1],vegLocations[2],flatVegData[:,2],scale_factor = .0005)

    #Now plot the ground tracks

    #Just plot the target ground track
    if baselines is None and separations is None:
        mlab.plot3d(targetGroundTrack[0],targetGroundTrack[1],targetGroundTrack[2],tube_radius=None,color=(0,0,0)) 
        mlab.title("Target Ground Track Over Orbit")
    #At least one of the features should be shown along the target ground track
    elif baselines is not None:
        mlab.plot3d(targetGroundTrack[0],targetGroundTrack[1],targetGroundTrack[2],baselines,tube_radius=None) 
        mlab.title("Target Ground Track Over Orbit, colored to indicate baseline")
        #Plot second figure!
        if separations is not None:
            mlab.figure()
            s = mlab.mesh(x, y, z,colormap="pink") #Pick a neutral color map for the planet
            mlab.plot3d(targetGroundTrack[0],targetGroundTrack[1],targetGroundTrack[2],separations,tube_radius=None) 
            #Replot veg data too if necessary
            if elevationData is not None: 
                mlab.points3d(vegLocations[0],vegLocations[1],vegLocations[2],flatVegData[:,2],scale_factor = .0005)
            mlab.title("Target Ground Track Over Orbit, colored to indicate separations")
    #Just the separations should be shown
    else:
        mlab.plot3d(targetGroundTrack[0],targetGroundTrack[1],targetGroundTrack[2],separations,tube_radius=None) 
        mlab.title("Target Ground Track Over Orbit, colored to indicate separations")
    
    
    
    #Plot the swarm ground track too (only on one of the figures, if two made)
    if latitude is not None:
        mlab.plot3d(groundTrack[0],groundTrack[1],groundTrack[2],tube_radius=None,color=(1,0,0))







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

def orbitPeriodComputation(orbParams,timeStepsPerOrbit):
    """
    This code computes the time for a full orbit, accounting 
    for J2 perturbations

    Parameters
    ----------
    orbParams : dict
        Dictionary of the orbital parameters used to compute
        the orbit. Required values are:
            NoRev, altitude, ecc, inc,
            num_deputy, lookAngle, mu, r_e, J2

            These are: 
            revolutions considered, orbit altitude, orbit 
            eccentricity, orbit inclination (deg),
            number of deputies, look angle (deg), 
            earth gravitational parameter, radius of earth,
            earth J2
    timeStepsPerOrbit : int
        Number of time steps to be computed per orbit, which
        determines the spacing of the time steps in the 
        returned time vector

    Returns
    -------
    time : array, shape(NoRev*timeStepsPerOrbit)
        Time steps along the orbit. Last element is the 
        computed timespan

    """

    # Energy Matched 
    # J2 disturbance 
    # No drag

    k_J2 = (3/2)*orbParams["J2"]*orbParams["mu"]*(orbParams["r_e"]**2)

    # Orbital Elements
    a = orbitParams["r_e"] + orbitParams["altitude"]           # semimajor axis [km] (Assumed circular orbit)
    inc = orbitParams["inc"]*np.pi/180             # inclination [rad]


    # Xu Wang Parameters
    h = np.sqrt(a*(1 - orbitParams["ecc"]**2)*orbitParams["mu"])           # angular momentum [km**2/s]


    # effective period of the orbit
    a_bar = a*(1 + 3*orbitParams["J2"]*orbitParams["r_e"]**2*(1-orbitParams["ecc"]**2)**0.5/(4*h**4/orbitParams["mu"]**2)*(3*np.cos(inc)**2 - 1))**(-2/3)  
    period = 2*np.pi/np.sqrt(orbitParams["mu"])*a_bar**(3/2)  
 

    # simulation time
    time = np.linspace(0,orbitParams["NoRev"]*period,int(period/(timeStepsPerOrbit)))  
    
    return time               # time vector with units of orbits instead of seconds

    


def animationTools(orbitState, time,azim=-100, elev=43, animate=False,frames=None,animationName="animation.mp4",sliders=False):
    """
    Helper method to animate or provide lightweight 
    visualization of the formation dynamics. Several
    optional parameters configure the type of visualization
    or animation displayed

    Parameters
    ----------
    orbitState : array, shape(len(time),6*(num_deputy+1)
        State vector of the orbit at each time specified.
        First 6 states are the chief orbital parameters.
        Each subsequent 6 states are a deputy's relative
        state in LVLH as (x,y,z,vx,vy,vz) 
    time : array
        Time at which each state is provided
    azim : double, (default=-100)
        Azimuth angle of the initial plot rendering
    elev : double, (default=43)
        Elevation angle of the initial plot rendering
    animate : Boolean, (default=False)
        Flag to animate the formation over the orbit
    frames : int, (default=None)
        If animating, how many frames to animate. 
        Default of none animates full orbit
    animationName : string, (default="animation.mp4")
        If animating, name of output file (in local
        directory). NOTE: No overwrite protection!
    sliders : boolean, (default=False)
        Flag to produce plot with interactive sliders
        for the formation over its orbit history.
    """
    

    #Plot the relative orbit tracks, at a provided or arbitrary view angle (found to work well for these visualizations)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_title("6 space craft formation in NISAR J2 dynamic orbit, LVLH frame")
    ax.set_xlabel("x, radial out from Earth (km)")
    ax.set_ylabel("y, along track (km)")
    ax.set_zlabel("z, cross track (km)")
    ax.set_xlim(-1, 1)
    ax.set_ylim(-1, 1)
    ax.set_zlim(-1, 1)
    ax.azim = azim
    ax.elev = elev

    #Loop through each deputy
    for i in range(int(len(orbitState[0])/6-1)):
    
        ax.plot(orbitState[:,6*(i+1)],orbitState[:,6*(i+1)+1],orbitState[:,6*(i+1)+2])
    
    ax.plot([0],[0],[0],"ko")
    plt.show()
    
    if animate or sliders:
        #Save the user selected "best" veiw for animation
        azimuth = ax.azim
        elevation = ax.elev
    
    #Show the orbit controlled by sliders (sort of works) if desired, so the user can manipulate the dynamics
    if sliders:
        fig = go.Figure()
        # Add traces, one for each slider step
        for state in orbitState:
            xs = state[6::6]
            ys = state[7::6]
            zs = state[8::6]
            fig.add_trace(
                go.Scatter3d(
                    visible=False,
                    x=xs,y=ys,z=zs),
                    range_x=[-1,1], range_y=[-1,1], range_z=[-1,1])
        
        # Make 0th trace visible
        fig.data[0].visible = True
        
        # Create and add slider
        steps = []
        for i in range(len(fig.data)):
            step = dict(
                method="update",
                args=[{"visible": [False] * len(fig.data)}],
            )
            step["args"][0]["visible"][i] = True  # Toggle i'th trace to "visible"
            steps.append(step)
        
        sliders = [dict(
            active=0,
            steps=steps
        )]
        
        fig.update_layout(
            sliders=sliders
        )
        
        
        fig.show()
    
    #Animate if desired
    if animate:
        #Check if user specified number of frames, or animate whole thing
        if frames is None:
            frames = len(time)


        #Only import when animating
        import matplotlib.animation as animation
        orbitData = np.array([orbitState[:,6::6],orbitState[:,7::6],orbitState[:,8::6]])
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)
        
        fig = plt.figure(figsize=(10,10))
        ax = fig.add_subplot(111, projection='3d')
        ax.set_xlim(-1, 1)
        ax.set_ylim(-1, 1)
        ax.set_zlim(-1, 1)
        ax.azim = azimuth
        ax.elev = elevation
        ani = animation.FuncAnimation(fig, animate, frames=frames, fargs=(orbitData,ax,"6 space craft formation in NISAR J2 dynamic orbit, LVLH frame"))
        
        ani.save(animationName, writer=writer)




if __name__ == "__main__":
    main()

   