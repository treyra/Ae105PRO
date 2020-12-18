import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.spatial.transform import Rotation
import pro_lib
import numpy as np

#Class constants (Earth parameters) can be overwritten locally or when calling the utility methods
mu = 398600.432896939164493230  #gravitational constant
r_e = 6378.136  # Earth Radius 
J2 = 0.001082627 #J2 Constant

#Currently the same as up to date DARTS, but used by coverage analysis in quick hack for visualization
def costFunction(t,stateVector,lookAngle,visualize=False,otherData=None):
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
    print("set up cost:")
    print(orbitCost)

    #Compute the science merit that offsets the cost
    scienceMerit = computeScienceMerit(t,stateVector,lookAngle,visualize,otherData)

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

    #Compute MAX delta-V to initialize orbit (as we want to know most expensive s/c to initialize)
    deltaV = 0
    for i in range(numSC - 1):
        initialV = stateVector[0,3+6*(i+1):6+6*(i+1)] * 1000 #Want in m/s
        print(initialV)
        deltaV = np.maximum(deltaV,np.linalg.norm(initialV))
        print(deltaV)

    #Exponentiate to estimate mass cost. Assume cold gas thrusters with ~150 ISP
    return np.exp(deltaV/150)


#Old Baseline and sepearations computation that seperatly computed these quantities

def computeBaseline(stateVector, lookAngle=0):
    """
    Return the baseline for the formation at the current time of the orbit.
    
    Parameters
    ----------
    stateVector : array
        State Vector of the swarm at the time to compute the baseline for.
        State should be (x,y,z,vx,vy,vz) stacked for each space craft. 
        x,y,z directions defined by the LVLH coordinate system, where
        x is radially out from the Earth to the s/c, y is along track (or 
        tangential velocity) and z is the orbit angular momentum or cross track.
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
            #Project onto the imaging plane (This is just removing the y component in this formulation)
            #If the look vector wasn't just a rotation of the nadir direction around the tangential veloctity, this would need to be more complicated
            candiadateBaseline[1] = 0
            #Project onto the look angle and subtract this off to get component perpendicular to the look angle
            candiadateLn = np.linalg.norm(candiadateBaseline - np.dot(candiadateBaseline,lookVector)*candiadateBaseline/np.linalg.norm(candiadateBaseline))
            if candiadateLn > bestLn:
                bestLn = candiadateLn
    return bestLn

def computeAverageSeperation(stateVector, lookAngle=0):
    """
    Return the average cross track separation to estimate the ambiguity.
    
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
        Average gap between any two space craft perpendicular to the look angle
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
        #This is projection onto the imaging plane
        #If the look vector wasn't just a rotation of the nadir direction around the tangential veloctity, this would need to be more complicated
        projectedPositions[i] = np.array([positions[i,0],positions[i,2]])

    #Now we know they are all along a straight line, so we can sort by our x axis! (Now each s/c is in order along the baseline)
    sortIndex = np.argsort(projectedPositions[:,0])
    sortedPositions = projectedPositions[sortIndex]

    #Now loop through recording seperations to take the average
    mus = np.zeros(len(sortedPositions)-1)
    for i in range(len(sortedPositions)-1):
        mus[i] = np.linalg.norm(sortedPositions[i+1]-sortedPositions[i])

    return np.median(mus)

#Old science Merit computation that took into account the targets.

def computeScienceMerit(t,stateVector,lookAngle=0,visualizeTrajectory=False,otherData=None):
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
        #Compute the median spacing
        separations[i] = computeAverageSeperation(stateVector[i],lookAngle)
    
    #Function for visualizing the ground tracks
    #TO DO: Roll in other visulaization methods?
    #TO DO: Add scales to the plots to quantify better
    #TO DO: Consider looking at the resolution and ambiguities (projected possibly) rather than baselines and spereations
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
    
    
    lam = .24 #Units in meters
    resolutions = np.zeros(len(t))
    ambiguities = np.zeros(len(t))
    numViolateRes = 0
    numViolateAmb = 0
    #Currently giving score as # of times better than constraint (1/5 of veg height)
    score = 0
    for i in range(len(t)):
        resolutions[i] = lam * r0s[i]*1000 / (2 * baselines[i]*1000) #Convert to meters
        ambiguities[i] = lam * r0s[i]*1000 / (2 * separations[i]*1000)
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
    
    #Compute additional data as requested to pass out (such as min resolution, max baseline, etc.)
    if otherData is not None:
        #Other data is a dictionary of desired parameters. Check what is requested and compute
        if "maxResolution" in otherData:
            #Global max resolution (worst case)
            otherData["maxResolution"] = np.max(resolutions)
        if "minBaseline" in otherData:
            #Global minimum baseline (worst case)
            otherData["minBaseline"] = np.min(baselines)
        if "minAmbiguity" in otherData:
            #Global minimum Ambiguity (worst case)
            otherData["minAmbiguity"] = np.min(ambiguities)
        if "maxSeparation" in otherData:
            #Global max separation (worst case)
            otherData["maxSeparation"] = np.max(separations)
            
    #Plot resolution and ambiguity over orbit
    if visualizeTrajectory:
        plt.figure()
        plt.plot(resolutions)
        plt.title("Resolution over the orbit")
        plt.xlabel("Time (s)")
        plt.ylabel("Resolution (m)")
        plt.figure()
        plt.plot(ambiguities)
        plt.title("Nearest ambiguity over the orbit")
        plt.xlabel("Time (s)")
        plt.ylabel("Nearest Ambiguity (m)")
        plt.show()


    
    print(resolutions)   
    print(f"Time over targets: {len(np.where(vegH > 0)[0])/len(vegH)*100}%")

    print("Violations of Scientific Constraints")
    print(f"Resolution Violations (>1/5 target height): {numViolateRes}")
    print("Percentage of orbit below 1m resolution")
    print(len(np.where(resolutions <= 1)[0])/len(resolutions))

    print(f"Ambiguity Violations (< target height): {numViolateAmb}")
    print("Percentage of orbit above 30m ambiguity")
    print(len(np.where(ambiguities >= 1)[0])/len(ambiguities))
    return score

###############
#Used for visulaization and computing the target location, not needed in lighterweight analysis.


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
    targetLatitude : array shape (len(t))
        list of target latitudes at each time to 
        visualize, targets are where the formation
        is looking
    targetLongitude : array shape (len(t))
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
