import sys 
import os
sys.path.append(os.path.abspath("C:\\Users\\trey_\\Documents\\GitHub\\PROpy"))

import matplotlib.pyplot as plt
import scipy
import pro_lib
import numpy as np




def main():
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
    stateVector : array, shape (len(t), len(state))
        State throughout the orbit at each time step. State should be 
        (x,y,z,vx,vy,vz) stacked for each space craft

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

 


def computeScienceMerit(t,stateVector):
    """
    Return a science merit score for the given orbit trajectory. 
    Intended to be over a single orbit 
    
    Parameters
    ----------
    t : array
        Time at which each state is provided
    stateVector : array, shape (len(t), len(state))
        State throughout the orbit at each time step. State should be 
        (x,y,z,vx,vy,vz) stacked for each space craft

    Returns
    -------
    J : double
        Numerical cost of the orbit. A lower cost means a more 
        efficient/desirable orbit.
    
    """

    #Pseudo code outline
    numSC = int(len(stateVector[0]/6))
    #Define targets/ s/c position to see them

    #Load the elevation data in
    #elevationData = 


    #Compute ground track over orbit, including the effect of earth rotating
    
    #Time,s/c,ground track
    groundTracks = np.zeros((len(t),numSC,2))
    
    for i in range(numSC):
        groundTracks = groundTrackComputation(stateVector[:,0+i*6:3+i*6])
        

    #When over target, compute baseline, ambiguity

    #Note if violate baseline/ambiguity constraint, and how often

    #Score

    
    #TODO, account for starting position of the Earth
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
    t0 : initial time 

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
    longitude = longitude - t* EARTH.rotD *pi/180;
    
    #stack up the output
    groundTrack = np.vstack(latitude,longitude)
    return groundTrack.transpose()

    

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

   