import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.spatial.transform import Rotation
import pro_lib
import numpy as np

#Test import for comparing against legacy functions
import legacyDARTSfunctions


#Class constants (Earth parameters) can be overwritten locally or when calling the utility methods
mu = 398600.432896939164493230  #gravitational constant
r_e = 6378.136  # Earth Radius 
J2 = 0.001082627 #J2 Constant

def main():
    visualizeOrbit(NoRev = 433,altitude = 741.846, days = 30,ecc = 0,inc = 98.3597,Om = 0,om = 0,f = 0) 



def visualizeOrbit(NoRev,altitude,inc,days,ecc = 0,Om = 0,om = 0,f = 0):
    """
    Visualizes the formation and target ground track for coverage anlysis.
    Resuses demo code and old Mayavi visulizations, so uses an arbitrary formation for now.
    """



    init_deputy = np.array([[-0.72184274, -1.79858705, -2.43919107],
                            [ 0.74070093,  2.26334727, -2.06676058],
                            [-2.76214932,  0.82574802,  1.22148641],
                            [-5.07317191,  1.92931908,  6.16591021],
                            [ 1.16866556,  3.48151115,  6.47005005]])

    num_deputy = len(init_deputy)

    # compute initial conditions for the chief and deputy
    ys = pro_lib.initial_conditions_deputy("nonlinear_correction_linearized_j2_invariant", 
                                           [NoRev,altitude,ecc,inc,Om,om,f,num_deputy], init_deputy, mu,r_e,J2)
    

    #We know the orbit is periodic over 12 days, so just compute over those 12 days, every 60 seconds
    time = np.arange(0,days*24*3600,60)
    print("start")
    print(ys)
 
    #Integrate the relative dynamics and chief orbital elements using pro_lib's dynamics function
    orbitState  = odeint(pro_lib.dyn_chief_deputies,ys,time,args=(mu,r_e,J2,num_deputy))

    #Plot the computed dynamics (set optional parameters to configure for animation)
    #animationTools(orbitState, time)

    #Create dictionary of other parameters to compute

    otherData = {"maxResolution" : None, "minBaseline" : None, "minAmbiguity" : None, "maxSeparation" : None }

    #Compute orbit cost

    #No longer dispay 3d visulaizations
    ##For displaying plots 3D visulaizations
    ##This is a package that provides matlab-like plotting capabilites, and better
    ##3D capabilities than matplot lib. Requires vtk to be installed for visualization, and may require a downgraded
    ##version of vtk to run (9.0 instead of 8.2 caused issues, developed with vtk 8.1 and Mayavi 4.7.1
    #from mayavi import mlab
    cost = legacyDARTSfunctions.costFunction(time,orbitState,30,visualize=True,otherData=otherData)
    print("Cost:")
    print(cost)

    print("Max Resolution")
    print(otherData["maxResolution"])
    print("Min Baseline")
    print(otherData["minBaseline"])
    print("Min Ambiguity")
    print(otherData["minAmbiguity"])
    print("Max Separation")
    print(otherData["maxSeparation"])


if __name__ == "__main__":
    main()