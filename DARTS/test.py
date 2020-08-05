import sys 
import os
sys.path.append(os.path.abspath("C:\\Users\\trey_\\Documents\\GitHub\\PROpy"))

import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.spatial.transform import Rotation
import pro_lib
import numpy as np



def main():
    r = np.array([1,0,0])
    print(np.dot(pro_lib.rotation_matrix_lvlh_to_eci(np.pi/4,np.pi/4,np.pi/4),r))




if __name__ == "__main__":
    main()

