import sys 
import os
sys.path.append(os.path.abspath("C:\\Users\\trey_\\Documents\\GitHub\\PROpy"))

import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.spatial.transform import Rotation
import pro_lib
import numpy as np
#For parallel evaluation
from multiprocessing import Pool


def main():
    with Pool(processes=10) as pool:
        #Applies each state to the evaluate method and submits to the pool,
        #then waits for them to all return
        #costs = pool.starmap(test,((1,),(2,),(3,))) #,4,5,6,7,8,9,0
        pool.map(test, range(10))
    print("done")
    print(costs)


def test(i):
    print("starting process")
    print(i)
    print("ending process")
    return i


if __name__ == "__main__":
    main()

