import numpy as np
import pandas as pd
import math
import itertools


# list of all the minutiae points
list_of_minutiae = [(3,2,0.5), (-4,1.2,3.4), (5,-2,4.4), (1,-1,1), (3.6,2,0.5), (1,1,1), (-2,-0.3, 0), (1,2, -4)]

#print(itertools.combinations_with_replacement(list_of_minutiae, 2))

for mi,mj in itertools.combinations_with_replacement(list_of_minutiae, 2):
    def distance(mi,mj):
        dist = math.sqrt((mj[0] - mi[0])**2 + (mj[1] - mi[1])**2)
        return dist
    #print(mi, mj, "Distance :", dist)
    def unit_vector(list_of_minutiae):
        return list_of_minutiae / np.linalg.norm(list_of_minutiae)  # this statement breaks when mj = (0,0,0)
    #print(unit_vector(mi))j
    mi_u = unit_vector(mi)
    mj_u = unit_vector(mj)

    def angle_between(mi,mj):
        angle = np.arccos(np.clip(np.dot(mi_u, mj_u), -1.0, 1.0))
        return angle
    #print (angle)

    print ("   ")
    print (mi, mj)
    print ("Distance: ", distance(mi,mj))
    print ("Angle: ", angle_between(mi,mj))

    # define X coordinate of Vicinity 1 (corresponding to mi and mj)
    X = distance(mi, mj) * (math.cos(angle_between(mi,mj) - mi[2]))
    print("X: ", X)

    # define Y coordinate
    Y = distance(mi, mj) * (math.sin(angle_between(mi,mj)-mi[2]))
    print("Y: ", Y)

    # define theta
    def theta(mi, mj):
        new_theta = mi[2] - mj[2]
        return new_theta
    theta = theta(mi, mj)
    print("Theta: ", theta)

    def vicinity(mi, mj):
        output = []
        output.append([X, Y, theta])
        return output
    print("Vicintiy: ", vicinity(mi, mj))
