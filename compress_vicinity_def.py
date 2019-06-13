import numpy as np
import pandas as pd
import math
import itertools

from yapf.yapflib.yapf_api import FormatFile
FormatFile(
    "/Users/ivysandberg/dmc/NXGBCC/dev/compress_vicinity_def.py",
    in_place=True)

# list of all the minutiae points
list_of_minutiae = [(3, 2, 0.5), (-4, 1.2, 3.4), (5, -2, 4.4), (1, -1, 1),
                    (3.6, 2, 0.5), (1, 1, 1), (-2, -0.3, 0), (1, 2, -4)]

# define a radius
center_minutiae = list_of_minutiae[0]
other_minutiae = list_of_minutiae[1:]


def vicinity_within_radius(minutiae):
    list_of_minutiae_in_r = []
    r = 2.5
    for i in other_minutiae:
        x = i[0]
        y = i[1]
        if x**2 + y**2 <= r**2:
            list_of_minutiae_in_r.append(other_minutiae)
    return list_of_minutiae_in_r


print(
    "This is a list of all the mjs within radius = 2.5 when mi = (3, 2, 0.5) : "
)
print(vicinity_within_radius(list_of_minutiae))

# Create a change of coords for every combination of mi (center minutiae) and mj (other minutiae)
# print(itertools.combinations_with_replacement(list_of_minutiae, 2)) = every combination of coordinates

for mi, mj in itertools.combinations_with_replacement(list_of_minutiae, 2):
    print("This is the (mi, mj) coordinate being changed")
    print(mi, mj)

    # Functions needed to determine the change of coordinates when one coord is translated & rotated to the origin


    def distance(mi, mj):
        dist = math.sqrt((mj[0] - mi[0])**2 + (mj[1] - mi[1])**2)
        return dist

    #print(mi, mj, "Distance :", dist)
    def unit_vector(list_of_minutiae):
        return list_of_minutiae / np.linalg.norm(
            list_of_minutiae)  # this statement breaks when mj = (0,0,0)

    #print(unit_vector(mi))j
    mi_u = unit_vector(mi)
    mj_u = unit_vector(mj)

    def angle_between(mi, mj):
        angle = np.arccos(np.clip(np.dot(mi_u, mj_u), -1.0, 1.0))
        return angle

    #print (angle)

    # print ("   ")
    # print (mi, mj)
    # print ("Distance: ", distance(mi,mj))
    # print ("Angle: ", angle_between(mi,mj))
    ''' Find the new X coordinate '''

    # define X coordinate of Vicinity 1 (corresponding to mi and mj)
    X = distance(mi, mj) * (math.cos(angle_between(mi, mj) - mi[2]))
    #print("X: ", X)
    ''' Find the new Y coordinate '''
    # define Y coordinate
    Y = distance(mi, mj) * (math.sin(angle_between(mi, mj) - mi[2]))

    #print("Y: ", Y)
    ''' Find the new theta coordinate '''

    # define theta
    def theta(mi, mj):
        new_theta = mi[2] - mj[2]
        return new_theta

    theta = theta(mi, mj)

    #print("Theta: ", theta)


    def vicinity(mi, mj):
        output = []
        output.append([X, Y, theta])  # creates the double brackets
        #output = np.array([X, Y, theta])    # has no commas
        return output

    #print("Vicintiy: ", vicinity(mi, mj))
    print("This is the change of coord for the mj point with mi at (0,0,0)")
    print(vicinity(mi, mj))
