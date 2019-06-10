import numpy as np
import pandas as pd
import math


# list of all the minutiae points
list_of_minutiae = [(3,2,0.5), (-4,1.2,3.4), (5,-2,4.4), (1,-1,1), (3.6,2,0.5), (1,1,1), (-2,-0.3, 0), (1,2, -4), (0,0,0)]

''' Create a Vicinity of the minutiae points within a certain radius oriented around one central minutiae at (0,0,0) '''


# define vicinity - set central minutia at (0,0,0) and rotate/translate other minutia relative to the central
# general structure
#def vicinity (list_of_minutiae):
    #mi = list_of_minutiae[0]
    #mj = list_of_minutiae[1:]
    #output = []
    # find change of coords

# define central menutiae
mi = list_of_minutiae[0]
#print(mi[0])
mj = list_of_minutiae[:] - mi
#print(mj[1])

''' Find new X1 coord '''

# distance 2D
def distance(mj):
    dist = math.sqrt((mj[0] - mi[0])**2 + (mj[1] - mi[1])**2)
    return dist
print("Disance between mi and mj1: ", distance(mj[1]))

# distance 3D (not sure which to use)
def distance_3D(mj):
    dist = math.sqrt((mj[0] - mi[0])**2 + (mj[1] - mi[1])**2 + (mj[2] - mi[2])**2)
    return dist
#print(distance_3D(mj[1]))

# define cos(alpha(mi,mj))
def unit_vector(list_of_minutiae):
    return list_of_minutiae / np.linalg.norm(list_of_minutiae)  # this statement breaks when mj = (0,0,0)
#print(unit_vector(mi))

def angle_between(mi, mj):
    mi_u = unit_vector(mi)
    mj_u = unit_vector(mj)
    return np.arccos(np.clip(np.dot(mi_u, mj_u), -1.0, 1.0))
    #return np.arccos(np.dot(mi,mj) / (math.sqrt(mi[0]**2 + mi[1]**2 + mi[2]**2) * math.sqrt(mj[0]**2 + mj[1]**2 + mj[2]**2)))
    # ^^ both of these functions for angle between two vectors returns the same output, but the second breaks on mj = (0,0,0)
print("Angle between mi and mj1: ", angle_between(mi,mj[1]))

# define X coordinate of Vicinity 1 (corresponds to mi and mj[1])
x_1 = distance(mj[1]) * (math.cos(angle_between(mi,mj[1]) - mi[2]))
print("X1: ", x_1)


''' Find new Y1 coord '''
y_1 = distance(mj[1]) * (math.sin(angle_between(mi,mj[1])-mi[2]))
print("Y1: ", y_1)

''' Find new theta1 '''
def theta(mj):
    new_theta = mi[2] - mj[2]
    return new_theta
theta_1 = theta(mj[1])
print("Theta1: ", theta_1)

def vicinity_1(mj):
    output = []
    output.append([x_1, y_1, theta_1])
    return output
print("Vicintiy_1: ", vicinity_1(mj[1]))


''' Generalize to all mj'''
for j in mj:
    #print("Distance between mi and mj: ", distance(j))
    #print("Angle between mi and mj: ", angle_between(mi,j))
    x = distance(j) * (math.cos(angle_between(mi,j) - mi[2]))
    y = distance(j) * (math.sin(angle_between(mi,j)-mi[2]))
    #print("x coord", x)
    #print("y coord", y)
    def theta(j):
        new_theta = mi[2] - j[2]
        return new_theta
    def vicinity(j):
        output = []
        output.append([x,y,theta(j)])
        return output
    print("Vicinity: ", vicinity(j))
