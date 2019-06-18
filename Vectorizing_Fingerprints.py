"""
We are choosing to characterize Fingerprints by neighborhoods of defining points.

These defining points are called Minuitae points and are tuples of length four: (X, Y, direction in radians, type).
Type is either a ridge ending or a bifrications.

The neighborhoods of these points are called vicinities. Multiple vicinities cover the space of a fingerprint.
Vicinities are represented as vectors, and a set of vicinities from one fingerprint are compared to the set of another to determine matches.

In order to define a vicinity, one sets one minutiae point as the vicinity's origin,
and orients the other minutiae accordingly that are within a specified radius from the origin.



"""

import numpy as np
import pandas as pd
import math
import itertools
from munkres import Munkres, print_matrix


# functions to compute the change of coordinates for a pair of Minutiae (mi, mj)
# new coordinates = [d(mi,mj)cos(ami,mj), -d(mi,mj)sin(ami,mj), Bmi,mj] computes the new X, Y, and θ.
# ami,mj = θmi,mj - θmi
# Bmi,mj = θmj - θmi


# Input: a pair of minutiae (mi, mj)
# Output: the new coordinates for mj when it's changed relative to mi changed to (0,0,0)
def new_x_coord(mi, mj):
    X = compute_distance(mi, mj) * (math.cos(angle_between(mi, mj) - mi[2]))
    return X

def new_y_coord(mi, mj):
    Y = compute_distance(mi, mj) * (math.sin(angle_between(mi, mj) - mi[2]))
    return Y

def new_theta_coord(mi, mj):
    θ = mi[2] - mj[2]
    return θ

def new_mj_coords(mi, mj):
    output = []
    output.append([new_x_coord(mi,mj), new_y_coord(mi,mj), new_theta_coord(mi,mj)])  # creates the double brackets
    return output

# the following functions (distance, unit_vector, and angle_between points)
# are needed to compute the above change of coordinates functions
def compute_distance(mi, mj):
    dist = math.sqrt((mj[0] - mi[0])**2 + (mj[1] - mi[1])**2)
    return dist

# function to compute unit vectors which are needed to compute the angles between vectors
# Input: any number of minutiae
# Output: the minutiae in unit vector form
def unit_vector(list_of_minutiae):
    return list_of_minutiae / np.linalg.norm(
        list_of_minutiae)  # this statement breaks when mj = (0,0,0)


def angle_between(mi, mj):
    mi_u = unit_vector(mi)
    mj_u = unit_vector(mj)
    angle = np.arccos(np.clip(np.dot(mi_u, mj_u), -1.0, 1.0))
    return angle



# THE FOLLOWING FUNCTION CONSTRUTS VICINITIES (LISTS OF MINUTIAE)
# each vicinity covers a certain area of a fingerprint & areas overlap

# this function takes in a minutia list and radius
# it returns all the un-normalizedvicinities it can find
# the un-normalized vicinities are a list of triples
def get_vicinities_from_minutia_list(minutia_list, radius):
    all_vicinities = []
    for i in minutia_list:
        this_vicinity = []
        this_vicinity.append(i)
        for j in minutia_list:
            if compute_distance(i,j) < radius:
                this_vicinity.append(j)
        all_vicinities.append(this_vicinity)
    return all_vicinities


# Function takes in a list of vicinities (as put out by getVicinitiesFromMinutiaList)
# Returns the same list, but removes vicinities with too few or too many minutia points
def filter_list_of_vinities_by_size(list_of_vicinities, min_len, max_len):
    output = []
    for i in list_of_vicinities:
        if len(i) > min_len and len(i) < max_len:
            output.append(i)
    return output

def normalize_this_vicinity(vicinity):
    mi = vicinity[0]
    mj = vicinity[1:]
    normalized_vicinity = []
    for mj in vicinity:
        new_mj = new_mj_coords(mi, mj)
        normalized_vicinity.append(new_mj)
    return normalized_vicinity

def normalize_vicinities(list_of_vicinities):
    output = []
    for i in list_of_vicinities:
        output.append(normalize_this_vicinity(i)) # note: you'll need to write 'normalizeThisVicinity'
    return output




# Conduct comparisons to find similar minutiae
# This function compares minutia points
def compare_minutiae(Mi, Mj, sigmaX, sigmaY, sigmatheta):
    # compute a pairing score between minutiae
    # score = 1 is two identical minutiae
    score = math.exp((-((Mi[0] - Mj[0])**2) / sigmaX**2)) * math.exp(
        (-((Mi[1] - Mj[1])**2) / sigmaY**2)) * math.exp(
            (-((Mi[2] - Mj[2])**2) / sigmatheta**2))
    output = np.array(score)
    return output


# Input two vicinities
# Output Constructs a matrix populated with the pairing scores of each minutiae in vic1 matched with minutiae in vic2
def minutiae_pairing_scores_matrix(vic_i, vic_j):
    output_matrix = np.zeros((len(vic_i), len(vic_j)))
    for i in range(len(vic_i)):
        for j in range(len(vic_j)):
            matrix_value = compare_minutiae(vic_i[i], vic_j[j], 2.5, 4.1, 1.2)
            output_matrix[i, j] = matrix_value
    return output_matrix



# Compare Vicinities

# Create matrix filled with values = 1 - pairing score, in order to find maximums
# Hungarian is designed to minimize cost (in the assignment problem)
def Hungarian(matrix):
    inv_matrix = []
    for row in matrix:
        #print(row)
        inv_row = []
        for col in row:
            inv_row += [1 - col]
            inv_matrix += [inv_row]
    return inv_matrix


# Derive final comparison score from the vicinities wanting to compare

# Input: the inverted minutiae_pairing_scores_matrix of two vicinities
# Inputed matrix computed by Hungarian(output_matrix from the minutiae_pairing_scores_matrix)

# Output: use the Hungarian alg to compute the best pairs in the two vicinities
# and an overall pair score of the two vicinities (identical = 0)
def compare_vicinities(inv_matrix):
    m = Munkres()  # Munkres = Hungarian
    indexes = m.compute(inv_matrix)
    #print_matrix(inv_matrix, msg='compute highest pairing scores')
    total = 0
    for row, col in indexes:
        val = inv_matrix[row][col]
        total += val
    comparison_score = f'sum of scores: {total}'
    return comparison_score





# TEST FUNCTIONS :

# example set of minutiae
list_of_minutiae = [(3, 2, 0.5), (-4, 1.2, 3.4), (5, -2, 4.4), (1, -1, 1),
                    (3.6, 2, 0.5), (1, 1, 1), (-2, -0.3, 0), (1, 2, -4)]



# Create a change of coords for every combination of mi (center minutiae) and mj (other minutiae)

# for mi, mj in itertools.combinations_with_replacement(list_of_minutiae, 2):
#     print("This is the (mi, mj) coordinate pair being changed: ", mi, mj)
#     print ("New mj coords with mi now at (0,0,0): ", new_mj_coords(mi,mj))
#     print ("Pairing Score: ", compare_minutiae(mi, mj, 2.5, 4.1, 1.2))
#

# Create a list of of vicinities
# each vicinity has a list of minutia points within a certain radius from mi
# loop through list_of_minutiae so each minutia is mi once
# print (get_vicinities_from_minutia_list(list_of_minutiae, 3))

# filter list
vicinities_with_radius_of_3 = get_vicinities_from_minutia_list(list_of_minutiae, 3)
# print (filter_list_of_vinities_by_size(vicinities_with_radius_of_3, 4, 10))

# Normalize the first vicinity in the list of vicinities
vicinity_1 = vicinities_with_radius_of_3[0]
# print (normalize_this_vicinity(vicinity_1))

# Normalize all the vicinities in a list of vicinities
# print (normalize_vicinities(vicinities_with_radius_of_3))


# populate a matrix with the pairing scores of all the minutiae in on vicinity to all the minutiae in another vicinity
vicinity_2 = vicinities_with_radius_of_3[1]
print ("Vicinity 1: ", vicinity_1)
print ("Vicinity 2: ", vicinity_2)
print (minutiae_pairing_scores_matrix(vicinity_1, vicinity_2))

output_matrix = minutiae_pairing_scores_matrix(vicinity_1, vicinity_2)
# print (Hungarian(output_matrix))

inv_matrix = Hungarian(output_matrix)
print ("Comparison Score: ", compare_vicinities(inv_matrix))
