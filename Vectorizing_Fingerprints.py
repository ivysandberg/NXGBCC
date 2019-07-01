"""
We are choosing to characterize Fingerprints by neighborhoods of defining points.

These defining points are called Minuitae points and are tuples of length four: (X, Y, direction in radians, type).
Type is either a ridge ending or a bifrication.

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
    Y = -compute_distance(mi, mj) * (math.sin(angle_between(mi, mj) - mi[2]))
    return Y

def new_theta_coord(mi, mj):
    theta = mj[2] - mi[2]
    return theta

def new_mj_coords(mi, mj):
    output = (new_x_coord(mi,mj), new_y_coord(mi,mj), new_theta_coord(mi,mj))
    return output

# the following functions (distance, unit_vector, and angle_between points)
# are needed to compute the above 'change of coordinate' functions
def compute_distance(mi, mj):
    dist = math.sqrt((mj[0] - mi[0])**2 + (mj[1] - mi[1])**2)
    return dist


# function to compute unit vectors which are needed to compute the angles between vectors
# Input: any number of minutiae
# Output: the minutiae in unit vector form
def unit_vector(list_of_minutiae):
    return list_of_minutiae / np.linalg.norm(
        list_of_minutiae)  # this statement breaks when minutiae = (0,0,0)


def angle_between(mi, mj):
    mi_u = unit_vector(mi)
    mj_u = unit_vector(mj)
    angle = np.arccos(np.clip(np.dot(mi_u, mj_u), -1.0, 1.0))
    return angle


# THE FOLLOWING FUNCTION CONSTRUTS VICINITIES (WHICH ARE LISTS OF MINUTIAE)
# each vicinity covers a certain area of a fingerprint & areas overlap

# Input: this function takes in a minutia list and a given radius
# Output: it returns all the un-normalizedvicinities it can find
# by taking each minutia in the list and returning all the minutiae within a given radius.
# the un-normalized vicinities are a list of triples
def get_vicinities_from_minutia_list(minutia_list, radius):
    all_vicinities = []
    for i in minutia_list:
        this_vicinity = []
        this_vicinity.append(i)
        for j in minutia_list:
            if compute_distance(i,j) < radius and i != j: #add if i != j
                this_vicinity.append(j)
        all_vicinities.append(this_vicinity)
    return all_vicinities


# Function takes in a list of vicinities (the output from get_vicinities_from_minutia_list)
# Returns the same list, but removes vicinities with too few or too many minutia points
def filter_list_of_vinities_by_size(list_of_vicinities, min_len, max_len):
    output = []
    for i in list_of_vicinities:
        if len(i) >= min_len and len(i) <= max_len:
            output.append(i)
    return output

# Input: a single vicinity
# (i.e. a single vicinity indexed from the list from the output from get_vicinities_from_minutia_list or filter_list_of_vinities_by_size)
# Output: a normalized_vicinity (the coordinates of the mjs have been changed relative to the set mi when moved to (0,0,0))
def normalize_this_vicinity(vicinity):
    mi = vicinity[0]
    mj = vicinity[1:]
    normalized_vicinity = []
    for mj in vicinity:
        new_mj = new_mj_coords(mi, mj)
        normalized_vicinity.append(new_mj)
    return normalized_vicinity


# Input: a list of vicinities (i.e. the output from get_vicinities_from_minutia_list or filter_list_of_vinities_by_size)
# Ouput: all the vicinities in the list get normalized
def normalize_vicinities(list_of_vicinities):
    output = []
    for i in list_of_vicinities:
        output.append(normalize_this_vicinity(i))
    return output




# Compare Minutiae

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
def create_minutiae_pairing_scores_matrix(vic_i, vic_j):
    output_matrix = np.zeros((len(vic_i), len(vic_j)))
    for i in range(len(vic_i)):
        for j in range(len(vic_j)):
            matrix_value = compare_minutiae(vic_i[i], vic_j[j], 135, 120, 1.85)
            '''
            Using a sample standard deviations of each x, y, & theta from one fingerprint image.
            Need more accurate measurements

            '''
            output_matrix[i, j] = matrix_value
    return output_matrix



# Compare Vicinities

# Create matrix filled with values = 1 - pairing score, in order to find maximums
# Do this b/c the Hungarian alg is designed to minimize cost (in the assignment problem)
def invert_matrix(matrix):
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
def get_comparison_score(inv_matrix):
    m = Munkres()  # Munkres = Hungarian
    indexes = m.compute(inv_matrix)
    #print_matrix(inv_matrix, msg='compute highest pairing scores')
    total = 0
    for row, col in indexes:
        val = inv_matrix[row][col]
        total += val
    comparison_score = total        # use this to print out just the comparison score and then populate the comparison score matrix
    # comparison_score = f'sum of scores: {total}'      # this prints out "Sum of Scores: total"
    return comparison_score


def compare_vicinities(vic1, vic2):
    pairing_scores_matrix = create_minutiae_pairing_scores_matrix(vic1, vic2)
    inv_matrix = invert_matrix(pairing_scores_matrix)
    score = get_comparison_score(inv_matrix)
    return score



def create_vicinity_comparison_scores_matrix(list_of_vicinities_1, list_of_vicinities_2):
    output_matrix = np.zeros((len(list_of_vicinities_1), len(list_of_vicinities_2)))
    for i in range(len(list_of_vicinities_1)):
        for j in range(len(list_of_vicinities_2)):
            matrix_value = compare_vicinities(list_of_vicinities_1[i], list_of_vicinities_2[j])
            output_matrix[i, j] = matrix_value
    return output_matrix



def filter_list_of_vinities_by_similarity(list1, list2, min_comparison_score):
    # if two vicinities have a comparison_score < the min_comparison_score then one must be removed
    output = []
    for i in list1:
        for j in list2:
            if compare_vicinities(i,j) > min_comparison_score:
                output.append(i)
    return output






# TEST FUNCTIONS :

# example set of minutiae
list_of_minutiae = [(3, 2, 0.5), (-4, 1.2, 3.4), (5, -2, 4.4), (1, -1, 1),
                    (3.6, 2, 0.5), (1, 1, 1), (-2, -0.3, 0), (1, 2, -4)]



# Create a change of coords for every combination of mi (center minutiae) and mj (other minutiae)

#for mi, mj in itertools.combinations_with_replacement(list_of_minutiae, 2):
    # print("This is the (mi, mj) coordinate pair being changed: ", mi, mj)
    # print ("New mj coords with mi now at (0,0,0): ", new_mj_coords(mi,mj))
    # print ("Pairing Score: ", compare_minutiae(mi, mj, 2.5, 4.1, 1.2))



# Create a list of of vicinities
# each vicinity has a list of minutia points within a certain radius from mi
# loop through list_of_minutiae so each minutia is mi once
# print ("List of Vicinities: ")
# print (get_vicinities_from_minutia_list(list_of_minutiae, 3))


# Normalize vicinities
# sets the first minutiae as mi = (0,0,0) and translates the other minutiae (mjs) accordingly
list_of_vicinities = get_vicinities_from_minutia_list(list_of_minutiae,3)
# print ("Normalized List of Vicinities: ")
# print (normalize_vicinities(list_of_vicinities))


# Compare two vicinities
# Create a matrix of pairing scores
normalized_list_of_vicinities = normalize_vicinities(list_of_vicinities)
vic1 = normalized_list_of_vicinities[0]
vic2 = normalized_list_of_vicinities[1]
# print("Pairing Scores Matrix between Vic1 and Vic2: ")
# print (create_minutiae_pairing_scores_matrix(vic1, vic2))

# Compute Pair score
# print ("Comparison Score of Vic1 and Vic2: ", compare_vicinities(vic1, vic2))



# Output the comparison scores for every combination of vicinities in a list of vicinities :
# for vic_i, vic_j in itertools.combinations_with_replacement(normalized_list_of_vicinities, 2):
#     print (vic_i, vic_j)
#     print (compare_vicinities(vic_i, vic_j))


# Build a matrix populated with vicinity comparison scores
# in this case I'm building it out by (normalized_list_of_vicinities X normalized_list_of_vicinities)
list = filter_list_of_vinities_by_size(normalized_list_of_vicinities, 2, 4)
# print (list)
mat = create_vicinity_comparison_scores_matrix(list, list)
# print (type(mat))





# FIND K VICINITIES THAT MAXIMIZE THE MINIMUM DISTANCE BETWEEN X AND ALL POINTS IN THE HEAP

# print ("Matrix of Vicinity Comparison Scores: ")
# print (mat)

# Input: matrix as a 2D numpy array and k = # of vicinities want in Representative set
# Ouput: indices of the vicinities for the Representative set

def k_distinct_vicinities(mat, k):

    heap = []
    r = np.random.choice(range(mat.shape[0]))       # randomly choose the first row
    heap.append(r)
    next_vic = np.argmax(mat[r])    # choose the row furthest from that one
    heap.append(next_vic)

    # supposed we have some points already, called 'heap'


    for j in range(k-2):
        min_distances_to_the_heap = []
        best_point = 'done'

        print ("heap: ", heap)

        for i in range(mat.shape[1]):
            if i in heap:
                continue
            this_point_min_dist_to_heap = np.min(mat[heap,i])
            print (mat[heap,i])
            min_distances_to_the_heap.append(this_point_min_dist_to_heap)
            current_largest_min = np.max(min_distances_to_the_heap)
            if this_point_min_dist_to_heap == current_largest_min:
                best_point = i

        print("min distances: ", min_distances_to_the_heap)

        heap.append(best_point)

    chosen_vicinities_from_heap = []
    for i in heap:
        chosen_vicinities = list[i]
        chosen_vicinities_from_heap.append(chosen_vicinities)
        # print ("chosen_vicinities: ", chosen_vicinities)

    return heap, ("chosen_vicinities: ", chosen_vicinities_from_heap)


# print (k_distinct_vicinities(mat, 4))
