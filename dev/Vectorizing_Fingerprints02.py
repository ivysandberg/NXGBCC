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
#import pandas as pd
import math
import itertools
from munkres import Munkres

from lapsolver import solve_dense


def new_mj_coords(mi, mj):
    X = compute_distance(mi, mj) * (math.cos(angle_between3(mi, mj) - mi[2]))
    Y = compute_distance(mi, mj) * (math.sin(angle_between3(mi, mj) - mi[2]))
    theta = mj[2] - mi[2]
    output = (X, Y, theta)
    return output


# the following functions (distance, unit_vector, and angle_between points)
# are needed to compute the above 'change of coordinate' functions
def compute_distance(mi, mj):
    dist = math.sqrt((mj[0] - mi[0])**2 + (mj[1] - mi[1])**2)
    return dist


# unit_vector function to compute unit vectors which are needed to compute the angles between vectors
# Input: any number of minutiae
# Output: the minutiae in unit vector form


###############################
### edited unit_norm function; only "norm" of the vector portion, not theta

# TODO: probably don't want norm of the triple; likely just the x,y portion...will test
def unit_vector2(triple):
    vector_portion = (triple[0], triple[1])
    norm = np.linalg.norm(vector_portion, ord=2)
    first_entry = triple[0] / norm
    second_entry = triple[1] / norm
    third_entry = triple[2]
    return (first_entry, second_entry, third_entry)





###############################

# TODO this.
# angle_between3 returns only positive angles. Is this a problem?
# If not, need to handle other cases
def angle_between3(mi, mj):
    # translate mj by mi
    point = (mj[0] - mi[0], mj[1]-mi[1], mj[2])
    # a point on the x-axis, which already is unit vector length
    axis_point = (1,0)
    mj_u = unit_vector2(point)
    mj_u_vec = (mj_u[0], mj_u[1])
    angle = np.arccos(np.clip(np.dot(axis_point, mj_u_vec), -1.0, 1.0))
    if mj_u[1] < 0:
        angle = (2*(math.pi - angle)) + angle
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
    all_mjs = vicinity[1:]
    normalized_vicinity = []
    # add in the mi:
    normalized_mi_triple = (0.0, 0.0, 0.0)
    normalized_vicinity.append(normalized_mi_triple)
    # loop over all possible mjs
    for mj in all_mjs:
        try:
            new_mj = new_mj_coords(mi, mj)
        except RuntimeWarning:
            print(vicinity)
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

#sigmaX = math.sqrt(21)
#sigmaY = math.sqrt(15)
#sigmatheta = math.sqrt(0.037)

# Conduct comparisons to find similar minutiae
# This function compares minutia points
# note: hardcoding some findings for sigmaX, sigmaY, and sigmatheta
'''
def compare_minutiae(Mi, Mj):
    #, sigmaX, sigmaY, sigmatheta):
    # compute a pairing score between minutiae
    # score = 1 is two identical minutiae
    score = math.exp((-((Mi[0] - Mj[0])**2) / sigmaX**2)) * math.exp(
        (-((Mi[1] - Mj[1])**2) / sigmaY**2)) * math.exp(
            (-((Mi[2] - Mj[2])**2) / sigmatheta**2))
    output = np.array(score)
    return output
'''
def compare_minutiae(Mi, Mj):
    #, sigmaX, sigmaY, sigmatheta):
    # compute a pairing score between minutiae
    # score = 1 is two identical minutiae
    score = math.exp((-((Mi[0] - Mj[0])**2) / s1**2)) * math.exp(
        (-((Mi[1] - Mj[1])**2) / s1**2)) * math.exp(
            (-((Mi[2] - Mj[2])**2) / s2**2))
    output = np.array(score)
    return output
#m1 = (3, 8, 0.5)
#m2 = (3, 9, 0.45)
#compare_minutiae(m1, m2)

# Input two vicinities
# Output Constructs a matrix populated with the pairing scores of each minutiae in vic1 matched with minutiae in vic2
def create_minutiae_pairing_scores_matrix(vic_i, vic_j):
    output_matrix = np.zeros((len(vic_i), len(vic_j)), dtype=np.float32)
    for i in range(len(vic_i)):
        for j in range(len(vic_j)):
            #matrix_value = compare_minutiae(vic_i[i], vic_j[j], 135, 120, 1.85)
            matrix_value = compare_minutiae(vic_i[i], vic_j[j])
            '''
            Using a sample standard deviations of each x, y, & theta from one fingerprint image.
            Need more accurate measurements

            '''
            output_matrix[i, j] = matrix_value
    return output_matrix



# Compare Vicinities

# Create matrix filled with values = 1 - pairing score, in order to find maximums
# Do this b/c the Hungarian alg is designed to minimize cost (in the assignment problem)
'''
def invert_matrix(matrix):
    inv_matrix = []
    for row in matrix:
        #print(row)
        inv_row = []
        for col in row:
            inv_row += [1 - col]
        inv_matrix += [inv_row]
    return inv_matrix
'''

def invert_matrix(matrix):
    output = 1-matrix
    return output


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
        # print (row, col)
        val = inv_matrix[row][col]
        total += val
    comparison_score = total        # use this to print out just the comparison score and then populate the comparison score matrix
    # comparison_score = f'sum of scores: {total}'      # this prints out "Sum of Scores: total"
    return comparison_score

# costs = np.array([[0.6,0.9,0.1], [0.10,0.3,0.2], [0.8,0.7,0.4]], dtype=np.float32)

def lapsolver_get_comparison_score(inv_matrix):
    rows, cols = solve_dense(inv_matrix)
    total = 0
    for row, col in zip(rows, cols):
        # print(row, col) # Row/column pairings
        val = inv_matrix[row][col]
        total += val
    comparison_score = total        # use this to print out just the comparison score and then populate the comparison score matrix
    # comparison_score = f'sum of scores: {total}'      # this prints out "Sum of Scores: total"
    return comparison_score



''' Greedy algorithm for solving the problem.
Input is 2d numpy array of costs, all > 0.
Output is total cost of GREEDY min assignment.
Assumes costs are >= zero
'''

import copy

def greedy_min_cost(cost_matrix, print_matches):
    cost_matrix = np.array(cost_matrix)
    # make sure there are fewer rows than cols
    if cost_matrix.shape[0] > cost_matrix.shape[1]:
        mat = copy.copy(cost_matrix.transpose())
    else:
        mat = copy.copy(cost_matrix)
    # make a variable that is over the max of the matrix (can never be a row-wise min)
    over_max = np.max(mat) + 1
    # initialize cost to zero
    cost = 0
    # for each row...
    for i in range(mat.shape[0]):
        # find the index of the row that minimizes the cost of assigning row i to col j
        row_min_ind = np.argmin(mat[i])
        # show which cells were the matches, if print_matches=True
        # if print_matches:
            # print((i, row_min_ind))
        # What was the min of that row
        row_min = mat[i,row_min_ind]
        # add that to costs
        cost += row_min
        # set the column that was used/matched to over_max so it is never a min again
        mat[:,row_min_ind] = over_max
    return cost





def compare_vicinities(vic1, vic2, type):
    pairing_scores_matrix = create_minutiae_pairing_scores_matrix(vic1, vic2)
    inv_matrix = invert_matrix(pairing_scores_matrix)
    if type == "greedy":
        score = greedy_min_cost(inv_matrix, False)
    else:
        score = lapsolver_get_comparison_score(inv_matrix)
    return score

def compare_vicinities_adjusted(vic1, vic2, type):
    percent_over = (np.max((len(vic1), len(vic2))) / np.min((len(vic1), len(vic2)))) - 1
    pairing_scores_matrix = create_minutiae_pairing_scores_matrix(vic1, vic2)
    inv_matrix = invert_matrix(pairing_scores_matrix)
    if type == "greedy":
        score = greedy_min_cost(inv_matrix, False)
    else:
        score = lapsolver_get_comparison_score(inv_matrix)
    # choice: normalize by smaller vicinity size
    score = score / np.min((len(vic1), len(vic2)))
    delta = 1.0 - score
    xb = -3.0 + (13.0*percent_over)
    adj_factor = np.exp(xb) / (1 + np.exp(xb))
    adj_score = score + (delta * adj_factor)
    return adj_score


def compare_vicinities_adjusted_debug(vic1, vic2, type):
    percent_over = (np.max((len(vic1), len(vic2))) / np.min((len(vic1), len(vic2)))) - 1
    print(percent_over)
    pairing_scores_matrix = create_minutiae_pairing_scores_matrix(vic1, vic2)
    print(pairing_scores_matrix)
    inv_matrix = invert_matrix(pairing_scores_matrix)
    print(inv_matrix)
    if type == "greedy":
        score = greedy_min_cost(inv_matrix, False)
    else:
        score = lapsolver_get_comparison_score(inv_matrix)
    print("score: " + str(score))
    # choice: normalize by smaller vicinity size
    score = score / np.min((len(vic1), len(vic2)))
    delta = 1.0 - score
    xb = -3.0 + (13.0*percent_over)
    adj_factor = np.exp(xb) / (1 + np.exp(xb))
    adj_score = score + (delta * adj_factor)
    print(delta)
    print(xb)
    print(adj_factor)
    print(adj_score)
    return adj_score

def quick_compare_vicinities_adjusted(vic1, vic2, type):
    #percent_over = (np.max((len(vic1), len(vic2))) / np.min((len(vic1), len(vic2)))) - 1
    lv1 = len(vic1) + 10
    lv2 = len(vic2) + 10
    percent_over = (np.max((lv1, lv2)) / np.min((lv1, lv2))) - 1
    if percent_over > 0.5:
        adj_score = 0.98
        #print("Savings!")
    else:
        #print("None")
        pairing_scores_matrix = create_minutiae_pairing_scores_matrix(vic1, vic2)
        inv_matrix = invert_matrix(pairing_scores_matrix)
        if type == "greedy":
            score = greedy_min_cost(inv_matrix, False)
        else:
            score = lapsolver_get_comparison_score(inv_matrix)
        # choice: normalize by smaller vicinity size
        score = score / np.min((len(vic1), len(vic2)))
        delta = 1.0 - score
        xb = -3.0 + (13.0*percent_over)
        adj_factor = np.exp(xb) / (1 + np.exp(xb))
        adj_score = score + (delta * adj_factor)
    return adj_score


from os import listdir
from os.path import isfile, join
import json
import re
#my_data_directory = "/Users/ivysandberg/dmc/NXGBCC/fingerprint_hashing/data/fingerprints"
# get all the files (not folders) in 'my_data_directory'
#onlyfiles = [f for f in listdir(my_data_directory) if isfile(join(my_data_directory, f))]

# select down to the files we want with regular expressions
#onlyfilesoftype = [f for f in onlyfiles if re.search("templatedata", f) != None]

def create_minutiae_coords(data):
    list_of_minutiae_coords = []
    minutiae = data['minutiae']
    for i in minutiae:
        # print (i)
        output = (i['x'], i['y'], i['direction'])
        # print (output)
        list_of_minutiae_coords.append(output)
    return list_of_minutiae_coords

def minutia_from_one_print(onlyfilesoftype):
    indiv_minutia_coords = []
    for i in onlyfilesoftype:
        with open(join(my_data_directory, i)) as json_file:
            data = json.load(json_file)
        # print (data)
        json_string = json.dumps(data, indent=4)
        # print (json_string)
        minutia_of_one_print = create_minutiae_coords(data)
        indiv_minutia_coords.append(minutia_of_one_print)
        # print (len(minutia_of_one_print))
    return indiv_minutia_coords

#indiv_prints = minutia_from_one_print(onlyfilesoftype)

#list_of_all_vicinities = []
#for i in indiv_prints:
#        vicinities = get_vicinities_from_minutia_list(i, 50)
#        for vic in vicinities:
#            list_of_all_vicinities.append(vic)


def create_vicinity_comparison_scores_matrix(list_of_vicinities_1, list_of_vicinities_2, type):
    #print("len lov1 = " + str(len(list_of_vicinities_1)))
    #print("len lov2 = " + str(len(list_of_vicinities_2)))
    output_matrix = np.zeros((len(list_of_vicinities_1), len(list_of_vicinities_2)), dtype=np.float32)
    for i in range(len(list_of_vicinities_1)):
        if i % 500 == 0:
            print("Computed " + str(i) + " vicinity matrix rows")
        for j in range(len(list_of_vicinities_2)):
            if i > j:
                matrix_value = quick_compare_vicinities_adjusted(list_of_vicinities_1[i], list_of_vicinities_2[j], type)
                output_matrix[i, j] = matrix_value
                output_matrix[j, i] = matrix_value
    return output_matrix

# this variant (_full) does not have the optimization that assumes a square matrix, because it won't always be
def create_vicinity_comparison_scores_matrix_full(list_of_vicinities_1, list_of_vicinities_2, type):
    #print("len lov1 = " + str(len(list_of_vicinities_1)))
    #print("len lov2 = " + str(len(list_of_vicinities_2)))
    output_matrix = np.zeros((len(list_of_vicinities_1), len(list_of_vicinities_2)), dtype=np.float32)
    for i in range(len(list_of_vicinities_1)):
        for j in range(len(list_of_vicinities_2)):
            # if i == 3 and j == 3:
            #     matrix_value = compare_vicinities_adjusted_debug(list_of_vicinities_1[i], list_of_vicinities_2[j])
            # else:
            matrix_value = compare_vicinities_adjusted(list_of_vicinities_1[i], list_of_vicinities_2[j], type)
            output_matrix[i, j] = matrix_value
    return output_matrix

'''
# Itertools variant
def create_vicinity_comparison_scores_matrix(vicinity_list):
    output_matrix = np.zeros((len(vicinity_list), len(vicinity_list)))
    
    for pair in itertools.combinations(range(len(vicinity_list)), r=2):
        i, j = pair
        output_matrix[i, j] = compare_vicinities(vicinity_list[i], vicinity_list[j])
    
    i_lower = np.tril_indices(output_matrix.shape[0], -1)
    output_matrix[i_lower] = output_matrix.T[i_lower]
    return output_matrix
'''

# if two vicinities have a comparison_score < the min_comparison_score then one must be removed
# this threshold approach is not being used 
def filter_list_of_vinities_by_similarity(list1, list2, min_comparison_score):
    output = []
    for i in list1:
        for j in list2:
            if compare_vicinities(i,j) > min_comparison_score:
                output.append(i)
    return output



# FIND K VICINITIES THAT MAXIMIZE THE MINIMUM DISTANCE BETWEEN X AND ALL POINTS IN THE HEAP

# print ("Matrix of Vicinity Comparison Scores: ")
# print (mat)

# Input: matrix as a 2D numpy array and k = # of vicinities want in Representative set
# Ouput: indices of the vicinities for the Representative set



# mat = np.random.rand(6,6)
# np.fill_diagonal(mat, 0)



# print (k_distinct_vicinities(mat, 4))
    
def optimized_k_distinct_vicinities(mat, k):
    heap = []
    # supposed we have some points already, called 'heap'


    # function
    r = np.random.choice(range(mat.shape[0]))  # randomly choose the first row
    heap.append(r)
    next_vic = np.argmax(mat[r])  # choose the row furthest from that one
    heap.append(next_vic)

    for j in range(k - 2):
        # make a mat of the rows in heap, all columns
        nmat = mat[heap,]
        # find min distance from the heap to each point
        # note, this will be zero for points in heap, assuming points have zero dist to themselves.
        # CHECK THIS ASSUMPTION WITH IVY!
        mins_heap_to_points = np.amin(nmat, axis=0)
        # Find the index of the point with the largest minimum distance
        index_of_max_of_mins = np.argmax(mins_heap_to_points)
        heap.append(index_of_max_of_mins)
   # note: if choose a k > the # of vicinties in mat then it will just append 0s
    print("Heap")
    print(sorted(heap))

    return heap

# print (optimized_k_distinct_vicinities(mat, 5))


# TEST

import json


#with open('/Users/ivysandberg/dmc/NXGBCC/fingerprint_hashing/data/fingerprints/f0001_01.png.templatedata') as json_file:
#    data = json.load(json_file)

def create_minutiae_coords(data):
    list_of_minutiae_coords = []
    minutiae = data['minutiae']
    for i in minutiae:
        # print (i)
        output = (i['x'], i['y'], i['direction'])
        # print (output)
        list_of_minutiae_coords.append(output)
    return list_of_minutiae_coords



#minutia_list = create_minutiae_coords(data)
# print (minutia_list)

