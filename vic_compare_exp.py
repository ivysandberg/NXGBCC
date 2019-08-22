''' This is an effort to identify a good vicinity hashing function.
It will also back-port and clean up some functions in V_F and exp_for.
'''

import argparse
from time import time
from multiprocessing import Process, Manager
import numpy as np
# from sklearn import datasets
# from sklearn.linear_model import LogisticRegression
# from sklearn.metrics import roc_auc_score
# from sklearn.model_selection import train_test_split

import scipy.stats

# for file/directory operations
from os import listdir
from os.path import isfile, join
import json

import csv

# for regular expressions
import re

# for matrix computations
#from Vectorizing_Fingerprints import *

#from experiment_for import *

import math
import time

import pickle
from pathlib import Path
import numpy as np
import os


# If hashing goes well, we should expect some things.
# 1) Hashes, individually, have high variance.
# 2) Within prints, hashes are "closer"
# 3) Between prints, hashes are "further"

def qdist(x,y):
    return math.sqrt(x **2 + y**2)


def mk_rep_set(num_vic_to_make, min_min_per_vic, max_min_per_vic, rad):
    rl = int(int(rad/2) *-1)
    rh = int(rad / 2)
    rep_set = []
    for i in range(num_vic_to_make):
        vic = [(0.0, 0.0, 0.0)]
        num_mins_this_vic = int(np.random.uniform(min_min_per_vic, max_min_per_vic))
        # generate thetas that are moderately correlated within vicinities
        # eye-ball guess from this figure:
        # https://sourceafis.machinezoo.com/template
        theta_min = np.random.uniform(0.01, 2*math.pi/3)
        theta_max = theta_min + np.random.uniform(.2, math.pi/3)
        #print(str(theta_min) + "," + str(theta_max))
        while len(vic) < num_mins_this_vic:
            x = int(np.random.uniform(rl, rh))
            y = int(np.random.uniform(rl,rh))
            if qdist(x,y) < rad:  # rejection sampling
                # Decide on angle:
                theta = np.random.uniform(theta_min, theta_max)
                if np.random.uniform() < 0.5:
                    theta = theta + math.pi
                this_vic = (x,y,theta)
                vic.append(this_vic)
        rep_set.append(vic)
    return rep_set

# usage:
# # rep_set = mk_rep_set(num_vic_to_make=3, min_min_per_vic=3, max_min_per_vic=5)

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

def new_mj_coords(mi, mj):
    X = compute_distance(mi, mj) * (math.cos(angle_between3(mi, mj) - mi[2]))
    Y = compute_distance(mi, mj) * (math.sin(angle_between3(mi, mj) - mi[2]))
    theta = mj[2] - mi[2]
    output = (X, Y, theta)
    return output

def unit_vector2(triple):
    vector_portion = (triple[0], triple[1])
    norm = np.linalg.norm(vector_portion, ord=2)
    first_entry = triple[0] / norm
    second_entry = triple[1] / norm
    third_entry = triple[2]
    return (first_entry, second_entry, third_entry)



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


# the following functions (distance, unit_vector, and angle_between points)
# are needed to compute the above 'change of coordinate' functions
def compute_distance(mi, mj):
    dist = math.sqrt((mj[0] - mi[0])**2 + (mj[1] - mi[1])**2)
    return dist

def parse_fp_files_to_list_of_list_of_minutiae(onlyfilesoftype, data_dir):
    indiv_minutia_coords = []
    for i in onlyfilesoftype:
        with open(join(data_dir, i)) as json_file:
            data = json.load(json_file)
        # print (data)
        json_string = json.dumps(data, indent=4)
        # print (json_string)
        minutia_of_one_print = create_minutiae_coords(data)
        indiv_minutia_coords.append(minutia_of_one_print)
        # print (len(minutia_of_one_print))
    return indiv_minutia_coords

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


# input: list of prints
def obtain_and_normalize_vicinities(list_of_prints, r):
    ### for specified radius, r, identify and normalize (i.e. rotate) vicinities
    normed_vics_by_fingerprint = []
    list_of_all_vicinities = []

    for i in range(len(list_of_prints)):
        if i % 100 == 0:
            print("Getting prints; have " + str(i))
        vicinities_in_print = get_vicinities_from_minutia_list(list_of_prints[i], r)
        normed_vics_in_print = normalize_vicinities(vicinities_in_print)
        normed_vics_by_fingerprint.append(normed_vics_in_print)

    # unlist the "by-fingerprint" part
    for i in normed_vics_by_fingerprint:
        for j in i:
            list_of_all_vicinities.append(j)

    return normed_vics_by_fingerprint, list_of_all_vicinities

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



def compare_minutiae(Mi, Mj):
    #, sigmaX, sigmaY, sigmatheta):
    # compute a pairing score between minutiae
    # score = 1 is two identical minutiae
    score = math.exp((-((Mi[0] - Mj[0])**2) / s1**2)) * math.exp(
        (-((Mi[1] - Mj[1])**2) / s1**2)) * math.exp(
            (-((Mi[2] - Mj[2])**2) / s2**2))
    output = np.array(score)
    return output

def invert_matrix(matrix):
    output = 1-matrix
    return output

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


def compare_vicinities_adjusted(vic1, vic2, type):
    lv1 = len(vic1) + 10
    lv2 = len(vic2) + 10
    percent_over = (np.max((lv1, lv2)) / np.min((lv1, lv2))) - 1
    #percent_over = (np.max((len(vic1), len(vic2))) / np.min((len(vic1), len(vic2)))) - 1
    pairing_scores_matrix = create_minutiae_pairing_scores_matrix(vic1, vic2)
    inv_matrix = invert_matrix(pairing_scores_matrix)
    if type == "greedy":
        score = greedy_min_cost(inv_matrix, False)
    else:
        score = lapsolver_get_comparison_score(inv_matrix)
    # choice: normalize by smaller vicinity size
    score = score / np.min((len(vic1), len(vic2)))
    delta = 1.0 - score
    xb = -3.0 + (4.0*percent_over)
    adj_factor = np.exp(xb) / (1 + np.exp(xb))
    adj_score = score + (delta * adj_factor)
    return adj_score

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

def create_print_hashed_vector2(list_of_normed_vics_from_print, rep_set):
    # print (len(minutia_of_one_print))
    #vics = get_vicinities_from_minutia_list(single_print, r)
    #norm_vics = normalize_vicinities(vics)
    # compare print_1 to the rep_set
    mat_compare_print_to_rep_set = create_vicinity_comparison_scores_matrix_full(list_of_normed_vics_from_print, rep_set, "greedy")
    # print (mat_print_to_rep_set)
    # choose the column-wise min = the vicinity closest to the representative vicinity
    print_vector = mat_compare_print_to_rep_set.min(axis=0)

    return print_vector

# Changed to put out the average minimum distance, rather than the average distance.
def get_dists(normed_alien, normed_within, rep_set):
    print("Making alien vectors")
    alien_vectors = []
    for counter, value in enumerate(normed_alien):
        vector = create_print_hashed_vector2(value, rep_set)
        alien_vectors.append(vector)
    print("Making within vectors")
    within_vectors = []
    for counter, value in enumerate(normed_within):
        vector = create_print_hashed_vector2(value, rep_set)
        within_vectors.append(vector)
    # save a random hashed vector
    with open("vectors.txt", 'a') as myfile:
        wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
        wr.writerow(vector)

    min_dists_within = []
    avg_dists_within = []
    max_dists_within = []
    for i in range(len(within_vectors)): # only within
        distances = []
        for j in range(len(within_vectors)):
            # note: don't compare i with i
            if i != j:
                dist = np.linalg.norm(within_vectors[i] - within_vectors[j])
                distances.append(dist)
        min_dists_within.append(np.min(distances))
        avg_dists_within.append(np.mean(distances))
        max_dists_within.append(np.max(distances))
    between_dists = []
    for i in range(len(within_vectors)):
        for j in range(len(alien_vectors)):
            dist = np.linalg.norm(within_vectors[i] - alien_vectors[j])
            between_dists.append(dist)
    mean_between = np.mean(between_dists)
    min_between = np.min(between_dists)
    # entropy of a hash:
    as_cat = [0,0,0,0,0,0,0,0,0,0]
    for k in range(9):
        as_cat[k] = len([i for i in vector if i > k*.1 and i < k*.1 + .1])
    entropy = scipy.stats.entropy(as_cat)


    output = (np.mean(min_dists_within), np.mean(avg_dists_within),
              np.mean(max_dists_within), mean_between, min_between, entropy)
    return output

if __name__ == "__main__":


    alien_dir = '/home/ubuntu/exp/alien'
    alien_files = [f for f in listdir(alien_dir) if isfile(join(alien_dir, f))]
    alien_templates = [f for f in alien_files if re.search("templatedata", f) != None]

    alien_minutia_from_file = parse_fp_files_to_list_of_list_of_minutiae(alien_templates, alien_dir)

    # for each in a range of parameters...
    for iter in range(1200):
        print("Iter" + str(iter))

        # sample a finger, get minutia
        identity_dirs = ['/home/ubuntu/exp/identity1', '/home/ubuntu/exp/identity2', '/home/ubuntu/exp/identity3',
                         '/home/ubuntu/exp/identity4', '/home/ubuntu/exp/identity5']

        hands = ["left", "right"]
        fingers = ["pinky", "index", "middle", "ring", "thumb"]

        id_dir = np.random.choice(identity_dirs)
        hand = np.random.choice(hands)
        finger = np.random.choice(fingers)

        id_files = id3_files = [f for f in listdir(id_dir) if isfile(join(id_dir, f))]
        id_templates = [f for f in id_files if re.search("templatedata", f) != None]
        id_hand = [f for f in id_templates if re.search(hand, f) != None]
        id_finger = [f for f in id_hand if re.search(finger, f) != None]

        for i in id_finger:
            print(i)

        within_minutia_from_file = parse_fp_files_to_list_of_list_of_minutiae(id_finger, id_dir)

        # set radii to choose from; choose one
        # this affects rep_set, and creating civs from minutia lists
        radii = [35, 40, 45, 50, 55, 60]
        r = np.random.choice(radii)

        # build representative set
        rep_set = mk_rep_set(num_vic_to_make=1000, min_min_per_vic=3, max_min_per_vic=25, rad=r)

        # build vicinities
        normed_alien_vics, all_alien_vics = obtain_and_normalize_vicinities(alien_minutia_from_file, r)
        normed_within_vics, all_within_vics = obtain_and_normalize_vicinities(within_minutia_from_file, r)


        global s1
        s1 = np.random.uniform(20, 160)
        global s2
        s2 = np.random.uniform(0.05, 3)
        # hash the vectors...calculate stats...
        stat = get_dists(normed_alien_vics, normed_within_vics, rep_set)

        line = str(iter) + ',' + str(stat[0]) + ',' + str(stat[1]) + ',' + \
               str(stat[2]) + ',' + str(stat[3]) + ',' + str(stat[4]) + ',' + str(stat[5]) + ',' + \
               str(s1) + ',' + str(s2) + ',' + str(r) + "\n"
        print(line)

        file1 = open("r_and_sigma_exp.txt", "a")  # append mode
        file1.write(line)
        file1.close()



