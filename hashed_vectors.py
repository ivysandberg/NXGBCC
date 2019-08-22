'''Goal: use average distance between hashed vectors as an outcome variable.
If avg dist is tiny, our hashing isn't discriminating enough.
Find vicinity radius, r, that maximizes this.'''


# for file/directory operations
from os import listdir
from os.path import isfile, join
import json

# for regular expressions
import re

# for matrix computations
import numpy as np
from Vectorizing_Fingerprints import *

import time



#start = time.time()


# a place to put outputs
my_output_file = "/Users/ivysandberg/dmc/NXGBCC/fingerprint_hashing/dev/vector_output_02_nonrandomheap.txt"



### Prepare the data: only need to do this once...

### read in all fingerprints
my_data_directory = "/Users/ivysandberg/dmc/NXGBCC/fingerprint_hashing/data/figs_0"
# get all the files (not folders) in 'my_data_directory'
onlyfiles = [f for f in listdir(my_data_directory) if isfile(join(my_data_directory, f))]

# select down to the files we want with regular expressions
onlyfilesoftype = [f for f in onlyfiles if re.search("templatedata", f) != None]
# print (onlyfilesoftype)


# create a list of all minutia but separate them by print

# Input: a directory of fingerprint images
# Output: lists of minutia coordinates separated by fingerprint
def minutia_by_indiv_print(onlyfilesoftype):
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

indiv_prints = minutia_by_indiv_print(onlyfilesoftype)


'''
# experiment to get run times

list_of_all_vicinities = []
for i in indiv_prints:
        vicinities = get_vicinities_from_minutia_list(i, 50)
        for vic in vicinities:
            list_of_all_vicinities.append(vic)

# norm_vicinities = normalize_vicinities(list_of_all_vicinities)

start = time.time()
list_of_all_vicinities = []
norm_vicinities = normalize_vicinities(list_of_all_vicinities)
end = time.time()
total_time_for_this_portion = end-start
print ("time: ", total_time_for_this_portion)


vic_lens = []
for i in list_of_all_vicinities:
    vic_lens.append(len(i))
avg_len = np.average(vic_lens)
print ("avg vic len: ", avg_len)


'''


# Input: minutia coords from ALL the fingerprints
# Output: k most distinct vicinities = the representative set
### note: takes about 5min to run
def create_representative_set(all_minutia_from_prints, r, k):
    ### for specified radius, r, identify and normalize (i.e. rotate) vicinities
    list_of_all_vicinities = []
    for i in indiv_prints:
        vicinities = get_vicinities_from_minutia_list(i, r)
        for vic in vicinities:
            list_of_all_vicinities.append(vic)    # ^ when scaled up to 48,000 minutia this took 10+ min
    norm_vicinities = normalize_vicinities(list_of_all_vicinities)
    # runs into the "RuntimeWarning: invalid value encountered in double_scalars" when data is scaled up
    mat = create_vicinity_comparison_scores_matrix(norm_vicinities, norm_vicinities)
    # ^^ takes 2min to run with 443 vicinities
    representative_set_indicies = optimized_k_distinct_vicinities(mat, k)
    # with optimized_k_distinct_vicinities this line runs in < 10 sec
    rep_set = []
    for i in representative_set_indicies:
        rep_set.append(norm_vicinities[i])
    return rep_set


### set the representative set: assign r and k


### hash all fingerprints against the k vicinities


#Input: a single print's minutia data
#Output: an identifying single hashed vector 
def create_print_hashed_vector(print, r):
    # print (len(minutia_of_one_print))
    vics = get_vicinities_from_minutia_list(print, r)
    norm_vics = normalize_vicinities(vics)
    # compare print_1 to the rep_set
    mat_compare_print_to_rep_set = create_vicinity_comparison_scores_matrix(norm_vics, rep_set)
    # print (mat_print_to_rep_set)
    # sum all the columns
    print_vector = mat_compare_print_to_rep_set.min(axis=0)

    return print_vector



### calculate the average distance between the vector representations of fingerprints
# want fingerprint vectors to be different/defining = large distances
# let's use L2 distance, aka Euclidean distance.
# if x and y are numpy arrays, you can do: dist = numpy.linalg.norm(x-y)


# Input: minutia data from many prints (use the function minutia_by_indiv_prints to point to)
# computes a hashed vector for each print (= each list of coordinates) 
# Output: the average distance between each of the hashed vectors. 
def find_avg_distance(prints, r):
    hashed_vectors = []
    for i in prints:
        vector = create_print_hashed_vector(i, r)
        print (vector)
        hashed_vectors.append(vector)
    distances = []
    for i in hashed_vectors:
        for j in hashed_vectors:
            # note: computing distance(i,j) and distance(j,i)
            dist = np.linalg.norm(i-j)
            distances.append(dist)
    avg_dist = np.average(distances)

    return avg_dist


end = time.time()
total_time_for_this_portion = end-start
print ("time: ", total_time_for_this_portion)


### write to file r and this average distance
my_file = open(my_output_file, 'a')  # the 'a' means "append" instead of "overwrite"
note = "I created k = 5, r = 5, and the avg_distance was"
full_note = note + "; " + str(find_avg_distance(indiv_prints, 10)) + "  time: " + str(total_time_for_this_portion)
my_file.write(full_note + "\n")
my_file.close()


print("Great success!!!") # always give yourself a high-five if your code runs!
# whole thing takes just over 2 min to run with 5 print files = 443 vicinities, k = 10
