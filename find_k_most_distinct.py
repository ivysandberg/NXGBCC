
'''
Pick k points that aren't too close together algo:

1) pick a point at random
2) pick the point furthest from the first point
for (3-k): pick the next point with the largest sum of distances to all points in the picked set (aka heap)

Requirements: This function should take in a symmetric or half-filled numpy array (with other half across the diagonal being all zeros), and have k be a parameter the user inputs as well.
The output should be the indexes of the points selected, e.g. [2, 13, 18, 35…]


'''

import numpy as np
import matplotlib.pyplot as plt
import math
import itertools

mat = np.zeros((4,4))
mat[1,2] = 1
mat[1,3] = 2
mat[2,1] = 3
mat[2,3] = 4
mat[3,1] = 4
mat[3,2] = 5
mat[3,3] = 7
mat[0,0] = 2
mat[0,1] = 3
mat[0,3] = 1

# print (mat)


# pick a row at random:
r = np.random.choice(range(mat.shape[0])) # randomly selects a row by index
# r = mat[np.random.choice(mat.shape[0], 2, replace=False), :] # randomly selects a whole row


# let's keep our selections in a np array called selections
heap = []
heap.append(r)


# find the point furthers from point r
second_point = np.argmax(mat[r])
# print(second_point)
heap.append(second_point)


# for all future points:
# make a vector that is sum of points already selected
# find the max
# make sure the argmax of that vector isn't already in heap
# if not, cool, if next...


# heap = []
# r = np.random.choice(range(mat.shape[0]))
# heap.append(r)
# next_vic = np.argmax(mat[r])
# heap.append(next_vic)
#
# # print ("Heap: ", heap)
#
# min_distances_to_the_heap = []
#
# # iterate k times
# # for _ in range(k):
#
# best_point = -1
#
# for i in range(mat.shape[1]):
#     if i in heap:
#         continue
#     this_point_min_dist_to_heap = np.min(mat[heap,i])
#     min_distances_to_the_heap.append(this_point_min_dist_to_heap)
#     current_largest_min = np.max(min_distances_to_the_heap)
#     if this_point_min_dist_to_heap == current_largest_min:
#         best_point = i
#
# heap.append(best_point)

# print ("Min Distances to Heap: ", min_distances_to_the_heap)
# print ("Current largest min: ", current_largest_min)
# print ("Point to append to heap: ", best_point)
#
# print ("New Heap: ", heap)


# ^^^ this works

# now dow it with the graphs

coords = [(1,1), (2,3), (4,1), (5,5), (3,9), (4,4)]
x = 1, 2, 4, 5, 3, 4
y = 1, 3, 1, 5, 9, 4
plt.scatter(x,y)
# plt.show()

def dist(i,j):
    dist = math.sqrt((j[0] - i[0])**2 + (j[1] - i[1])**2)
    output = np.array(dist)
    return output

def matrix(coords):
    output_matrix = np.zeros((len(coords), len(coords)))
    for i in range(len(coords)):
        for j in range(len(coords)):
            matrix_value = dist(coords[i], coords[j])
            output_matrix[i, j] = matrix_value
    return output_matrix

matrix = matrix(coords)
print (matrix)



def k_distinct_coords(matrix, k):
    heap = []
    r = np.random.choice(range(matrix.shape[0]))
    heap.append(r)
    next_vic = np.argmax(matrix[r])
    heap.append(next_vic)

    print ("Starting heap: ", heap)

    # iterate k times
    for _ in range(k):
        min_not_in_heap = []
        best_point = 'done'

        for i in range(matrix.shape[1]):
            if i in heap:
                continue
            this_point_min_dist_to_heap = np.min(matrix[heap,i])    # min distance to a point in the heap
            # print ("Heap, i: ", np.min(matrix[heap,i]))
            min_not_in_heap.append(this_point_min_dist_to_heap)
            current_largest_min = np.argmax(min_not_in_heap)
            if this_point_min_dist_to_heap >= current_largest_min:
                best_point = i
        heap.append(best_point)
    return heap
