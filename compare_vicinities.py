import numpy as np
import pandas as pd
import itertools
import math
from munkres import Munkres, print_matrix

# from yapf.yapflib.yapf_api import FormatFile
# FormatFile(
#     "/Users/ivysandberg/dmc/NXGBCC/dev/compare_vicinities.py", in_place=True)
''' Compare Vicinities '''

V1 = [(3, 2, 0.5), (-4, 1.2, 3.4), (5, -2, 4.4), (1, -1, 1), (3.6, 2, 0.5),
      (1, 1, 1), (-2, -0.3, 0), (1, 2, -4), (0, 0, 0)]

V2 = [(2.2, 3.1, -1.2), (-2, 3, 1), (1, 1, 1), (1, -1, -1), (4, 0, 0), (0, 0,
                                                                        0)]

#for i in list_of_minutiae:
#print (i[2])

for p, q in itertools.product(V1, V2):
    #print(p, q)

    #for p, q in zip(V1, V2):   # doesn't match ALL possible combinations

    def pairing_score(p, q):
        #output = []
        sd_x = 2.3
        sd_y = 4.1
        sd_theta = 1.2
        score = math.exp((-((p[0] - q[0])**2) / sd_x)) * math.exp(
            (-((p[1] - q[1])**2) / sd_y)) * math.exp(
                (-((p[2] - q[2])**2) / sd_theta))
        output = np.array(score)

        #output.append(score)
        return output

    #print(pairing_score(p, q))

    # determine max pairing score (it's 1)
    # if pairing_score(p, q) >= 0.9:
    #     print(pairing_score(p, q))
    # if pairing_score(p,q) >= 0.1:
    #     print(p,q)

    # Create matrix P x Q of all the pairing scores


def vicinities_matrix(Vi, Vj):
    output_matrix = np.zeros((len(Vi), len(Vj)))
    for i in range(len(Vi)):
        for j in range(len(Vj)):
            matrix_value = pairing_score(Vi[i], Vj[j])
            output_matrix[i, j] = matrix_value
    return output_matrix


print("This is the comparison matrix Vicinity 1 to Vicinity 2: ")
print(vicinities_matrix(V1, V2))

# create a new vicinity to comparison

V3 = [(0,0,0),(4, 2, 0.3), (2, 3.3, -2), (-2, 1, 1), (-4, 0, 1)]
V4 = [(1.1, 2.2, 3.3), (-1.1, -2.2, -3.3), (1.1, 2.2, 3.3), (1,1,1)]

print("This is the comparison matrix of Viciinty 3 to Vicinity 4: ")
print(vicinities_matrix(V3, V4))

#for pi, qj in itertools.product(V3, V4):
#print (pairing_score(pi,qj))

# Implement Hungarian Method to determine the best pairings between minutiae


def Hungarian(matrix):
    # Create matrix filled with values = 1 - pairing score, in order to find maximums
    # Hungarian is designed to minimize cost (in the assignment problem)
    inv_matrix = []
    for row in matrix:
        #print(row)
        inv_row = []
        for col in row:
            inv_row += [1 - col]
        inv_matrix += [inv_row]
        #print(inv_matrix)

    m = Munkres()  # Munkres = Hungarian
    indexes = m.compute(inv_matrix)
    #print_matrix(inv_matrix, msg='compute highest pairing scores')
    total = 0
    for row, col in indexes:
        val = inv_matrix[row][col]
        total += val
        print(f'({row}{col}) -> {val}')
    return (f'sum of scores: {total}')
    #return inv_matrix, total


# compute total pairing score for vicinities 3 and 4
M2 = vicinities_matrix(V3, V4)
print("This is the Hungarian alg on V3 and V4")
print(Hungarian(M2))

# compute total pairing score for vicinities 1 and 2
M1 = (vicinities_matrix(V1, V2))
print("This is the Hungarian alg on V1 and V2")
print(Hungarian(M1))
'''
Total Pairing Score is used to compare vicinities of the same fingerprint or of different fingerprint
'''


# run test to compare to vicinities that are the same (result: sum of pair scores = 0)
# this is bc we do 1-pairing score to find maxs (best pairs)
V4 = [(1.1, 2.2, 3.3), (-1.1, -2.2, -3.3), (1.1, 2.2, 3.3)]
V5 = [(1.1, 2.2, 3.3), (-1.1, -2.2, -3.3), (1.1, 2.2, 3.3)]

M3 = (vicinities_matrix(V4, V5))
print("This is the Hungarian alg on V4 and V5 (identical)")
print (Hungarian(M3))
