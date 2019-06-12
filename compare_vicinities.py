import numpy as np
import pandas as pd
import itertools
import math

# from yapf.yapflib.yapf_api import FormatFile
# FormatFile(
#     "/Users/ivysandberg/dmc/NXGBCC/dev/compare_vicinities.py", in_place=True)
# ''' Compare Vicinities '''

V1 = [(3, 2, 0.5), (-4, 1.2, 3.4), (5, -2, 4.4), (1, -1, 1), (3.6, 2, 0.5),
      (1, 1, 1), (-2, -0.3, 0), (1, 2, -4)]

V2 = [(2.2, 3.1, -1.2), (-2, 3, 1), (1, 1, 1), (1, -1, -1), (4, 0, 0)]

#for i in list_of_minutiae:
#print (i[2])

for p, q in itertools.product(V1, V2):
    print(p, q)

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

    print(pairing_score(p, q))

    # Create matrix P x Q of all the pairing scores

    
