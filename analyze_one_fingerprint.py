'''
Analyzing real minutia data from fingerprint images

'''

import numpy as np
import math
from Vectorizing_Fingerprints import *



# test with json file of minutiae data
import json

with open('/Users/ivysandberg/dmc/NXGBCC/data/fingerprints/f0001_01.png.templatedata') as json_file:
    data = json.load(json_file)

json_string = json.dumps(data, indent=4)
# print (json_string)


# print (data)

# for key, value in data.items():
# 	print(key + ':' , value)


# minutiae = data['minutiae']


# Input: the data from one json file = data from one fingerprint image
# Output: the minutiae data as a list of coordinates
def create_minutiae_coords(data):
    list_of_minutiae_coords = []
    minutiae = data['minutiae']
    for i in minutiae:
        # print (i)
        output = (i['x'], i['y'], i['direction'])
        # print (output)
        list_of_minutiae_coords.append(output)
    return list_of_minutiae_coords



minutia_list = create_minutiae_coords(data)
# print (minutia_list)
# print (len(minutia_list))

x_coords = []
y_coords = []
theta_coords = []
for i in minutia_list:
    x_coords.append(i[0])
    y_coords.append(i[1])
    theta_coords.append(i[2])

print ("X Coord Mean: ", np.mean(x_coords), "St dev: ", np.std(x_coords), "Min: ", np.min(x_coords), "Max: ", np.max(x_coords))
print ("Y Coord Mean: ", np.mean(y_coords), "St dev: ", np.std(y_coords), "Min: ", np.min(y_coords), "Max: ", np.max(y_coords))
print ("Theta Mean: ", np.mean(theta_coords), "St dev: ", np.std(theta_coords), "Min: ", np.min(theta_coords), "Max: ", np.max(theta_coords))


vicinities = get_vicinities_from_minutia_list(minutia_list, 50)
# print ("Vicinities: ", vicinities)


nom_vicinities = normalize_vicinities(vicinities)
# print (nom_vicinities)

list = filter_list_of_vinities_by_size(nom_vicinities, 4, 7)
# print (list)

comp_scores_matrix = create_vicinity_comparison_scores_matrix(list, list)
print (comp_scores_matrix)
