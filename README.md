# NXGBCC

construct-vicinities: defines the change of coordinates for all the minuitae points for a single minuitae point set at the origin

compress-vicinities: returns the change of coordinates (for mj) for every pair (mi, mj) in the list of minuitae points. 

compare-vicinities: takes two given vicinities (which will contain a list of the new coordinates from prior code ^^ but this has yet to be implemented). Defines matrix (Vicinity1, Vicinity2) and each cell contains a pairing score between each pair of minuitae points in each of the vicinities. Perfect pair score (identical coords) = 1.  
Implements the Hungarian algorithm. Creates a matrix filled with values = 1 - pairing score (in order to find maximums bc, Hungarian is designed to minimize cost (in the assignment problem)) --> thus identical vicinities after implementation of the Hungarian alg would result in sum of pairing scores = 0.  
