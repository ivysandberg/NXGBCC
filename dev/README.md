# NXGBCC

For an indepth explanation to the Fingerprint hash approach here: https://docs.google.com/document/d/1OkhVMNrllXu0T6eg3PnfpK556qRn8zn8bwPBMbIStD0/edit?usp=sharing


Vectorizing_Fingerprints: contains functions for vectorizing fingerprints into rotationally and translationally invariant representative sets of vicinities. 

Vectorizing_Fingerprints02: optimizes previous functions 


Files to test the Fingerprint hash algo: 
hashed_vectors
vic_compare_exp



The following files were scratch files to build the hash algorithm:


construct_vicinities: defines the change of coordinates for all the minuitae points for a single minuitae point set at the origin


compare_vicinities: Defines matrix (Vicinity1 x Vicinity2) and each cell contains a pairing score between each pair of minuitae points in each of the vicinities. Perfect pair score (identical coords) = 1.  
Implements the Hungarian algorithm. Creates a matrix filled with values = 1 - pairing score (in order to find maximums bc, Hungarian is designed to minimize cost (in the assignment problem)) --> thus identical vicinities after implementation of the Hungarian alg would result in sum of pairing scores = 0.  

find_k_most_distinct
