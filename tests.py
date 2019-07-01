import unittest
import math
import numpy as np
from Vectorizing_Fingerprints import new_mj_coords, compute_distance, angle_between


# test change of coordinates
# input = (mi, mj)

# check by hand
print ("dist: ", compute_distance((0,1,1.5708), (0,2,1.5708)))
print ("angle: ", angle_between((0,1,1.5708), (0,2,1.5708)))
print(new_mj_coords((0,1,1.5708), (0,2,1.5708)))


class TestCoords(unittest.TestCase):    # TestCoords = class which inherits unittest.TestCase
    def test_coords_same_direction(self):       # test_change_coords = method
        """
        Test that it can change a set of coordinates correctly
        """
        mi = (0,1,1.5708)
        mj = (0,2,1.5708)
        result = new_mj_coords(mi, mj)
        self.assertEqual(result, (0.3317021649341794, 0.9433841602327115, 0.0))

        '''
        the method .assertEqual(a,b) is equivalent to a == b
        other methods include: .assertIs(a,b) = a is b, .assertIsNone(x) = x is None,
        .assertIn(a,b) = a in b, and .assertIsInstance(a,b) = isinstance(a, b)


        '''
    def test_int(self):
        '''
        Test coordinates with integers
        '''

        mi = (1,2,3)
        mj = (2,3,4)
        result = new_mj_coords(mi,mj)
        self.assertEqual(result, (-1.3654156128202577, 0.3682936386454162, 1))



# if __name__ == '__main__':
#     unittest.main()


# check graphically
import matplotlib.pyplot as plt

x = 0, 0, 0, 0.3317021649341794
y = 1, 2, 0, 0.9433841602327115
fig = plt.figure()
ax = fig.add_subplot(111)
ax.scatter(x, y)
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)

plt.show()


# test compare minutia
