import unittest
import math
import numpy as np

from triangleClass import *
# from utils import *
# from mathUtils import *

""" These are our 0th order functions, base functions unused by others"""

# The following functions are tests for packing functions

class packingTest(unittest.TestCase):

    # TODO: complete this test
    def testAll(self):

        # TODO: complete the creation of the test packings
        # First we create the packings that we want to operate on
        angle = toRadians(75)
        point1 = [0, 0]
        point2 = [0, 1]
        point3 = [round3(1/(2*math.tan(angle / 2))), 1/2]

        triangle_coordinates = np.asarray([point1, point2, point3])
        first_triangle = triangle(triangle_coordinates, 1)
        first_packing = randomPacking(first_triangle.coordinates, 1, [first_triangle])

        # randomEdgeTest
        assert(first_packing.boundaryDist[first_packing.generateRandomEdge('proposals')] in first_packing.boundaryDist)

        # TODO: complete this test
        # newRadiusTest
        first_packing.packingCenter = np.sum(first_triangle.coordinates, axis=0)
        first_packing.radiusOfGyration =

        # TODO: complete this test
        # updatePackingTest


""" These are our first order functions"""

""" Run the tests"""

if __name__ == '__main__':
    unittest.main()