import unittest
import math
import numpy as np

from utils import *
from mathUtils import *

class UtilsTest(unittest.TestCase):

    def test_generateIndividualTriangle(self):
        # boundary, edge, angle
        boundary = np.array([[0, 0],
                             [1, 0]])
        edge = 1
        angle = 60
        triangleCoordinates = generateIndividualTriangle(boundary, edge, angle)
        self.assertEqual(triangleCoordinates, [[0,0],
                                               [1,0],
                                               [0.5, 0.866]])
        angle = 90
        triangleCoordinates = generateIndividualTriangle(boundary, edge, angle)
        self.assertEqual(triangleCoordinates, [[0,0],
                                               [1,0],
                                               [0.5, 0.707]])

class mathUtilsTest(unittest.TestCase):

    def testToRadians(self):
        self.assertEqual(toRadians(90), math.pi/2)
        self.assertEqual(toRadians(0), 0)
        self.assertEqual(toRadians(180), math.pi)

    def testNorm(self):
        self.assertEqual(norm([1, 1]), math.sqrt(2))
        self.assertEqual(norm([1, 0]), 1)
        self.assertEqual(norm([0, 1]), 1)

    #def testIntersection(self):
        # boundary, triangle, edge

class triangleClassTest(unittest.TestCase):

    def testTriangle(self):
        newTriangle = triangle([[0,0],
                               [1,0],
                               [1,1]])
        self.assertEqual(newTriangle.coordinates, [[0,0],
                                                   [1,0],
                                                   [1,1]])

        newBoundary = boundaryObject(newTriangle.coordinates, 1)
        self.assertEqual(newBoundary.boundary, newTriangle.coordinates)
        self.assertEqual(newBoundary.triangleCount, 1)

if __name__ == '__main__':
    unittest.main()