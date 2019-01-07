import unittest
import math
import numpy as np

from utils import *
from mathUtils import *


class UtilsTest(unittest.TestCase):

    @staticmethod
    def test_generateIndividualTriangle():
        # boundary, edge, angle
        boundary = np.array([[0.0, 0.0],
                             [1.0, 0.0]])
        edge = 1
        angle = 60
        newBoundary = generateIndividualTriangle(boundary, edge, angle)
        areEqual = np.array_equal(newBoundary, [[0.0, 0.0],
                                                [1.0, 0.0],
                                                [0.5, 0.866]])
        print("The first boundary for the 60 degree triangles is: ", areEqual)

        angle = 90
        boundary = generateIndividualTriangle(boundary, edge, angle)
        areEqual = np.array_equal(boundary, [[0.0, 0.0],
                                             [1.0, 0.0],
                                             [0.5, 0.5]])
        print("The first boundary for the 90 degree triangles is: ", areEqual)


class mathUtilsTest(unittest.TestCase):

    def testToRadians(self):
        self.assertEqual(toRadians(90), math.pi / 2)
        self.assertEqual(toRadians(0), 0)
        self.assertEqual(toRadians(180), math.pi)
        print('all radians test passed')

    def testNorm(self):
        self.assertEqual(norm([1, 1]), math.sqrt(2))
        self.assertEqual(norm([1, 0]), 1)
        self.assertEqual(norm([0, 1]), 1)
        print('all norm test passed')

    def testIntersection(self):
        newBoundary = boundaryObject()
        newBoundary.setBoundary([[0.0, 0.0],
                                 [1.0, 0.0],
                                 [0.5, 0.5]])


class triangleClassTest(unittest.TestCase):

    def testTriangle(self):
        newTriangle = triangle([[0, 0],
                                [1, 0],
                                [1, 1]])
        self.assertEqual(newTriangle.coordinates, [[0, 0],
                                                   [1, 0],
                                                   [1, 1]])

        newBoundary = boundaryObject(newTriangle.coordinates, 1)
        self.assertEqual(newBoundary.boundary, newTriangle.coordinates)
        self.assertEqual(newBoundary.triangleCount, 1)


if __name__ == '__main__':
    unittest.main()
