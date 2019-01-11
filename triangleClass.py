import numpy as np
import growthMethods

from mathUtils import *

class randomPacking:
    def __init__(self, boundary, triangleCount, triangleList):
        """
        Creates an object that defines the entire random packing
        :param boundary: the boundary for the region, a numpy array defining a set of points
        :param triangleCount: the number of triangles in the packing, an int
        :param triangleList: a list of all of the individual triangles, a list
        """
        self.boundary = boundary
        self.triangleCount = triangleCount
        self.triangleList = triangleList

    def generateRandomEdge(self, method):
        """
        Generates the random edge according to the method
        :param method: A string describing the growth method for the randomPacking
        :return: a random edge on the boundary of the randomPacking
        """
        if method == 'uniform':
            edge = growthMethods.uniformDist(self.boundary)
            return edge

    def insertTriangle(self, addedTriangle, edge):
        """
        Inserts a valid triangle onto the object and updates the boundary accordingly
        :param addedTriangle: triangle to insert
        :param edge: edge to insert onto
        :return: the updated boundary and triangle count
        """
        # Need to handle the case where we're not actually adding on any points
        print('adding to edge: ', edge)
        edgePoint1, edgePoint2 = self.boundary[edge], self.boundary[(edge + 1) % len(self.boundary)]
        addedPoint = None
        for i in range(0, 3):
            if not np.array_equal(addedTriangle.coordinates[i], edgePoint1) and \
                    not np.array_equal(addedTriangle.coordinates[i], edgePoint2):
                addedPoint = addedTriangle.coordinates[i]

        newBoundary = self.boundary[:edge] + [addedPoint] + self.boundary[(edge + 1):]

        self.triangleCount += 1
        self.boundary = newBoundary


class triangle:
    def __init__(self, coordinates, serialNumber):
        """
        The triangle class
        :param coordinates: coordinates of the triangle, a numpy list
        :param serialNumber: the index of the triangle added onto the packing, an int
        """
        self.coordinates = coordinates
        self.serialNumber = serialNumber

    def setCoordinates(self, coordinates):
        self.coordinates = coordinates

    def setSerialNumber(self, serialNumber):
        self.serialNumber = serialNumber
