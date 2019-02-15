import numpy as np
import growthMethods
from mathUtils import *


class randomPacking:
    def __init__(self, boundary, triangleCount, triangleList, boundaryDist=None, radius=0, center=None):
        """
        Creates an object that defines the entire random packing
        :param boundary: the boundary for the region, a numpy array defining a set of points
        :param triangleCount: the number of triangles in the packing, an int
        :param triangleList: a list of all of the individual triangles, a list
        """

        self.boundary = boundary
        self.boundaryDist = boundaryDist
        self.triangleCount = triangleCount
        self.triangleList = triangleList
        self.radiusOfGyration = radius
        self.packingCenter = center

    def generateRandomEdge(self, method):
        """
        Generates the random edge according to the method
        :param method: A string describing the growth method for the randomPacking
        :return: a random edge on the boundary of the randomPacking
        """
        if method == 'uniform':
            edge = growthMethods.uniformDist(self.boundary)
            return edge
        if method == 'proposals':
            return growthMethods.uniformAcrossProposals(self.boundaryDist)

    def insertTriangle(self, addedTriangle, edge):
        """
        Inserts a valid triangle onto the object and updates the boundary accordingly
        :param addedTriangle: triangle to insert
        :param edge: edge to insert onto
        :return: the updated boundary and triangle count
        """
        # Need to handle the case where we're not actually adding on any points
        # print('adding to edge: ', edge)
        edgePoint1, edgePoint2 = self.boundary[edge], self.boundary[(edge + 1) % len(self.boundary)]
        addedPoint = None
        for i in range(0, 3):
            if not np.array_equal(addedTriangle.coordinates[i], edgePoint1):
                if not np.array_equal(addedTriangle.coordinates[i], edgePoint2):
                    addedPoint = addedTriangle.coordinates[i]

        newBoundary = np.concatenate((self.boundary[:(edge + 1)], [addedPoint], self.boundary[(edge + 1):]), axis=0)

        self.triangleCount += 1
        self.boundary = newBoundary
        self.triangleList.append(addedTriangle)

    def computeNewRadiusGyration(self, new_proposal):
        """
        Computes the new radius of gyration given a new proposal triangle via an update equation
        :param new_proposal: new triangle that will determine the new radius of gyration
        :return new_radius_gyration
        """
        new_center = 1 / (self.triangleCount + 1) * (self.triangleCount * self.packingCenter + new_proposal.center)
        new_radius_gyration = -np.dot(new_center, new_center) + self.triangleCount / (self.triangleCount + 1) * \
                              (self.radiusOfGyration + np.dot(self.packingCenter, self.packingCenter)) + \
                              np.dot(new_proposal.center, new_proposal.center) / (self.triangleCount + 1)
        return new_radius_gyration


class triangle:
    def __init__(self, coordinates, serialNumber):
        """
        The triangle class
        :param coordinates: coordinates of the triangle, a numpy list
        :param serialNumber: the index of the triangle added onto the packing, an int
        """
        self.coordinates = coordinates
        self.serialNumber = serialNumber
        self.center = self.computeCenter(coordinates)

    def setCoordinates(self, coordinates):
        self.coordinates = coordinates

    def setSerialNumber(self, serialNumber):
        self.serialNumber = serialNumber

    def computeCenter(self):
        x_accum = 0
        y_accum = 0
        for i in range(3):
            y_accum += self.coordinates[i, 0] / 3
            x_accum += self.coordinates[i, 1] / 3
        return np.asarray([y_accum, x_accum])
