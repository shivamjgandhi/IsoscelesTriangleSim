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
        """
        if boundaryDist is None:
            self.boundaryDist = np.asarray([0, 1, 1, 2, 2])
        """
        self.boundary = boundary
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
            # this number will be in range(0, len(boundarydist) - 1)
            proposition = growthMethods.uniformAcrossProposals(self.boundaryDist)
            return proposition

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

    def updatePacking(self, proposal, growth_edge, match_point=None):
        """
        :param proposal: the proposal triangle being added to the packing, a triangle object
        :param match_point: the index of the boundary point with the match for the new point, an int or bool
        :param growth_edge: the index of the boundary point being added on to, an int
        :return: the updated packing, a packing object
        """
        # update packing center and radius of gyration
        new_center = 1 / (self.triangleCount + 1) * (self.triangleCount * self.packingCenter + proposal.center)
        new_radius_gyration = self.computeNewRadiusGyration(proposal)
        self.packingCenter = new_center
        self.radiusOfGyration = new_radius_gyration

        # update the boundary
        # return the new point
        first_point = self.boundary[growth_edge]
        second_point = self.boundary[(growth_edge + 1) % len(self.boundary)]
        new_point = None
        for point in proposal.coordinates:
            if not comparePoints(first_point, point) and not comparePoints(second_point, point):
                new_point = point

        # check if there is a matchPoint. If there is none, then the update is much easier
        if match_point:
            # in case the match point is ahead of the growth edge in the packing
            if match_point > growth_edge:
                self.boundary = np.concatenate((self.boundary[0:growth_edge], self.boundary[match_point:]))
            # in case the match point is before the growth edge in the packing
            else:
                self.boundary = np.concatenate((self.boundary[0:match_point], self.boundary[growth_edge+1:]))
        else:
            self.boundary = np.concatenate((self.boundary[0:growth_edge], [new_point], self.boundary[growth_edge:]),
                                          axis=0)

        """
        # As for updating the boundary dist, this depends on what length edges were added. It also depends on whether
        # there was a match point or not
        if match_point:
            # if the match comes before
            if match_point < growth_edge:
                k = growth_edge - match_point
                side_length = norm(second_point - match_point)
                # If we're adding to a side of length 1
                if (side_length > 0.98) and (side_length < 1.02):

                else:

            # if the match comes after
            else:

        else:
            side_length = norm(first_point - second_point)
            # check if we're adding to a side of length 1
            if (side_length > 0.98) and (side_length < 1.02):
                # update distribution accordingly
                a = self.boundaryDist[random_edge]
                self.boundaryDist = np.concatenate((self.boundaryDist[0:random_edge] , [a, a+1, a+1] , (1 + self.boundaryDist[random_edge:])))
            else:
                # added to a side of length not equal to 1. Thus check if right or left triangle
                # if right side
                if self.boundaryDist[random_edge] == self.boundaryDist[(random_edge - 1) % len(self.boundaryDist)]:
                    a = self.boundaryDist[random_edge]
                    self.boundaryDist = np.concatenate((self.boundaryDist[0:random_edge-1],
                                                        [a + 1, a + 1],
                                                        1 + self.boundaryDist[random_edge + 1:]))
                # if left side
                else:
                    a = self.boundaryDist[random_edge]
                    self.boundaryDist = np.concatenate((self.boundaryDist[0:random_edge],
                                                        [a + 1],
                                                        1 + self.boundaryDist[random_edge + 1:]))
        """
        return self


class triangle:
    def __init__(self, coordinates, serialNumber):
        """
        The triangle class
        :param coordinates: coordinates of the triangle, a numpy list
        :param serialNumber: the index of the triangle added onto the packing, an int
        """
        self.coordinates = coordinates
        self.serialNumber = serialNumber
        self.center = self.computeCenter()

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
