import numpy as np

class randomPacking:
    def __init__(self, boundary, triangleCount, triangleList):
        """
        Creates an object that defines the entire random packing
        :param boundary: the boundary for the region
        :param triangleCount: the number of triangles in the packing
        :param triangleList: a list of all of the individual triangles
        """
        self.boundary = boundary
        self.triangleCount = triangleCount
        self.triangleList = triangleList

class triangle:
    def __init__(self, coordinates):
        self.coordinates = coordinates

    def setCoordinates(self, coordinates):
        self.coordinates = coordinates

class boundaryObject:
    def __init__(self, boundary, triangleCount):
        """
        Creates an object that defines the boundaries of the packing
        :param boundary: the boundary is a numpy array defining the external boundary of the object
        :param triangleCount: the triangleCount tells how many triangles it was made from
        """
        self.boundary = boundary
        self.triangleCount = triangleCount

    def insertTriangle(self, triangle, edge):
        """
        Inserts a valid triangle onto the object and updates the boundary accordingly
        :param triangle: triangle to insert
        :param edge: edge to insert onto
        :return: the updated boundary and triangle count
        """
        # Return the actual edge
        edgePoint1, edgePoint2 = self.boundary[edge-1], self.boundary[edge % (self.triangleCount + 2)]
        self.triangleCount += 1

    def removeEdge(self, edge):
        self.boundary = np.array(self.boundary[:edge], self.boundary[edge+1:])