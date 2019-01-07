import numpy as np
import math
from triangleClass import *

def toRadians(angle):
    """
    Returns radians instead of degrees for angle
    :param angle: angle in degrees
    :return radians: angle in radians
    """
    return angle*math.pi/180

def norm(vector):
    """
    Returns the norm of a vector
    :param vector: input vector
    :return: norm
    """
    return math.sqrt(vector[0]*vector[0] + vector[1]*vector[1])

def pointSlopeForm(point1, point2):
    """
    Returns the point slope form of a line between 2 points
    :param point1: The first point
    :param point2: The second point
    :return: The slope, m, and the y intercept, b
    """
    m = (point2[1] - point1[1])/(point2[0] - point1[0])
    b = point1[1] - m*point1[0]
    return round3(m), round3(b)

def intersection(addedBoundary, addedTriangle):
    """
    Tells whether or not the new triangle intersects the boundary

    :param addedBoundary: the boundary object of the shape created so far
    :param addedTriangle: a triangle object
    :return isIntersection: tells whether or not the new triangle intersects the boundary
    """
    isIntersection = False
    boundary = addedBoundary.boundary
    b = np.asarray(boundary).shape[0]
    sameEdges = []
    for i in range(0, b):
        point1 = boundary[i]
        point2 = boundary[(i + 1) % b]
        # Compute the point slope form of the boundary that we're adding onto
        mBoundary, bBoundary = pointSlopeForm(point1, point2)
        for j in range(0, 3):
            # Check if line 1 intersects
            m, b = pointSlopeForm(addedTriangle.coordinates[j], addedTriangle.coordinates[(j+1) % 3])
            if m == mBoundary:
                print('Matching slopes')
                if bBoundary == b:
                    sameEdges = sameEdges + [1]
            else:
                x = - (bBoundary - b)/(mBoundary - m)
                if (x < addedTriangle.coordinates[j][0]) or (x > addedTriangle.coordinates[(j+1) % 3][0]):
                    isIntersection = True
        if sum(sameEdges) == 3:
            isIntersection = True
        if isIntersection:
            break
    return isIntersection

def round3(number):
    """
    Returns a number rounded to 3 digits
    :param number:
    :return: rounded number
    """
    if type(number) is float:
        return round(number, 3)
    elif type(number) is list:
        return [round(elem, 3) for elem in number]