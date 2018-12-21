import numpy as np
import math


def toRadians(angle):
    """
    Returns radians instead of degrees for angle
    :param angle: angle in degrees
    :return radians: angle in radians
    """
    return angle*math.pi/180

def intersection(boundary, triangle):
    """
    Tells whether or not the new triangle intersects the boundary

    :param boundary: the boundary of the shape created so far
    :param triangle: 3x2 array with the 3 points of the triangle
    :return isIntersection: tells whether or not the new triangle intersects the boundary
    """
    isIntersection = False
    [a, b] = boundary.shape
    for i in range(0, b):
        point1 = boundary[0:2, b]
        point2 = boundary[0:2, b+1]
        mBoundary = (point2[1]-point1[1])/(point2[0] - point1[0])
        bBoundary = point1[1] - mBoundary*point1[0]
        # Check if line 1 intersects


        # Check if line 2 intersects

        # Check if line 3 intersects
    return isIntersection