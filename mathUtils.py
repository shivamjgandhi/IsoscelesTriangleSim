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

    return isIntersection