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
    return round3(math.sqrt(vector[0]*vector[0] + vector[1]*vector[1]))

def pointSlopeForm(point1, point2):
    """
    Returns the point slope form of a line between 2 points
    :param point1: The first point
    :param point2: The second point
    :return: The slope, m, and the y intercept, b
    """
    if point2[0] == point1[0]:
        m = (point2[1] - point1[1])/0.001
        b = point1[1] - m*point1[0]
    else:
        m = (point2[1] - point1[1]) / (point2[0] - point1[0])
        b = point1[1] - m * point1[0]
    return round3(m), round3(b)

def intersection(proposals, packing):
    """
    Tells whether or not the new triangle intersects any of the other triangles in the packing

    :param proposals: the proposal triangles that could be added to the packing, a list of triangle objects
    :param packing: the existing packing with all of its triangles, a packing object
    :return intersections: tells which of the new proposal triangles intersects the existing
    packing triangles, a list of booleans
    """
    intersections = [False]*len(proposals)

    # Go through each triangle in the packing and check to make sure there isn't an intersection
    for i in range(0, packing.triangleCount):
        for j in range(0, len(intersections)):
            if not intersections[j]:
                intersections[j] = triangleIntersection(packing.triangleList[i], proposals[j])

    return intersections

def fastIntersection(proposal, packing, randomEdge):
    """
    A method that checks the intersections of the proposal triangles in O(1) time. Used as an alternative to
    the above method

    :param proposals: the proposal triangles that could be added to the packing, a list of triangle objects
    :param packing: the existing packing with all of its triangles, a packing object
    :param randomEdge: the edge that we're adding onto
    :return intersections: tells which of the new proposal triangles intersects the existing
    packing triangles, a list of booleans
    """
    # If the length of the packing boundary is less than 7, then just brute force it for not
    if len(packing.boundary) < 7:
        for i in range(0, len(packing.boundary)):
            m_edge, b_edge = pointSlopeForm()
            for j in range(3):
                m1, b1 = pointSlopeForm(proposal.coordinates[i], proposal.coordinates[(i+1) % 3])

    else:


def triangleIntersection(triangle1, triangle2):
    """

    :param triangle1: a triangle object
    :param triangle2: a triangle object
    :return intersection: whether or not the triangles intersect
    """
    individualIntersection = False
    sameEdges = 0
    for i in range(0,3):
        for j in range(0, 3):
            # Compute the point slope form of the first line
            point11 = triangle1.coordinates[i]
            point12 = triangle1.coordinates[(i+1) % 3]
            m1, b1 = pointSlopeForm(point11, point12)

            # Compute the point slope form of the second line
            point21 = triangle2.coordinates[j]
            point22 = triangle2.coordinates[(j+1) % 3]
            m2, b2 = pointSlopeForm(point21, point22)

            # Check if they're the same lines
            if m1 == m2 and b1 == b2:
                sameEdges += 1
            elif m1 != m2:
                x = - (b2 - b1)/(m2 - m1)
                if (x < point11[0]) and (x > point12[0]):
                    individualIntersection = True
                    break
    if sameEdges == 3:
        individualIntersection = True

    return individualIntersection

def triangleLineIntersection(triangle, point1, point2):

    intersection = False


    return intersection

def round3(number):
    """
    Returns a number rounded to 3 digits
    :param number:
    :return: rounded number
    """
    if type(number) is np.float64:
        return np.round_(number, 3)
    elif type(number) is float:
        return round(number, 3)
    elif type(number) is list and type(number[0]) is np.float64:
        return [np.round_(elem, 3) for elem in number]
    elif type(number) is list:
        return [round(elem, 3) for elem in number]

def floatMult(number, array):
    """
    Multiplies an array elementwise by a float
    """
    return [number*i for i in array]