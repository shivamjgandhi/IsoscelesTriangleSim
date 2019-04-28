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
    return round3(np.linalg.norm(vector))

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

    :param proposal: the proposal triangles that could be added to the packing, a list of triangle objects
    :param packing: the existing packing with all of its triangles, a packing object
    :param randomEdge: the edge that we're adding onto
    """
    # If the length of the packing boundary is less than 7, then just brute force it for not
    if len(packing.boundary) < 7:
        for i in range(0, len(packing.boundary)):
            m_edge, b_edge = pointSlopeForm(packing.boundary[i], packing.boundary[(i+1) % len(packing.boundary)])
            for j in range(3):
                m1, b1 = pointSlopeForm(proposal.coordinates[i], proposal.coordinates[(i+1) % 3])
            # Check if they're the same lines
            if m1 == m_edge and b1 == b_edge:
                sameLine = i
            elif m1 != m_edge:
                x = - (b_edge - b1) / (m_edge - m1)
                if (x < packing.boundary[i, 0]) and (x > packing.boundary[(i+1) % len(packing.boundary), 0]):
                    isIntersection = True
                    break
    else:
        oldPoint1 = packing.boundary[randomEdge]
        oldPoint2 = packing.boundary[(randomEdge + 1)]
        # Forward direction
        # for i in range(4):

    # Backward direction
    # for i in range(3):

    return isIntersection, sameLine

def boundaryIntersection(packing, proposal_triangle, growthEdge):
    """
    This algorithm checks if any of the lines on the boundary of the packing intersect with the proposal triangle. It
    also returns points that match the new point in the triangle.

    :param packing: the packing with the boundary, a packing object
    :param proposal_triangle: the proposal triangle, a triangle object
    :param growthEdge: the edge the triangle is growing off of, an int
    :return:
    """

    # Figure out which point is the new point
    first_point = packing.boundary[growthEdge]
    second_point = packing.boundary[(growthEdge + 1) % len(packing.boundary)]
    new_point = None
    for point in proposal_triangle.coordinates:
        if not comparePoints(first_point, point) and not comparePoints(second_point, point):
            new_point = point

    # Put the new lines being added from new_point in slope intercept form
    m1, b1 = pointSlopeForm(first_point, new_point)
    m2, b2 = pointSlopeForm(second_point, new_point)

    # Go through each edge in the boundary and check if any of the new lines on the triangle intersect. Also check if
    # new_point matches any of these points
    matchPoint = None
    isIntersection = False
    for edge_pt in range(len(packing.boundary)):
        m_edge, b_edge = pointSlopeForm(packing.boundary[edge_pt],
                                        packing.boundary[(edge_pt + 1) % len(packing.boundary)])
        # check if new_point is the same as edge_pt
        if comparePoints(packing.boundary[edge_pt], new_point):
            matchPoint = edge_pt

        # check if the edge intersects with either side of the triangle
        if m_edge == m1 and b_edge == b1:
            pass
        elif m_edge != m1:
            x = - (b1 - b_edge) / (m1 - m_edge)
            if (x < first_point[0]) and (x > new_point[0]):
                isIntersection = True
                break

        if m_edge == m2 and b_edge == b2:
            pass
        elif m_edge != m2:
            x = - (b2 - b_edge) / (m2 - m_edge)
            if (x < second_point[0]) and (x > new_point[0]):
                isIntersection = True
                break

    return isIntersection, matchPoint


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

def comparePoints(point1, point2):
    """
    Checks whether two points are the same approximately. Returns true if the two points are under the
    threshold distance, returns false otherwise
    :param point1:
    :param point2:
    :return:
    """
    return norm(point2 - point1) < 0.05