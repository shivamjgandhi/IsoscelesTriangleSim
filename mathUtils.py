import numpy as np
import math

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

def intersection(boundary, triangle):
    """
    Tells whether or not the new triangle intersects the boundary

    :param boundary: the boundary of the shape created so far
    :param triangle: 3x2 array with the 3 points of the triangle
    :return isIntersection: tells whether or not the new triangle intersects the boundary
    """
    isIntersection = False
    b = boundary.shape[1]
    for i in range(0, b):
        print('look here: ', boundary)
        point1 = boundary[0:2, b]
        point2 = boundary[0:2, b+1]
        mBoundary = (point2[1]-point1[1])/(point2[0] - point1[0])
        bBoundary = point1[1] - mBoundary*point1[0]
        # Check if line 1 intersects
        m1 = (triangle[0][1] - triangle[1][1])/(triangle[0][0] - triangle[1][0])
        b1 = triangle[0][1] - m1*triangle[0][0]
        if m1 == mBoundary:
            print('Matching slopes')
            # initiate sequence that checks whether it's the same line or not
        else:
            x = - (bBoundary - b1)/(mBoundary - m1)
            if (x < triangle[0][0]) or (x > triangle[1][0]):
                isIntersection = True

        # Check if line 2 intersects
        m2 = (triangle[1][1] - triangle[2][1])/(triangle[1][0] - triangle[2][0])
        b2 = triangle[1][1] - m2*triangle[1][0]
        if m2 == mBoundary:
            print('Matching slope')
        else:
            x = - (bBoundary - b2)/(mBoundary - m2)
            if (x < triangle[1][0]) or (x > triangle[2][0]):
                isIntersection = True

        # Check if line 3 intersects
        m3 = (triangle[2][1] - triangle[0][1])/(triangle[2][0] - triangle[0][0])
        b3 = triangle[2][1] - m3*triangle[2][0]
        if m3 == mBoundary:
            print('Matching slope')
        else:
            x = -(bBoundary - b3)/(mBoundary - m3)
            if (x < triangle[2][0]) or (x > triangle[0][0]):
                isIntersection = True

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