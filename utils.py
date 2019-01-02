import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import math
from random import randint
from triangleClass import *

from mathUtils import *


def newline(p1, p2):
    """
    This function creates a new line through p1 and p2 spanning the whole plot

    :param p1: first point
    :param p2: second point
    :return: a line that spans the whole plot through p1 and p2
    """
    ax = plt.gca()
    xmin, xmax = ax.get_xbound()

    if p2[0] == p1[0]:
        xmin = xmax = p1[0]
        ymin, ymax = ax.get_ybound()
    else:
        ymax = p1[1] + (p2[1] - p1[1]) / (p2[0] - p1[0]) * (xmax - p1[0])
        ymin = p1[1] + (p2[1] - p1[1]) / (p2[0] - p1[0]) * (xmin - p1[0])

    l = mlines.Line2D([xmin, xmax], [ymin, ymax])
    ax.add_line(l)
    return l


def generateTriangles(angle, N):
    """
    Generates N isosceles triangles in a monte-carlo fashion to create a random packing
    :param angle: the main angle of the isosceles triangles being simulated
    :param N: The number of triangles to be generated
    :return: triangleCoordinates: the coordinates of the N acute triangle's points
    """
    angle = toRadians(angle)
    # generate the first triangle
    point1 = [0, 0]
    point2 = [round3(2*math.sin(angle/2)), 0]
    point3 = [round3(math.sin(angle/2)), round3(math.cos(angle/2))]

    triangleCoordinates = np.zeros((3,2*N))
    boundary = np.zeros((3, 2))

    triangleCoordinates[0] = point1
    triangleCoordinates[1] = point2
    triangleCoordinates[2] = point3
    print(triangleCoordinates[2], point3)
    """
    firstTriangle = triangle(triangleCoordinates)

    # Generate the other N-1 triangles
    for i in range(2, N):
        # Select line on boundary on which to add new triangle, try to add new triangle
        b = boundary.shape[0]
        addedTriangle = False
        while not addedTriangle:
            edge = randint((0, b))
            triangle = generateIndividualTriangle(boundary, edge)

            if not intersection(boundary, triangle):
                addedTriangle = True


        # Update boundary
        
    """

    return triangleCoordinates


def drawTriangles(coordinates):
    """

    :param coordinates: The triples containing the coordinates of the triangles's points.
    An 3x2*N array where N is the number of triangles
    :return:
    """
    N = int(coordinates.shape[1]/2)
    print(N)
    for j in range(0, N):
        print(coordinates)
        plt.plot([coordinates[0][2*j], coordinates[1][2*j]],
                 [coordinates[0][2*j+1], coordinates[1][2*j+1]],
                 marker='o')
        plt.plot([coordinates[1][2*j], coordinates[2][2*j]],
                 [coordinates[1][2*j+1], coordinates[2][2*j+1]],
                 marker='o')
        plt.plot([coordinates[2][2*j], coordinates[0][2*j]],
                 [coordinates[2][2*j+1], coordinates[0][2*j+1]],
                 marker='o')
    return plt


def generateIndividualTriangle(boundary, edge, angle):
    """
    Generates an individual triangle based on the existing boundary and the edge we want to add the triangle to

    :param boundary: The existing shape's boundary, a list
    :param edge: The edge that we're adding the triangle onto, an int
    :param angle: the angle of the isosceles triangle
    :return: the coordinates of the new triangle object
    """
    b = boundary.shape[0]
    angle = toRadians(angle)
    point1 = boundary[edge]
    point2 = boundary[(edge+1) % b]
    AB = point2 - point1
    Midpoint = point1 + 1/2*point2
    alpha = math.acos(AB[0]/norm(AB))
    e = [round3(math.cos(math.pi/2 + alpha)), round3(math.sin(math.pi/2 + alpha))]
    magnitude = round3(1/2 * norm(AB) / math.tan(angle/2))
    v = round3(magnitude*e)

    # try first orientation
    point3 = Midpoint + v
    coordinates = np.zeros((3,2))
    coordinates[0] = point1
    coordinates[1] = point2
    coordinates[2] = point3
    newTriangle = triangle(coordinates)
    added = intersection(boundary, newTriangle, edge)

    # try second orientation
    if not added:
        point3 = Midpoint - v
        coordinates[2] = point3
        newTriangle.setCoordinates(coordinates)
        added = intersection(boundary, newTriangle, edge)

    # In the case that neither side works, return false
    if not added:
        return False
