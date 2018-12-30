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

    :param angle: the main angle of the isosceles triangles being simulated
    :param N: The number of triangles to be generated
    :return: triangleCoordinates: the coordinates of the N acute triangle's points
    """
    angle = toRadians(angle)
    # generate the first triangle
    point1 = [0, 0]
    point2 = [2*math.sin(angle/2), 0]
    point3 = [math.sin(angle/2), math.cos(angle/2)]

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
        b = boundary.shape[1]
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
    for i in range(0,3):
        for j in range(0, N):
            plt.plot(coordinates[i][2*j], coordinates[i][2*j+1], marker='o')
    return plt


def generateIndividualTriangle(boundary, edge):
    """
    Generates an individual triangle based on the existing boundary and the edge we want to add the triangle to

    :param boundary: The existing shape's boundary
    :param edge: The edge that we're adding the triangle onto
    :return: the coordinates of the new triangle object
    """

