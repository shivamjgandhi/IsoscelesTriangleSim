import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import math
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

    triangleCoordinates[0, 0:2], boundary[0, 0:2] = point1
    triangleCoordinates[1, 0:2], boundary[0, 0:2] = point2
    triangleCoordinates[2, 0:2], boundary[0, 0:2] = point3

    # Generate the other N-1 triangles
    for i in range(2, N):
        # Select line on boundary on which to add new triangle, try to add new triangle
        [a, b] = boundary.shape
        addedTriangle = False
        while not addedTriangle:

            if not intersection(boundary, triangle):
                addedTriangle = True

        # Update boundary

    return triangleCoordinates


def drawTriangeles(coordinates):
    """

    :param coordinates: The triples containing the coordinates of the triangles's points
    :return:
    """
