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


def generateTriangles(angle, N, method):
    """
    Generates N isosceles triangles in a monte-carlo fashion to create a random packing
    :param method: the method by which we're generating the triangles, a string
    :param angle: the main angle of the isosceles triangles being simulated
    :param N: The number of triangles to be generated
    :return packing: The randomPacking of the total object, a randomPacking object
    """
    angle = toRadians(angle)
    # generate the first triangle
    point1 = [0, 0]
    point2 = [round3(2 * math.sin(angle / 2)), 0]
    point3 = [round3(math.sin(angle / 2)), round3(math.cos(angle / 2))]

    triangleCoordinates = np.asarray([point1, point2, point3])
    firstTriangle = triangle(triangleCoordinates, 1)
    packing = randomPacking(firstTriangle.coordinates, 1, [firstTriangle])

    # Generate the other N-1 triangles
    for i in range(2, N):
        packing = generateIndividualTriangle(packing, angle, method)

    return packing


def drawTriangles(packing):
    """
    Returns a matplotlib object with the drawing of all of the triangles
    :param packing: The entire packing object
    :return plt: the image with the triangles drawn, a plt from matplotlib
    """
    for individualTriangle in packing.triangleList:
        coordinates = individualTriangle.coordinates
        plt.plot([coordinates[0][0], coordinates[1][0]],
                 [coordinates[0][1], coordinates[1][1]],
                 marker='o')
        plt.plot([coordinates[1][0], coordinates[2][0]],
                 [coordinates[1][1], coordinates[2][1]],
                 marker='o')
        plt.plot([coordinates[2][0], coordinates[0][0]],
                 [coordinates[2][1], coordinates[0][1]],
                 marker='o')

    return plt

def generateIndividualTriangle(packing, angle, method):
    """
    Generates an individual triangle based on the existing boundary and the edge we want to add the triangle to

    :param packing: The existing packing, a packing object
    :param method: The method by which we're growing the triangle, a string
    :param angle: the angle of the isosceles triangle, a float
    :return packing: the random packing, a packing object
    """
    notAdded = True
    while notAdded:
        randomEdge = packing.generateRandomEdge(method)
        proposals = generateProposalCoordinates(randomEdge, packing, angle)
        intersections = intersection(proposals, packing, randomEdge)
        if False in intersections:
            packing.insertTriangle(proposals[intersections.index(False)], randomEdge)
            notAdded = False
    return packing

def generateProposalCoordinates(edge, packing, angle):
    """

    :param edge: the edge that we're adding onto, an int
    :param packing: the total packing, a packing object
    :param angle: the angle of the triangle, a float
    :return proposal1, proposal2: the two proposal triangles, triangle objects
    """
    boundary = packing.boundary
    point1 = boundary[edge - 1]
    point2 = boundary[edge % len(boundary)]
    AB = point2 - point1
    if norm(AB) == 1:
        # When we're adding onto the edge opposite the angle of interest
        Midpoint = point1 + 1 / 2 * point2
        alpha = math.acos(AB[0])
        e = [round3(math.cos(math.pi / 2 + alpha)), round3(math.sin(math.pi / 2 + alpha))]
        magnitude = round3(1 / 2 / math.tan(angle / 2))
        v = round3([magnitude * elem for elem in e])

        # first proposal
        point3 = Midpoint + v
        coordinates = np.asarray([point1, point2, point3])
        proposal1 = triangle(coordinates, None)
        print(point3)

        # second proposal
        point3 = Midpoint - v
        coordinates = np.asarray([point1, point2, point3])
        proposal2 = triangle(coordinates, None)
        print(point3)

        print('look here: ', boundary, v, magnitude, edge, AB, alpha)

        proposals = [proposal1, proposal2]
    else:
        # When we're adding onto one of the isosceles sides
        if AB[0] == 0:
            gamma = math.pi/2
        else:
            gamma = math.atan(AB[1]/AB[0])
        theta1 = math.pi - angle + gamma
        theta2 = math.pi - angle - gamma
        C = [[AB + norm(AB) * [math.cos(theta1), math.sin(theta1)],
             AB + norm(AB) * [math.cos(theta1), -math.sin(theta1)],
             AB + norm(AB) * [-math.cos(theta2), math.sin(theta2)],
             AB + norm(AB) * [-math.cos(theta2), -math.sin(theta2)]]]

        proposals = []
        for i in range(0, 4):
            coordinates = np.asarray([point1, point2, C[i]])
            newTriangle = triangle(coordinates, None)
            proposals.append(newTriangle)

    return proposals