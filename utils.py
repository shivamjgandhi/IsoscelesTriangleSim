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
    packing = randomPacking(firstTriangle, 1, [firstTriangle])

    # Generate the other N-1 triangles
    for i in range(2, N):
        packing = generateIndividualTriangle(packing, angle, method)

    return packing


def drawTriangles(coordinates):
    """
    Returns a matplotlib object with the drawing of all of the triangles
    :param coordinates: The triples containing the coordinates of the triangles's points.
    An 3x2*N array where N is the number of triangles
    :return:
    """
    # TODO
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


def generateIndividualTriangle(packing, angle, method):
    """
    Generates an individual triangle based on the existing boundary and the edge we want to add the triangle to

    :param boundary: The existing shape's boundary, a list
    :param edge: The edge that we're adding the triangle onto, an int
    :param angle: the angle of the isosceles triangle
    :return: the coordinates of the new triangle object
    """
    notAdded = True
    while notAdded:
        randomEdge = packing.generateRandomEdge(method)
        proposal1, proposal2 = packing.generateProposalCoordinates(randomEdge, packing, angle)
        intersection1 = intersection(proposal1, packing, randomEdge)
        intersection2 = intersection(proposal2, packing, randomEdge)


    # added = intersection(boundary, newTriangle, edge)
    added = True
    for i in range(0, 3):
        point = coordinates[i]
        if point not in boundary:
            boundary = np.append(boundary, [point], axis=0)

    # try second orientation
    if not added:
        point3 = Midpoint - v
        coordinates[2] = point3
        newTriangle.setCoordinates(coordinates)
        added = intersection(boundary, newTriangle)

    # In the case that neither side works, return false
    if not added:
        return False
    else:
        # add to the boundary
        return boundary

def generateProposalCoordinates(edge, packing, angle):
    """

    :param edge: the edge that we're adding onto, an int
    :param packing: the total packing, a packing object
    :param angle: the angle of the triangle, a float
    :return proposal1, proposal2: the two proposal triangles, triangle objects
    """
    boundary = packing.boundary
    angle = toRadians(angle)
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

        # second proposal
        point3 = Midpoint - v
        coordinates = np.asarray([point1, point2, point3])
        proposal2 = triangle(coordinates, None)

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