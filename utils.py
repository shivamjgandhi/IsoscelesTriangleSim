import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import math
from math import cos, sin, acos, asin
import random
from random import randint
from triangleClass import *
from numpy.linalg import norm

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
    point2 = [1, 0]
    point3 = [1/2, 0.5/round3(math.tan(angle / 2))]

    triangleCoordinates = np.asarray([point1, point2, point3])
    firstTriangle = triangle(triangleCoordinates, 1)
    packing = randomPacking(firstTriangle.coordinates, 1, [firstTriangle])
    packing.packingCenter = firstTriangle.center
    packing.computeNewRadiusGyration(firstTriangle)
    print('added triangle 1')

    # Generate the other N-1 triangles
    for i in range(1, N):
        print('added triangle ', str(i))
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

    """
    notAdded = True
    while notAdded:
        randomEdge = packing.generateRandomEdge(method)
        proposal = generateProposalCoordinates(randomEdge, packing, angle)
        intersections = intersection(proposal, packing)
        print(proposal[0].coordinates)
        if False in intersections:
            packing.insertTriangle(proposal[0][intersections.index(False)], randomEdge)
            notAdded = False
    return packing
    """
    notAdded = True
    proposal = None
    growthEdge = None
    alpha = 1
    while notAdded:
        randomEdge = packing.generateRandomEdge(method='uniform')
        proposal = generateUniformProposal(randomEdge, packing, angle)
        """
        if method == 'uniform':
            proposal = generateUniformProposal(randomEdge, packing, angle)
        elif method == 'proposals':
            # in this case, randomEdge is between 0 and len(boundarydist)
            growthEdge = packing.boundary[packing.boundaryDist[randomEdge]]
            proposal = generateSingleProposal(randomEdge, growthEdge, packing, angle)
        """
        isIntersection, matchPoint = boundaryIntersection(packing, proposal, randomEdge)
        deltaRadiusGyration = np.max(packing.computeNewRadiusGyration(proposal) - packing.radiusOfGyration, 0)
        t = random.uniform(0, 1)
        if math.exp(-alpha * deltaRadiusGyration) >= t:
            # Add in the new triangle
            notAdded = False
            packing.updatePacking(proposal, randomEdge, matchPoint)
        """
        if not isIntersection:
            deltaRadiusGyration = np.max(packing.computeNewRadiusGyration(proposal) - packing.radiusOfGyration, 0)
            t = random.uniform(0, 1)
            if math.exp(-alpha*deltaRadiusGyration) >= t:
                # Add in the new triangle
                notAdded = False
                packing.updatePacking(proposal, randomEdge, matchPoint)
        """

    return packing

def generateProposalCoordinates(edge, packing, angle):
    """

    :param edge: the edge that we're adding onto, an int
    :param packing: the total packing, a packing object
    :param angle: the angle of the triangle, a float
    :return proposal1, proposal2: the two proposal triangles, triangle objects
    """
    boundary = packing.boundary
    point1 = boundary[edge]
    point2 = boundary[(edge+1) % len(boundary)]
    AB = point2 - point1
    if norm(AB) == 1:
        # When we're adding onto the edge opposite the angle of interest
        Midpoint = point1 + 1 / 2 * AB
        e = [-AB[1], AB[0]]
        magnitude = round3(1 / 2 / math.tan(angle / 2))
        v = round3([magnitude * elem for elem in e])

        # proposal coordinates
        point3 = Midpoint - v
        coordinates = np.asarray([point1, point2, point3])
        proposal = triangle(coordinates, None)

        proposals = [proposal]
    else:
        # When we're adding onto one of the isosceles sides
        if AB[0] == 0:
            gamma = math.pi/2
        else:
            gamma = math.atan(AB[1]/AB[0])
        theta1 = math.pi - angle + gamma
        theta2 = math.pi - angle - gamma
        C = [AB + floatMult(norm(AB), [round3(math.cos(theta1)), round3(math.sin(theta1))]),
             AB + floatMult(norm(AB), [round3(math.cos(theta1)), -round3(math.sin(theta1))]),
             AB + floatMult(norm(AB), [-round3(math.cos(theta2)), round3(math.sin(theta2))]),
             AB + floatMult(norm(AB), [-round3(math.cos(theta2)), -round3(math.sin(theta2))])]
        print('look here: ', AB, [round3(math.cos(theta1)), round3(math.sin(theta1))])
        proposals = []
        for i in range(4):
            coordinates = np.asarray([point1, point2, C[i]])
            newTriangle = triangle(coordinates, None)
            proposals.append(newTriangle)

    return proposals

def generateSingleProposal(proposal_edge, edge, packing, angle):
    """
    This function generates a single Proposal triangle when using the weighted distribution method

    :param proposal_edge: the proposal triangle index in packing.dist, an int
    :param edge: the edge of growth, an int
    :param packing: the random packing, a packing object
    :param angle: the angle of the triangle, a float
    :return proposal_triangle: the proposal triangle
    """
    point1 = packing.boundary[edge]
    point2 = packing.boundary[(edge + 1) % len(packing.boundary)]
    AB = point2 - point1

    # First we check which orientation triangle we're creating based on the boundary dist, and then create
    if packing.boundaryDist[proposal_edge] == packing.boundaryDist[(proposal_edge + 1) % len(packing.boundaryDist)]:
        # If this is the case, then we have left orientation and we make the new point to be closer to the edge_point
        # Compute theta
        theta = acos(point1[0]/norm(point1))
        # Compute rotation matrix R
        c, s = np.cos(theta), np.sin(theta)
        R = np.array(((c, -s), (s, c)))
        # Compute v'
        v_prime = np.asarray([-sin(toRadians(90-angle/2)), cos(toRadians(90-angle/2))])
        # Compute v
        v = R.dot(v_prime)
        # Compute the new_point and proposal triangle
        new_point = point1 + v
        proposal_triangle = np.asarray([point1, point2, new_point])

    elif packing.boundaryDist[proposal_edge] == packing.boundaryDist[(proposal_edge - 1) % len(packing.boundaryDist)]:
        # If this is the case, then we have right orientation and we make the new point to be farther from the
        # edge_point
        theta = acos(point1[0] / norm(point1))
        # Compute rotation matrix R1, R2
        c1, s1 = np.cos(theta), np.sin(theta)
        c2, s2 = np.cos(-toRadians(angle), -toRadians(angle))
        R1 = np.array(((c1, -s1), (s1, c1)))
        R2 = np.array(((c2, -s2), (s2, c2)))
        # Compute v'
        v_prime = np.asarray([-sin(toRadians(90 - angle / 2)), cos(toRadians(90 - angle / 2))])
        # Compute v
        v = R2.dot(R1.dot(v_prime))
        # Compute the new_point and proposal triangle
        new_point = point2 + v
        proposal_triangle = np.asarray([point1, point2, new_point])
    else:
        # When we're adding onto the edge opposite the angle of interest
        Midpoint = point1 + 1 / 2 * AB
        e = [-AB[1], AB[0]]
        magnitude = round3(1 / 2 / math.tan(angle / 2))
        v = round3([magnitude * elem for elem in e])

        # proposal coordinates
        point3 = Midpoint - v
        coordinates = np.asarray([point1, point2, point3])
        proposal_triangle = triangle(coordinates, None)

    return proposal_triangle

def generateUniformProposal(edge, packing, angle):
    point1 = packing.boundary[edge]
    point2 = packing.boundary[(edge + 1) % len(packing.boundary)]
    AB = point2 - point1
    length = norm(AB)
    mode = None
    if (length > 1.02) or (length < 0.98):
        t = random.uniform(0, 1)
        if t >= 0.5:
            mode = "right"
        else:
            mode = "left"

    scaling_value = 1

    # First we check which orientation triangle we're creating based on the boundary dist, and then create
    if mode == "left":
        print('type 1')
        # If this is the case, then we have left orientation and we make the new point to be closer to the edge_point
        # Compute theta, compute v, rotate v by theta, and add to the first point
        theta = recoverAngle(AB)
        v = scaling_value*np.asarray([cos(math.pi/2 - angle/2), -sin(math.pi/2 - angle/2)])
        R = np.asarray([[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]])
        new_point = point1 + np.matmul(R, v)
        proposal_triangle = triangle(np.asarray([point1, point2, new_point]), None)

    elif mode == "right":
        print('type 2')
        # If this is the case, then we have right orientation and we make the new point to be farther from the
        # edge_point
        # Compute theta, compute v, rotate v by theta, and add to the second point
        theta = recoverAngle(AB)
        v = scaling_value*np.asarray([-cos(math.pi/2 - angle/2), -sin(math.pi/2 - angle/2)])
        R = np.asarray([[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]])
        new_point = point2 + np.matmul(R, v)
        proposal_triangle = triangle(np.asarray([point1, point2, new_point]), None)

    else:
        print('type 3')
        # When we're adding onto the edge opposite the angle of interest
        Midpoint = point1 + 1 / 2 * AB
        e = [-AB[1], AB[0]]
        magnitude = round3(1 / 2 / math.tan(angle / 2))
        v = round3([magnitude * elem for elem in e])

        # proposal coordinates
        point3 = Midpoint - v
        coordinates = np.asarray([point1, point2, point3])
        proposal_triangle = triangle(coordinates, None)

    # print('boundary: ', packing.boundary)
    # print('proposal triangle: ', proposal_triangle.coordinates)
    return proposal_triangle

def recoverAngle(vector):
    """
    Returns the angle the vector makes to the vector [1, 0]
    :param vector: The vector we want to find the angle of
    :return theta: The angle that the vector makes with the horizontal [1, 0] in radians
    """
    theta = acos(vector[0]/norm(vector))
    if vector[1] < 0:
        theta = 2*math.pi - theta
    return theta