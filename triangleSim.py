import numpy as np
import math
import argparse
import os

from utils import *
from mathUtils import *


def arg_parse():
    """
    :return: Arguments for triangle simulation
    """
    parser = argparse.ArgumentParser(description='Isosceles Triangle Simulation')
    parser.add_argument("--angle", dest="angle", help="Primary angle of the isoscles triangle",
                        default=60, type=int)
    parser.add_argument("--number", dest="N", help="Number of triangles to simulate",
                        default=100, type=int)
    return parser.parse_args()


def main():
    args = arg_parse()
    angle = args['angle']
    N = args['N']
    triangleCoodinates = generateTriangles(angle, N)
