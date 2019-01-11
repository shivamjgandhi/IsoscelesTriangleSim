import numpy as np
from random import randint

def uniformDist(boundary):
    """
    Gives each edge on the boundary a uniform chance of being selected to have growth off and returns an edge
    :param boundary: The boundary of tbe random packing, a list
    :return edge: The edge to grow off of, an int
    """
    return randint(0, len(boundary)-1)
