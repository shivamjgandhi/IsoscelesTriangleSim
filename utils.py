import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np
import math

def newline(p1, p2):
    """
    This function creates a new line through p1 and p2 spanning the whole plot

    :param p1: first point
    :param p2: second point
    :return: a line that spans the whole plot through p1 and p2
    """
    ax = plt.gca()
    xmin, xmax = ax.get_xbound()

    if (p2[0] == p1[0]):
        xmin = xmax = p1[0]
        ymin, ymax = ax.get_ybound()
    else:
        ymax = p1[1] + (p2[1] - p1[1])/(p2[0] - p1[0])*(xmax-p1[0])
        ymin = p1[1] + (p2[1] - p1[1])/(p2[0] - p1[0])*(xmin-p1[0])

    l = mlines.Line2D([xmin, xmax], [ymin, ymax])
    ax.add_line(l)
    return l