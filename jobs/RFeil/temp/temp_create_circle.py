import math
import numpy as np

import matplotlib.pyplot as plt

import csv

def rotate(origin, point, angle):
    """
    Rotate a point counterclockwise by a given angle around a given origin.

    The angle should be given in radians.
    """
    ox, oy = origin
    px, py = point

    qx = ox + math.cos(angle) * (px - ox) - math.sin(angle) * (py - oy)
    qy = oy + math.sin(angle) * (px - ox) + math.cos(angle) * (py - oy)
    return qx, qy



def circle_shape(x):
    # Def of circle x2 + y2 = r2
    r = np.max(x)  # define radius
    y = np.sqrt(r**2 - x**2)

    return y



if __name__ == "__main__":

    D = 2 # diameter in inches

    x = np.round(np.arange(-D/2, D/2 + 0.00001, 0.01), 2)
    print(str(x))

    y = circle_shape(x)
    print(y)
    plt.plot(x, y)
    plt.plot(x, -y)
    plt.show()

    with open('circle_shape.csv', mode='w') as egg_shape:
        egg_shape_writer = csv.writer(egg_shape, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)


        for i in range(len(x)):
            egg_shape_writer.writerow([str(-x[i]), str(y[i])])

        for i in range(len(x)):
            egg_shape_writer.writerow([str(-x[i]), str(-y[i])])

    # EOF