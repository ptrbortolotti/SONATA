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



def egg_shape(x):


    y = np.sqrt((144 - 9*x**2) / (17 + 2*x))


    return y



if __name__ == "__main__":

    # origin = [0.0, 0.0]
    # point = [0.24, -0.12]
    # rot_angle = 30 * np.pi /180
    #
    # qx, qy = rotate(origin, point, rot_angle)

    x = np.round(np.arange(-4, 4.00001, 0.05), 2)
    print(str(x))

    y = egg_shape(x)
    print(y)
    plt.plot(x, y)
    plt.plot(x, -y)
    plt.show()

    with open('egg_shape.csv', mode='w') as egg_shape:
        egg_shape_writer = csv.writer(egg_shape, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)


        for i in range(len(x)):
            egg_shape_writer.writerow([str(-x[i]), str(y[i])])

        for i in range(len(x)):
            egg_shape_writer.writerow([str(-x[i]), str(-y[i])])

    # EOF