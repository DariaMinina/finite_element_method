import matplotlib.pyplot as plt
import csv

import numpy as np

plt.rcParams.update({
    'font.size': 16,
    "text.usetex": True
})

if __name__ == "__main__":

    X1 = []
    Y1 = []
    X2 = []
    Y2 = []
    X3 = []
    Y3 = []

    with open('data/analSol20.txt', 'r') as datafile:
        plotting = csv.reader(datafile, delimiter=' ')

        for ROWS in plotting:
            X1.append(float(ROWS[0]))
            Y1.append(float(ROWS[1]))

    with open('data/cubicSol20.txt', 'r') as datafile:
        plotting = csv.reader(datafile, delimiter=' ')

        for ROWS in plotting:
            X2.append(float(ROWS[0]))
            Y2.append(float(ROWS[1]))

    with open('data/linearSol20.txt', 'r') as datafile:
        plotting = csv.reader(datafile, delimiter=' ')

        for ROWS in plotting:
            X3.append(float(ROWS[0]))
            Y3.append(float(ROWS[1]))

    fig = plt.figure(figsize=(9, 12))
    ax = fig.add_subplot(1, 1, 1)
    plt.xlabel(r"$x$")
    plt.ylabel(r"$u$")
    plt.grid(True)
    major_ticks = np.arange(-70, 21, 5)
    minor_ticks = np.arange(3, 5.5, 0.5)

    ax.set_xticks(minor_ticks)
    ax.set_yticks(major_ticks)

    ax.plot(X1, Y1, label=r'analytical 20 nodes')
    ax.plot(X2, Y2, label=r'cubic 20 nodes')
    ax.plot(X3, Y3, label=r'lineal 20 nodes')
    ax.legend()

    fig.savefig('20.png')

