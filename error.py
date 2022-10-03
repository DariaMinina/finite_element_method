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

    with open('data/errors_linearSolution_Nnodes.txt', 'r') as datafile:
        plotting = csv.reader(datafile, delimiter=',')

        for ROWS in plotting:
            X1.append(float(ROWS[0]))
            Y1.append(float(ROWS[1]))

    fig = plt.figure(figsize=(20, 12))
    ax = fig.add_subplot(1, 1, 1)
    plt.xlabel(r"$N$")
    plt.ylabel(r"$abs(lin-an)$")
    plt.grid(True)
    major_ticks = np.arange(100, 1901, 200)
    minor_ticks = np.arange(0.0000001, 0.02, 0.01)

    #ax.set_yticks(minor_ticks)
    ax.set_xticks(major_ticks)

    ax.plot(X1, Y1, label=r'error')
    ax.legend()

    fig.savefig('error.png')