import matplotlib.pyplot as plt
from glob import glob
import os
import numpy as np
from matplotlib.collections import LineCollection

def read_edges(filename):
    edges = set()
    X = []
    Y = []
    print("Starting reading edges and points")
    with open(filename, 'r') as file:
        lines = file.readlines()
        a = lines[0].split()
        N = int(a[0])
        N_pbc = int(a[1])
        for i in range(1, N_pbc + 1):
            x, y = map(float, lines[i].split())
            X.append(x)
            Y.append(y)
        for i in range(N_pbc + 2, len(lines), 2):
            x1, y1 = map(float, lines[i - 1].split())
            x2, y2 = map(float, lines[i].split())
            edge = tuple(sorted([(x1, y1), (x2, y2)]))
            edges.add(edge)
    print("Finished reading edges and points")
    return N, X, Y, edges

def plot_edges(N, X, Y, edges):
    lines = [ [edge[0], edge[1]] for edge in edges ]
    lc = LineCollection(lines, colors='b', linewidths=1)

    # Create figure and axis
    fig, ax = plt.subplots()
    ax.add_collection(lc)
    plt.scatter(X[:N], Y[:N], c = np.linspace(0, len(X[:N]), len(X[:N])))
    plt.scatter(X[N:], Y[N:], c = "black")
    plt.axis("off")
    plt.show()

if __name__ == "__main__":
    if 1:
        N, X, Y, edges = read_edges('dump/a.txt')
        plot_edges(N, X, Y, edges)
        
    else:
        liste = glob("dump/*.txt")
        for i in liste:
            X, Y, edges = read_edges(i)
            plot_edges(X, Y, edges)
            plt.savefig(i.replace("dump", "images").replace(".txt", ".png"))
            plt.close()
