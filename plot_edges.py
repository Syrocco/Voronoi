import matplotlib.pyplot as plt

def read_edges(filename):
    edges = set()
    X = []
    Y = []
    print("Starting reading edges and points")
    with open(filename, 'r') as file:
        lines = file.readlines()
        N = int(lines[0])
        for i in range(1, N + 1):
            x, y = map(float, lines[i].split())
            X.append(x)
            Y.append(y)
        for i in range(N + 2, len(lines), 2):
            x1, y1 = map(float, lines[i - 1].split())
            x2, y2 = map(float, lines[i].split())
            edge = tuple(sorted([(x1, y1), (x2, y2)]))
            edges.add(edge)
    print("Finished reading edges and points")
    return X, Y, edges

def plot_edges(X, Y, edges):
    plt.figure()
    for edge in edges:
        x_values = [edge[0][0], edge[1][0]]
        y_values = [edge[0][1], edge[1][1]]
        plt.plot(x_values, y_values, 'b-')
    plt.scatter(X, Y, color = "black")
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Voronoi Diagram Edges')
    plt.show()

if __name__ == "__main__":
    X, Y, edges = read_edges('edges.txt')
    plot_edges(X, Y, edges)