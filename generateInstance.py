# The Prize Collecting Steiner Tree Problem: Theory and Practice'' (D.S. Johnson, M. Minkoff, S. Phillips)
import operator
import argparse
from dataclasses import dataclass
from random import random, seed
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx


@dataclass(order=True)
class vertex:
    x: float
    y: float


@dataclass(order=True)
class edge:
    v1: int
    v2: int
    length: float
    required: bool


def on_segment(p, q, r):
    if (
        r[0] <= max(p[0], q[0])
        and r[0] >= min(p[0], q[0])
        and r[1] <= max(p[1], q[1])
        and r[1] >= min(p[1], q[1])
    ):
        return True
    return False


def orientation(p, q, r):
    val = ((q[1] - p[1]) * (r[0] - q[0])) - ((q[0] - p[0]) * (r[1] - q[1]))
    if val == 0:
        return 0
    return 1 if val > 0 else -1


def intersects(seg1, seg2):
    p1, q1 = seg1
    p2, q2 = seg2
    o1 = orientation(p1, q1, p2)
    o2 = orientation(p1, q1, q2)
    o3 = orientation(p2, q2, p1)
    o4 = orientation(p2, q2, q1)
    if o1 != o2 and o3 != o4:
        return True
    if o1 == 0 and on_segment(p1, q1, p2):
        return True
    if o2 == 0 and on_segment(p1, q1, q2):
        return True
    if o3 == 0 and on_segment(p2, q2, p1):
        return True
    if o4 == 0 and on_segment(p2, q2, q1):
        return True
    return False


def generateGraph(numV, requireds):
    vertices = []
    edges = []
    inUse = []

    for _ in range(numV):
        vertices.append(vertex(random(), random()))

    criticalLength = 1.6 / pow(numV, 0.5)
    for i in range(numV):
        for j in range(i + 1, numV):
            length = pow(
                pow(vertices[i].x - vertices[j].x, 2)
                + pow(vertices[i].y - vertices[j].y, 2),
                0.5,
            )
            if length <= criticalLength:
                edges.append(edge(i, j, length, True))

    edges.sort(key=operator.attrgetter("length"))
    for e in edges:
        notCrossed = True
        for i in inUse:
            if e.v1 != i.v1 and e.v1 != i.v2 and e.v2 != i.v1 and e.v2 != i.v2:
                myRet = intersects(
                    ((vertices[e.v1].x, vertices[e.v1].y), (vertices[e.v2].x, vertices[e.v2].y)),
                    ((vertices[i.v1].x, vertices[i.v1].y), (vertices[i.v2].x, vertices[i.v2].y)),
                )
                if myRet:
                    notCrossed = False
                    break
        if notCrossed:
            if random() > requireds:
                e.required = False
            inUse.append(e)

    return vertices, inUse


def testGraph(vertices, edges):
    connected = []
    nV = len(vertices)
    for i in range(nV):
        connected.append(i)
    for e in edges:
        if connected[e.v1] != connected[e.v2]:
            eMin = min(connected[e.v1], connected[e.v2])
            eMax = max(connected[e.v1], connected[e.v2])
            for i in range(nV):
                if connected[i] == eMax:
                    connected[i] = eMin
    for i in range(nV):
        if connected[i] != 0:
            return "No"

    connected = []
    numGroups = 0
    for i in range(nV):
        connected.append(i)
    for e in edges:
        if e.required and connected[e.v1] != connected[e.v2]:
            eMin = min(connected[e.v1], connected[e.v2])
            eMax = max(connected[e.v1], connected[e.v2])
            for i in range(nV):
                if connected[i] == eMax:
                    connected[i] = eMin
    for i in range(nV):
        counter = 0
        for ii in range(nV):
            if connected[ii] == i:
                counter = counter + 1
                if counter > 1:
                    break
        if counter > 1:
            numGroups = numGroups + 1
        if numGroups > 1:
            break
    if numGroups < 1:
        return "No"
    if numGroups == 1:
        return "Chinese"
    return "Rural"


def to_networkx(vertices, edges):
    g = nx.Graph()
    for i, v in enumerate(vertices):
        g.add_node(i, x=v.x, y=v.y)
    for e in edges:
        g.add_edge(e.v1, e.v2, length=e.length, required=int(e.required))
    return g


def has_single_strong_component(g):
    directed = g.to_directed()
    return nx.number_strongly_connected_components(directed) == 1


def plot_graph(vertices, edges, output_path):
    x_values = [v.x for v in vertices]
    y_values = [v.y for v in vertices]

    for e in edges:
        color = "red" if e.required else "black"
        width = 2 if e.required else 1
        plt.plot(
            (vertices[e.v1].x, vertices[e.v2].x),
            (vertices[e.v1].y, vertices[e.v2].y),
            color=color,
            linewidth=width,
        )

    plt.scatter(x_values, y_values, color="black", s=10)
    plt.savefig(output_path)
    plt.close()


def parse_args():
    parser = argparse.ArgumentParser(description="Generate a graph instance and export it as GraphML.")
    parser.add_argument("--graph", required=True, help="Output GraphML path.")
    parser.add_argument(
        "--required",
        required=True,
        type=float,
        help="Required-edge percentage (0-100) or ratio (0-1).",
    )
    parser.add_argument("--vertices", required=True, type=int, help="Number of vertices.")
    parser.add_argument("--seed", default="0", help="Random seed.")
    parser.add_argument(
        "--plot",
        nargs="?",
        const="",
        default=None,
        help="Optional PNG output path. If used without a value, uses <graph_basename>.png.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    output_path = args.graph
    required_percent = args.required
    numV = args.vertices
    seed_value = args.seed

    requireds = required_percent / 100.0 if required_percent > 1.0 else required_percent
    requireds = max(0.0, min(1.0, requireds))

    seed(seed_value)

    while True:
        vertices, edges = generateGraph(numV, requireds)
        if testGraph(vertices, edges) != "No":
            candidate = to_networkx(vertices, edges)
            if has_single_strong_component(candidate):
                g = candidate
                break

    if not has_single_strong_component(g):
        raise RuntimeError("Generated graph does not have a single strongly connected component.")

    nx.write_graphml(g, output_path)

    if args.plot is not None:
        if args.plot == "":
            plot_output = str(Path(output_path).with_suffix(".png"))
        else:
            plot_output = args.plot
        plot_graph(vertices, edges, plot_output)


if __name__ == "__main__":
    main()
