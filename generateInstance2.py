import argparse
from dataclasses import dataclass
from pathlib import Path
from random import random, seed

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import networkx as nx
import numpy as np


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


def generate_vertices(num_vertices):
    return [vertex(random(), random()) for _ in range(num_vertices)]


def euclidean_distance(v1, v2):
    dx = v1.x - v2.x
    dy = v1.y - v2.y
    return (dx * dx + dy * dy) ** 0.5


def build_voronoi_adjacency_edges(vertices):
    # Two seed points share a Voronoi edge iff they are adjacent in the Delaunay triangulation.
    x = np.array([v.x for v in vertices])
    y = np.array([v.y for v in vertices])
    triangulation = mtri.Triangulation(x, y)

    adjacency = set()
    if triangulation.triangles is None:
        return []

    for triangle in triangulation.triangles:
        a, b, c = int(triangle[0]), int(triangle[1]), int(triangle[2])
        adjacency.add(tuple(sorted((a, b))))
        adjacency.add(tuple(sorted((a, c))))
        adjacency.add(tuple(sorted((b, c))))

    edges = []
    for u, v in sorted(adjacency):
        dist = euclidean_distance(vertices[u], vertices[v])
        edges.append(edge(u, v, dist, True))
    return edges


def assign_required(edges, required_ratio):
    for e in edges:
        e.required = random() <= required_ratio


def to_networkx(vertices, edges):
    g = nx.Graph()
    for i, v in enumerate(vertices):
        g.add_node(i, x=v.x, y=v.y)
    for e in edges:
        g.add_edge(e.v1, e.v2, length=e.length, required=int(e.required))
    return g


def has_single_strong_component(g):
    return nx.number_strongly_connected_components(g.to_directed()) == 1


def plot_graph(vertices, edges, output_path):
    for e in edges:
        color = "red" if e.required else "black"
        width = 2 if e.required else 1
        plt.plot(
            (vertices[e.v1].x, vertices[e.v2].x),
            (vertices[e.v1].y, vertices[e.v2].y),
            color=color,
            linewidth=width,
        )

    plt.scatter([v.x for v in vertices], [v.y for v in vertices], color="black", s=10)
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.gca().set_aspect("equal", adjustable="box")
    plt.savefig(output_path)
    plt.close()


def parse_args():
    parser = argparse.ArgumentParser(description="Generate Voronoi-based graph instance and export as GraphML.")
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

    requireds = args.required / 100.0 if args.required > 1.0 else args.required
    requireds = max(0.0, min(1.0, requireds))

    if args.vertices < 3:
        raise ValueError("--vertices must be >= 3 to build Voronoi adjacency.")

    seed(args.seed)

    while True:
        vertices = generate_vertices(args.vertices)
        edges = build_voronoi_adjacency_edges(vertices)
        assign_required(edges, requireds)
        g = to_networkx(vertices, edges)
        if has_single_strong_component(g):
            break

    nx.write_graphml(g, args.graph)

    if args.plot is not None:
        plot_output = str(Path(args.graph).with_suffix(".png")) if args.plot == "" else args.plot
        plot_graph(vertices, edges, plot_output)


if __name__ == "__main__":
    main()
