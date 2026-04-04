import argparse
import csv
import os
from dataclasses import dataclass, fields
from itertools import product
from pathlib import Path
from random import seed as set_seed
import networkx as nx
import numpy as np
from scipy.spatial import ConvexHull
from generateInstance import generateGraph, has_single_strong_component, testGraph, to_networkx
from generateInstance2 import (
    assign_required,
    build_voronoi_adjacency_edges,
    generate_vertices,
    has_single_strong_component as has_single_strong_component2,
    to_networkx as to_networkx2,
)


@dataclass
class Metrics:
    num_nodes: int
    num_edges: int
    vertices_param: int
    seed: str
    required_ratio: float
    convex_hull_area: float
    convex_hull_perimeter: float
    width: float
    height: float
    avg_dist_mean: float
    avg_dist_median: float
    avg_dist_bbox: float
    graphml_path: str

    def to_string(self):
        values = [str(getattr(self, f.name)) for f in fields(self)]
        return ";".join(values)


def generate_instance_1(num_vertices, required_ratio, seed_value) -> nx.Graph:
    """Generate a graph using the original proximity-based method. Same as generateInstance.py"""
    set_seed(seed_value)
    while True:
        vertices, edges = generateGraph(num_vertices, required_ratio)
        if testGraph(vertices, edges) != "No":
            candidate = to_networkx(vertices, edges)
            if has_single_strong_component(candidate):
                return candidate


def generate_instance_2(num_vertices, required_ratio, seed_value) -> nx.Graph:
    """Generate a graph using the Voronoi/Delaunay-based method. Same as generateInstance2.py"""
    set_seed(seed_value)
    while True:
        vertices = generate_vertices(num_vertices)
        edges = build_voronoi_adjacency_edges(vertices)
        assign_required(edges, required_ratio)
        g = to_networkx2(vertices, edges)
        if has_single_strong_component2(g):
            return g


def mean_center(points: np.ndarray) -> np.ndarray:
    """Center as the mean of all node positions."""
    return points.mean(axis=0)


def median_center(points: np.ndarray) -> np.ndarray:
    """Center as the median of all node positions."""
    return np.median(points, axis=0)


def bbox_center(points: np.ndarray) -> np.ndarray:
    """Center as the midpoint of the bounding box."""
    return (points.min(axis=0) + points.max(axis=0)) / 2.0


def avg_distance_to_center(points: np.ndarray, center: np.ndarray) -> float:
    """Average Euclidean distance from all nodes to a given center point."""
    diffs = points - center
    distances = np.sqrt((diffs ** 2).sum(axis=1))
    return float(distances.mean())


def extract_metrics(graphml_path, vertices_param, seed_value, required_ratio) -> Metrics:
    """Extract geometric metrics from a generated graph instance."""
    g = nx.read_graphml(graphml_path)

    positions = []
    for _, data in g.nodes(data=True):
        positions.append((float(data["x"]), float(data["y"])))
    points = np.array(positions)

    hull = ConvexHull(points)

    min_x, max_x = points[:, 0].min(), points[:, 0].max()
    min_y, max_y = points[:, 1].min(), points[:, 1].max()

    return Metrics(
        num_nodes=g.number_of_nodes(),
        num_edges=g.number_of_edges(),
        vertices_param=vertices_param,
        seed=str(seed_value),
        required_ratio=required_ratio,
        convex_hull_area=hull.volume,
        convex_hull_perimeter=hull.area,
        width=max_x - min_x,
        height=max_y - min_y,
        avg_dist_mean=avg_distance_to_center(points, mean_center(points)),
        avg_dist_median=avg_distance_to_center(points, median_center(points)),
        avg_dist_bbox=avg_distance_to_center(points, bbox_center(points)),
        graphml_path=os.path.basename(graphml_path),
    )

def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for graph generation and analysis."""
    parser = argparse.ArgumentParser(description="Generate multiple graph instances and extract geometric metrics.")
    parser.add_argument("--vertices", nargs="+", type=int, required=True, help="List of vertex counts.")
    parser.add_argument("--seeds", nargs="+", default=["0"], help="List of random seeds.")
    parser.add_argument("--required", nargs="+", type=float, required=True, help="List of required-edge percentages (0-100).")
    parser.add_argument("--generator", type=int, choices=[1, 2], default=1,
                        help="1 = proximity/planar (generateInstance), 2 = Delaunay (generateInstance2).")
    parser.add_argument("--output-dir", default="instances", help="Directory for generated GraphML files.")
    parser.add_argument("--csv", default=None, help="Output CSV path for metrics.")
    return parser.parse_args()


def main() -> None:
    """Main function to generate graph instances and extract metrics.
    Missing: plotting functionality.
    Missing features
    """
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    ratios = []
    for r in args.required:
        ratio = r / 100.0 if r > 1.0 else r
        ratios.append(max(0.0, min(1.0, ratio)))

    results = []
    combos = list(product(args.vertices, args.seeds, ratios))
    total = len(combos)

    for i, (nv, sd, req) in enumerate(combos, 1):
        gen_label = "delaunay" if args.generator == 2 else "planar"
        filename = f"graph_{gen_label}_v{nv}_s{sd}_r{int(req * 100)}.graphml"
        graphml_path = output_dir / filename

        print(f"[{i}/{total}] Generating ({gen_label}) v={nv}, seed={sd}, required={req:.0%} ... ", end="", flush=True)
        generate_fn = generate_instance_2 if args.generator == 2 else generate_instance_1
        g = generate_fn(nv, req, sd)
        nx.write_graphml(g, str(graphml_path))

        metrics = extract_metrics(str(graphml_path), nv, sd, req)
        results.append(metrics)
        print(metrics.to_string())

    print(f"\nGenerated {len(results)} instances in '{output_dir}/'")

    header = ";".join(f.name for f in fields(Metrics))
    print(f"\n{header}")
    for m in results:
        print(m.to_string())

    if args.csv:
        field_names = [f.name for f in fields(Metrics)]
        with open(args.csv, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=field_names, delimiter=";")
            writer.writeheader()
            for m in results:
                writer.writerow({fn: getattr(m, fn) for fn in field_names})
        print(f"\nMetrics exported to '{args.csv}'")


if __name__ == "__main__":
    main()
