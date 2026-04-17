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
from scipy.spatial.distance import pdist
from generateInstance import generateGraph, has_single_strong_component, testGraph, to_networkx
from generateInstance2 import (
    assign_required,
    build_voronoi_adjacency_edges,
    generate_vertices,
    has_single_strong_component as has_single_strong_component2,
    to_networkx as to_networkx2,
)
from computations import (
    avg_distance_to_center,
    bbox_center,
    build_distance_matrix,
    calculate_avg_dist_between_centroids,
    calculate_avg_dist_depot_to_active_cells,
    calculate_avg_internal_dist_cells,
    calculate_dist_depot_to_hottest_cell,
    calculate_node_density,
    calculate_req_edge_length_stats,
    calculate_req_node_degrees,
    calculate_sammon_error,
    compute_circuity,
    compute_mst_odd_weight,
    mean_center,
    median_center,
    sammon_mapping,
)


@dataclass
class Metrics:
    num_nodes: int
    sqrt_num_nodes: float
    num_edges: int
    vertices_param: int
    seed: str
    required_ratio: float
    convex_hull_area: float
    bbox_area: float
    bbox_perimeter: float
    convex_hull_perimeter: float
    width: float
    height: float
    avg_pairwise_dist_req: float
    avg_dist_depot_req: float
    dist_depot_bbox_center: float
    nodes_within_50p_radius: int
    nodes_within_75p_radius: int
    req_nodes_within_50p_radius: int
    req_nodes_within_75p_radius: int
    nodes_within_50p_bbox: int
    nodes_within_75p_bbox: int
    avg_dist_mean: float
    avg_dist_median: float
    avg_dist_bbox: float
    dist_depot_to_hottest_cell_10: float
    dist_depot_to_hottest_cell_15: float
    avg_dist_depot_to_active_cells_10: float
    avg_dist_depot_to_active_cells_15: float
    avg_internal_dist_cells_10: float
    avg_internal_dist_cells_15: float
    avg_dist_between_centroids_10: float
    avg_dist_between_centroids_15: float
    node_density: float
    num_req_nodes: int
    num_even_req_nodes: int
    num_odd_req_nodes: int
    circuity_avg: float
    prop_dead_ends: float
    dist_min: float
    dist_max: float
    dist_mean: float
    mst_odd_weight: float
    num_req_edges: int
    sammon_path: str
    sammon_error: float
    req_edges_mean: float
    req_edges_median: float
    req_edges_std: float
    graphml_path: str

    def to_string(self) -> str:
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


def extract_metrics(graphml_path, vertices_param, seed_value, required_ratio) -> Metrics:
    """Extract geometric and topological metrics from a generated graph instance."""
    g = nx.read_graphml(graphml_path)

    node_list = []
    positions = []
    for node, data in g.nodes(data=True):
        node_list.append(node)
        positions.append((float(data["x"]), float(data["y"])))
    points = np.array(positions)
    depot = points[0]
    depot_pos = points[0]

    hull = ConvexHull(points)
    min_x, max_x = points[:, 0].min(), points[:, 0].max()
    min_y, max_y = points[:, 1].min(), points[:, 1].max()
    val_width = max_x - min_x
    val_height = max_y - min_y

    node_ids = list(g.nodes())
    node_to_idx = {n: i for i, n in enumerate(node_ids)}

    req_indices = set()
    for u, v, data in g.edges(data=True):
        if int(data.get("required", 0)) == 1:
            req_indices.add(node_to_idx[u])
            req_indices.add(node_to_idx[v])

    dist_to_depot = np.linalg.norm(points[1:] - depot, axis=1)
    if len(dist_to_depot) > 0:
        max_radius = float(np.max(dist_to_depot))
        nodes_in_50p = int(np.sum(dist_to_depot <= 0.5 * max_radius))
        nodes_in_75p = int(np.sum(dist_to_depot <= 0.75 * max_radius))
    else:
        nodes_in_50p = 0
        nodes_in_75p = 0

    if len(req_indices) > 0:
        req_points = points[list(req_indices)]
        if len(req_indices) > 1:
            avg_pairwise_req = float(pdist(req_points).mean())
        else:
            avg_pairwise_req = 0.0
        avg_depot_req = float(np.linalg.norm(req_points - depot, axis=1).mean())
        c_req = req_points.mean(axis=0)
        dist_req_to_c = np.linalg.norm(req_points - c_req, axis=1)
        max_radius_req = float(np.max(dist_req_to_c))
        req_in_50p = int(np.sum(dist_req_to_c <= 0.5 * max_radius_req))
        req_in_75p = int(np.sum(dist_req_to_c <= 0.75 * max_radius_req))
    else:
        avg_pairwise_req = 0.0
        avg_depot_req = 0.0
        req_in_50p = 0
        req_in_75p = 0

    c_bbox = bbox_center(points)
    c_mean = mean_center(points)
    c_median = median_center(points)

    dist_depot_bbox = float(np.linalg.norm(depot - c_bbox))

    dist_to_c_bbox = np.linalg.norm(points - c_bbox, axis=1)
    if len(dist_to_c_bbox) > 0:
        max_radius_bbox = float(np.max(dist_to_c_bbox))
        nodes_in_50p_bbox = int(np.sum(dist_to_c_bbox <= 0.5 * max_radius_bbox))
        nodes_in_75p_bbox = int(np.sum(dist_to_c_bbox <= 0.75 * max_radius_bbox))
    else:
        nodes_in_50p_bbox = 0
        nodes_in_75p_bbox = 0

    D = build_distance_matrix(g, node_list)

    circuity_avg = compute_circuity(D, points)

    dead_ends = [n for n, d in g.degree() if d == 1]
    prop_dead_ends = len(dead_ends) / g.number_of_nodes()

    upper = D[np.triu_indices(len(node_list), k=1)]
    dist_min = float(upper.min()) if len(upper) > 0 else 0.0
    dist_max = float(upper.max()) if len(upper) > 0 else 0.0
    dist_mean_val = float(upper.mean()) if len(upper) > 0 else 0.0

    mst_odd = compute_mst_odd_weight(g)

    num_req_edges = sum(1 for _, _, d in g.edges(data=True) if int(d.get('required', 0)) == 1)

    sammon_coords = sammon_mapping(D)
    sammon_file = str(Path(graphml_path).with_suffix('')) + '_sammon.npy'
    np.save(sammon_file, sammon_coords)

    dist_hottest_10 = calculate_dist_depot_to_hottest_cell(g, points, depot_pos, grid_size=10)
    dist_hottest_15 = calculate_dist_depot_to_hottest_cell(g, points, depot_pos, grid_size=15)
    avg_dist_active_10 = calculate_avg_dist_depot_to_active_cells(g, points, depot_pos, grid_size=10)
    avg_dist_active_15 = calculate_avg_dist_depot_to_active_cells(g, points, depot_pos, grid_size=15)
    avg_internal_10 = calculate_avg_internal_dist_cells(points, grid_size=10)
    avg_internal_15 = calculate_avg_internal_dist_cells(points, grid_size=15)
    avg_dist_centroids_10 = calculate_avg_dist_between_centroids(points, grid_size=10)
    avg_dist_centroids_15 = calculate_avg_dist_between_centroids(points, grid_size=15)
    node_density_val = calculate_node_density(points)
    num_req_val, num_even_req_val, num_odd_req_val = calculate_req_node_degrees(g)
    sammon_val = calculate_sammon_error(g, points)
    req_mean_val, req_median_val, req_std_val = calculate_req_edge_length_stats(g)

    return Metrics(
        num_nodes=g.number_of_nodes(),
        sqrt_num_nodes=float(np.sqrt(g.number_of_nodes())),
        num_edges=g.number_of_edges(),
        vertices_param=vertices_param,
        seed=str(seed_value),
        required_ratio=required_ratio,
        convex_hull_area=hull.volume,
        bbox_area=float(val_width * val_height),
        bbox_perimeter=float(2 * (val_width + val_height)),
        convex_hull_perimeter=hull.area,
        width=val_width,
        height=val_height,
        avg_pairwise_dist_req=avg_pairwise_req,
        avg_dist_depot_req=avg_depot_req,
        dist_depot_bbox_center=dist_depot_bbox,
        nodes_within_50p_radius=nodes_in_50p,
        nodes_within_75p_radius=nodes_in_75p,
        req_nodes_within_50p_radius=req_in_50p,
        req_nodes_within_75p_radius=req_in_75p,
        nodes_within_50p_bbox=nodes_in_50p_bbox,
        nodes_within_75p_bbox=nodes_in_75p_bbox,
        avg_dist_mean=avg_distance_to_center(points, c_mean),
        avg_dist_median=avg_distance_to_center(points, c_median),
        avg_dist_bbox=avg_distance_to_center(points, c_bbox),
        dist_depot_to_hottest_cell_10=dist_hottest_10,
        dist_depot_to_hottest_cell_15=dist_hottest_15,
        avg_dist_depot_to_active_cells_10=avg_dist_active_10,
        avg_dist_depot_to_active_cells_15=avg_dist_active_15,
        avg_internal_dist_cells_10=avg_internal_10,
        avg_internal_dist_cells_15=avg_internal_15,
        avg_dist_between_centroids_10=avg_dist_centroids_10,
        avg_dist_between_centroids_15=avg_dist_centroids_15,
        node_density=node_density_val,
        num_req_nodes=num_req_val,
        num_even_req_nodes=num_even_req_val,
        num_odd_req_nodes=num_odd_req_val,
        circuity_avg=circuity_avg,
        prop_dead_ends=prop_dead_ends,
        dist_min=dist_min,
        dist_max=dist_max,
        dist_mean=dist_mean_val,
        mst_odd_weight=mst_odd,
        num_req_edges=num_req_edges,
        sammon_path=os.path.basename(sammon_file),
        sammon_error=sammon_val,
        req_edges_mean=req_mean_val,
        req_edges_median=req_median_val,
        req_edges_std=req_std_val,
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
    """Main function to generate graph instances and extract metrics."""
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
