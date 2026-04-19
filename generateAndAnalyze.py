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
from computations import (
    avg_distance_to_center,
    bbox_center,
    build_distance_matrix,
    calculate_avg_dist_between_centroids,
    calculate_avg_dist_depot_to_active_cells,
    calculate_avg_dist_point_to_centroids,
    calculate_avg_internal_dist_cells,
    calculate_dist_depot_to_hottest_cell,
    calculate_distance_matrix_stats,
    calculate_edge_centroid_stats,
    calculate_edge_density,
    calculate_min_depot_to_req,
    calculate_node_density,
    calculate_nodes_within_radius_fractions,
    calculate_prop_dead_ends,
    calculate_prop_odd_req_nodes,
    calculate_req_bbox_area,
    calculate_req_convex_hull_area,
    calculate_req_distance_stats,
    calculate_req_edge_length_stats,
    calculate_req_node_degrees,
    calculate_req_spatial_stats,
    calculate_sammon_error,
    calculate_std_bearings,
    calculate_std_distances,
    calculate_std_dist_odd_nodes,
    calculate_xy_spread,
    compute_circuity,
    compute_mst_odd_weight,
    get_edge_centroid_distance_matrices,
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
    # Geometry
    convex_hull_area: float
    bbox_area: float
    req_bbox_area: float
    req_convex_hull_area: float
    bbox_perimeter: float
    convex_hull_perimeter: float
    width: float
    height: float
    # Node/edge spatial distances
    avg_pairwise_dist_req: float
    avg_dist_depot_req: float
    avg_pairwise_dist_req_edges: float
    avg_dist_depot_req_edges: float
    avg_dist_req_edges_bbox_center: float
    dist_depot_bbox_center: float
    dist_depot_to_node_centroid: float
    avg_dist_bbox_center_to_req_nodes: float
    avg_dist_node_centroid_to_req_nodes: float
    # Radius counts
    nodes_within_50p_radius: int
    nodes_within_75p_radius: int
    req_nodes_within_50p_radius: int
    req_nodes_within_75p_radius: int
    nodes_within_50p_bbox: int
    nodes_within_75p_bbox: int
    # Center distances
    avg_dist_mean: float
    avg_dist_median: float
    avg_dist_bbox: float
    # Bearing / spread stats
    std_bearing_depot: float
    std_bearing_node_centroid: float
    std_bearing_bbox_centroid: float
    avg_std_xy: float
    prod_std_xy: float
    std_dist_depot: float
    std_dist_node_centroid: float
    std_dist_bbox_centroid: float
    # Grid (10x10 and 15x15)
    dist_depot_to_hottest_cell_10: float
    dist_depot_to_hottest_cell_15: float
    avg_dist_depot_to_active_cells_10: float
    avg_dist_depot_to_active_cells_15: float
    avg_internal_dist_cells_10: float
    avg_internal_dist_cells_15: float
    avg_dist_between_centroids_10: float
    avg_dist_between_centroids_15: float
    # Density and degree
    node_density: float
    edge_density: float
    num_req_nodes: int
    num_even_req_nodes: int
    num_odd_req_nodes: int
    circuity_avg: float
    prop_dead_ends: float
    prop_odd_req_nodes: float
    # Full distance matrix stats
    dist_min: float
    dist_max: float
    dist_mean: float
    dist_total_median: float
    dist_total_std: float
    # Required subset distance stats
    dist_req_min: float
    dist_req_max: float
    dist_req_mean: float
    dist_req_median: float
    dist_req_std: float
    # Depot and MST
    min_dist_depot_to_req: float
    mst_odd_weight: float
    std_dist_odd_nodes: float
    num_req_edges: int
    # Sammon
    sammon_error: float
    sammon_layout_error: float
    # Required edge length stats
    req_edges_mean: float
    req_edges_median: float
    req_edges_std: float
    # Edge centroid stats
    avg_dist_all_edge_centroids: float
    avg_dist_req_edge_centroids: float
    std_dist_req_edge_centroids: float
    graphml_path: str

    def to_string(self) -> str:
        values = [str(getattr(self, f.name)) for f in fields(self)]
        return ";".join(values)


def generate_instance_1(num_vertices, required_ratio, seed_value) -> nx.Graph:
    """Genera un grafo usando el método de proximidad original."""
    set_seed(seed_value)
    while True:
        vertices, edges = generateGraph(num_vertices, required_ratio)
        if testGraph(vertices, edges) != "No":
            candidate = to_networkx(vertices, edges)
            if has_single_strong_component(candidate):
                return candidate


def generate_instance_2(num_vertices, required_ratio, seed_value) -> nx.Graph:
    """Genera un grafo usando el método Voronoi/Delaunay."""
    set_seed(seed_value)
    while True:
        vertices = generate_vertices(num_vertices)
        edges = build_voronoi_adjacency_edges(vertices)
        assign_required(edges, required_ratio)
        g = to_networkx2(vertices, edges)
        if has_single_strong_component2(g):
            return g


def extract_metrics(graphml_path, vertices_param, seed_value, required_ratio) -> Metrics:
    """Extrae métricas geométricas y topológicas de una instancia GraphML."""
    g = nx.read_graphml(graphml_path)

    # Node positions
    node_list = []
    positions = []
    for node, data in g.nodes(data=True):
        node_list.append(node)
        positions.append((float(data["x"]), float(data["y"])))
    points = np.array(positions)
    depot = points[0]
    depot_pos = points[0]

    # Required node indices
    node_ids = list(g.nodes())
    node_to_idx = {n: i for i, n in enumerate(node_ids)}
    req_indices = set()
    for u, v, data in g.edges(data=True):
        if int(data.get("required", 0)) == 1:
            req_indices.add(node_to_idx[u])
            req_indices.add(node_to_idx[v])

    # Geometry
    hull = ConvexHull(points)
    min_x, max_x = points[:, 0].min(), points[:, 0].max()
    min_y, max_y = points[:, 1].min(), points[:, 1].max()
    val_width = max_x - min_x
    val_height = max_y - min_y

    # Centers
    c_bbox = bbox_center(points)
    c_mean = mean_center(points)
    c_median = median_center(points)

    # Required node geometry
    req_bbox_area_val = calculate_req_bbox_area(points, req_indices)
    req_convex_hull_area_val = calculate_req_convex_hull_area(points, req_indices)

    # Edge centroids
    _, _, _, _, centroids_all, centroids_req = get_edge_centroid_distance_matrices(g)

    # Spatial distances
    dist_depot_bbox = float(np.linalg.norm(depot - c_bbox))
    dist_depot_node_centroid = float(np.linalg.norm(depot - c_mean))

    req_pts = points[list(req_indices)] if req_indices else np.empty((0, 2))
    avg_dist_bbox_to_req = float(np.linalg.norm(req_pts - c_bbox, axis=1).mean()) if len(req_pts) > 0 else 0.0
    avg_dist_centroid_to_req = float(np.linalg.norm(req_pts - c_mean, axis=1).mean()) if len(req_pts) > 0 else 0.0

    if len(centroids_req) >= 2:
        diff_c = centroids_req[:, np.newaxis, :] - centroids_req[np.newaxis, :, :]
        d_c = np.sqrt(np.sum(diff_c ** 2, axis=2))
        avg_pairwise_dist_req_edges_val = float(d_c[np.triu_indices(len(centroids_req), k=1)].mean())
    else:
        avg_pairwise_dist_req_edges_val = 0.0
    avg_dist_depot_req_edges_val = calculate_avg_dist_point_to_centroids(depot, centroids_req)
    avg_dist_req_edges_bbox_center_val = calculate_avg_dist_point_to_centroids(c_bbox, centroids_req)

    # Radius fractions
    nodes_in_50p, nodes_in_75p = calculate_nodes_within_radius_fractions(points, depot, [0.5, 0.75])
    nodes_in_50p_bbox, nodes_in_75p_bbox = calculate_nodes_within_radius_fractions(points, c_bbox, [0.5, 0.75])

    # Required node spatial stats
    avg_pairwise_req, avg_depot_req, req_in_50p, req_in_75p = calculate_req_spatial_stats(
        points, req_indices, depot
    )

    # Bearing / spread
    std_bearing_depot_val = calculate_std_bearings(points, depot)
    std_bearing_node_centroid_val = calculate_std_bearings(points, c_mean)
    std_bearing_bbox_centroid_val = calculate_std_bearings(points, c_bbox)
    avg_std_xy_val, prod_std_xy_val = calculate_xy_spread(points)
    std_dist_depot_val = calculate_std_distances(points, depot)
    std_dist_node_centroid_val = calculate_std_distances(points, c_mean)
    std_dist_bbox_centroid_val = calculate_std_distances(points, c_bbox)

    # Distance matrix
    D = build_distance_matrix(g, node_list)
    circuity_avg = compute_circuity(D, points)
    prop_dead_ends = calculate_prop_dead_ends(g)

    dist_min, dist_max, dist_mean_val, dist_total_median, dist_total_std = calculate_distance_matrix_stats(
        D, node_list
    )

    req_indices_sorted = sorted(req_indices)
    dist_req_min, dist_req_max, dist_req_mean, dist_req_median, dist_req_std = calculate_req_distance_stats(
        D, req_indices_sorted
    )

    min_dist_depot_to_req = calculate_min_depot_to_req(D, req_indices_sorted)

    # Odd nodes / MST
    req_edge_list = [(u, v) for u, v, d in g.edges(data=True) if int(d.get('required', 0)) == 1]
    if req_edge_list:
        sub_req = g.edge_subgraph(req_edge_list)
        odd_nodes = [n for n, deg in sub_req.degree() if deg % 2 != 0]
        odd_indices = [node_to_idx[n] for n in odd_nodes]
    else:
        odd_indices = []
    mst_odd = compute_mst_odd_weight(D, odd_indices)
    std_dist_odd_nodes_val = calculate_std_dist_odd_nodes(D, odd_indices)

    num_req_edges = len(req_edge_list)

    # Sammon
    _, sammon_error_val = sammon_mapping(D)
    sammon_layout_error_val = calculate_sammon_error(g, points)

    # Grid metrics
    dist_hottest_10 = calculate_dist_depot_to_hottest_cell(g, points, depot_pos, grid_size=10)
    dist_hottest_15 = calculate_dist_depot_to_hottest_cell(g, points, depot_pos, grid_size=15)
    avg_dist_active_10 = calculate_avg_dist_depot_to_active_cells(g, points, depot_pos, grid_size=10)
    avg_dist_active_15 = calculate_avg_dist_depot_to_active_cells(g, points, depot_pos, grid_size=15)
    avg_internal_10 = calculate_avg_internal_dist_cells(points, grid_size=10)
    avg_internal_15 = calculate_avg_internal_dist_cells(points, grid_size=15)
    avg_dist_centroids_10 = calculate_avg_dist_between_centroids(points, grid_size=10)
    avg_dist_centroids_15 = calculate_avg_dist_between_centroids(points, grid_size=15)

    # Density
    node_density_val = calculate_node_density(points)
    edge_density_val = calculate_edge_density(g, points)

    # Degree stats
    num_req_val, num_even_req_val, num_odd_req_val = calculate_req_node_degrees(g)
    prop_odd_req_val = calculate_prop_odd_req_nodes(g)

    # Required edge lengths
    req_mean_val, req_median_val, req_std_val = calculate_req_edge_length_stats(g)

    # Edge centroid stats
    avg_all_centroids, avg_req_centroids, std_req_centroids = calculate_edge_centroid_stats(g)

    return Metrics(
        num_nodes=g.number_of_nodes(),
        sqrt_num_nodes=float(np.sqrt(g.number_of_nodes())),
        num_edges=g.number_of_edges(),
        vertices_param=vertices_param,
        seed=str(seed_value),
        required_ratio=required_ratio,
        convex_hull_area=hull.volume,
        bbox_area=float(val_width * val_height),
        req_bbox_area=req_bbox_area_val,
        req_convex_hull_area=req_convex_hull_area_val,
        bbox_perimeter=float(2 * (val_width + val_height)),
        convex_hull_perimeter=hull.area,
        width=val_width,
        height=val_height,
        avg_pairwise_dist_req=avg_pairwise_req,
        avg_dist_depot_req=avg_depot_req,
        avg_pairwise_dist_req_edges=avg_pairwise_dist_req_edges_val,
        avg_dist_depot_req_edges=avg_dist_depot_req_edges_val,
        avg_dist_req_edges_bbox_center=avg_dist_req_edges_bbox_center_val,
        dist_depot_bbox_center=dist_depot_bbox,
        dist_depot_to_node_centroid=dist_depot_node_centroid,
        avg_dist_bbox_center_to_req_nodes=avg_dist_bbox_to_req,
        avg_dist_node_centroid_to_req_nodes=avg_dist_centroid_to_req,
        nodes_within_50p_radius=nodes_in_50p,
        nodes_within_75p_radius=nodes_in_75p,
        req_nodes_within_50p_radius=req_in_50p,
        req_nodes_within_75p_radius=req_in_75p,
        nodes_within_50p_bbox=nodes_in_50p_bbox,
        nodes_within_75p_bbox=nodes_in_75p_bbox,
        avg_dist_mean=avg_distance_to_center(points, c_mean),
        avg_dist_median=avg_distance_to_center(points, c_median),
        avg_dist_bbox=avg_distance_to_center(points, c_bbox),
        std_bearing_depot=std_bearing_depot_val,
        std_bearing_node_centroid=std_bearing_node_centroid_val,
        std_bearing_bbox_centroid=std_bearing_bbox_centroid_val,
        avg_std_xy=avg_std_xy_val,
        prod_std_xy=prod_std_xy_val,
        std_dist_depot=std_dist_depot_val,
        std_dist_node_centroid=std_dist_node_centroid_val,
        std_dist_bbox_centroid=std_dist_bbox_centroid_val,
        dist_depot_to_hottest_cell_10=dist_hottest_10,
        dist_depot_to_hottest_cell_15=dist_hottest_15,
        avg_dist_depot_to_active_cells_10=avg_dist_active_10,
        avg_dist_depot_to_active_cells_15=avg_dist_active_15,
        avg_internal_dist_cells_10=avg_internal_10,
        avg_internal_dist_cells_15=avg_internal_15,
        avg_dist_between_centroids_10=avg_dist_centroids_10,
        avg_dist_between_centroids_15=avg_dist_centroids_15,
        node_density=node_density_val,
        edge_density=edge_density_val,
        num_req_nodes=num_req_val,
        num_even_req_nodes=num_even_req_val,
        num_odd_req_nodes=num_odd_req_val,
        circuity_avg=circuity_avg,
        prop_dead_ends=prop_dead_ends,
        prop_odd_req_nodes=prop_odd_req_val,
        dist_min=dist_min,
        dist_max=dist_max,
        dist_mean=dist_mean_val,
        dist_total_median=dist_total_median,
        dist_total_std=dist_total_std,
        dist_req_min=dist_req_min,
        dist_req_max=dist_req_max,
        dist_req_mean=dist_req_mean,
        dist_req_median=dist_req_median,
        dist_req_std=dist_req_std,
        min_dist_depot_to_req=min_dist_depot_to_req,
        mst_odd_weight=mst_odd,
        std_dist_odd_nodes=std_dist_odd_nodes_val,
        num_req_edges=num_req_edges,
        sammon_error=sammon_error_val,
        sammon_layout_error=sammon_layout_error_val,
        req_edges_mean=req_mean_val,
        req_edges_median=req_median_val,
        req_edges_std=req_std_val,
        avg_dist_all_edge_centroids=avg_all_centroids,
        avg_dist_req_edge_centroids=avg_req_centroids,
        std_dist_req_edge_centroids=std_req_centroids,
        graphml_path=os.path.basename(graphml_path),
    )


def normalize_graph(g: nx.Graph) -> nx.Graph:
    """Normaliza atributos de aristas para compatibilidad con extract_metrics.

    Los grafos reales (ej. mapas OSM) usan 'demand' en lugar de 'required'.
    Esta función agrega el atributo 'required' a cada arista que no lo tenga,
    derivándolo de 'demand': demand > 0 → required=1, demand == 0 → required=0.
    Si tampoco existe 'demand', se asigna required=0 por defecto.
    """
    for _, _, data in g.edges(data=True):
        if 'required' not in data:
            data['required'] = 1 if int(data.get('demand', 0)) > 0 else 0
    return g


def analyze_existing(input_dir: Path, output_dir: Path, csv_path: str | None) -> None:
    """Analiza instancias GraphML existentes en input_dir sin generar nuevas instancias.

    Antes de llamar a extract_metrics, cada grafo pasa por normalize_graph,
    que traduce el atributo 'demand' al atributo 'required' esperado por extract_metrics
    (demand > 0 → required=1). El grafo normalizado se escribe en output_dir
    para que extract_metrics pueda leerlo desde disco.
    """
    graphml_files = sorted(input_dir.glob('*.graphml'))
    if not graphml_files:
        print(f"No se encontraron archivos .graphml en '{input_dir}/'")
        return

    output_dir.mkdir(parents=True, exist_ok=True)
    results = []
    total = len(graphml_files)

    for i, src_path in enumerate(graphml_files, 1):
        dst_path = output_dir / src_path.name
        print(f"[{i}/{total}] Analizando {src_path.name} ... ", end="", flush=True)

        g = nx.read_graphml(str(src_path))
        normalize_graph(g)
        nx.write_graphml(g, str(dst_path))

        num_edges = g.number_of_edges()
        req_edges = sum(1 for _, _, d in g.edges(data=True) if int(d.get('required', 0)) == 1)
        required_ratio = req_edges / num_edges if num_edges > 0 else 0.0

        metrics = extract_metrics(str(dst_path), g.number_of_nodes(), src_path.stem, required_ratio)
        results.append(metrics)
        print(metrics.to_string())

    header = ";".join(f.name for f in fields(Metrics))
    print(f"\n{header}")
    for m in results:
        print(m.to_string())

    if csv_path:
        field_names = [f.name for f in fields(Metrics)]
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=field_names, delimiter=";")
            writer.writeheader()
            for m in results:
                writer.writerow({fn: getattr(m, fn) for fn in field_names})
        print(f"\nMétricas exportadas a '{csv_path}'")


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for graph generation and analysis."""
    parser = argparse.ArgumentParser(description="Generate multiple graph instances and extract geometric metrics.")
    parser.add_argument("--mode", choices=["generate", "analyze"], default="generate",
                        help="'generate' crea nuevas instancias; 'analyze' analiza GraphMLs existentes.")
    parser.add_argument("--input-dir", default="mapas", help="Directorio con GraphMLs existentes (modo analyze).")
    parser.add_argument("--vertices", nargs="+", type=int, help="List of vertex counts (modo generate).")
    parser.add_argument("--seeds", nargs="+", default=["0"], help="List of random seeds (modo generate).")
    parser.add_argument("--required", nargs="+", type=float, help="List of required-edge percentages 0-100 (modo generate).")
    parser.add_argument("--generator", type=int, choices=[1, 2], default=1,
                        help="1 = proximity/planar (generateInstance), 2 = Delaunay (generateInstance2).")
    parser.add_argument("--output-dir", default="instances", help="Directory for generated/copied GraphML files.")
    parser.add_argument("--csv", default=None, help="Output CSV path for metrics.")
    return parser.parse_args()


def main() -> None:
    """Main function to generate graph instances and extract metrics."""
    args = parse_args()
    output_dir = Path(args.output_dir)

    if args.mode == "analyze":
        analyze_existing(Path(args.input_dir), output_dir, args.csv)
        return

    if not args.vertices or not args.required:
        print("Error: --vertices y --required son obligatorios en modo generate.")
        return

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
