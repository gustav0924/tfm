import argparse
import csv
import os
from dataclasses import dataclass, fields
from itertools import product
from pathlib import Path
from random import seed as set_seed
import networkx as nx
import numpy as np
import pandas as pd
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
    calculate_avg_internal_dist_cells,
    calculate_dist_depot_to_hottest_cell,
    calculate_distance_matrix_stats,
    calculate_min_depot_to_req,
    calculate_node_density,
    calculate_nodes_within_radius_fractions,
    calculate_prop_dead_ends,
    calculate_req_distance_stats,
    calculate_req_edge_length_stats,
    calculate_req_node_degrees,
    calculate_edge_centroid_stats,
    calculate_req_spatial_stats,
    calculate_sammon_error,
    compute_circuity,
    compute_mst_odd_weight,
    mean_center,
    median_center,
    sammon_mapping,
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
    num_edges: int
    vertices_param: int
    seed: str
    required_ratio: float
    convex_hull_area: float
    bbox_area: float
    req_bbox_area: float
    req_convex_hull_area: float
    bbox_perimeter: float
    convex_hull_perimeter: float
    width: float
    height: float
    avg_pairwise_dist_req: float
    avg_dist_depot_req: float
    avg_pairwise_dist_req_edges: float  # Nueva característica: Distancia por pares (aristas)
    avg_dist_depot_req_edges: float     # Nueva característica: Distancia depot a aristas
    avg_dist_req_edges_bbox_center: float # Nueva característica: Distancia aristas a centro BBox
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
    dist_total_median: float
    dist_total_std: float
    dist_req_min: float
    dist_req_max: float
    dist_req_mean: float
    dist_req_median: float
    dist_req_std: float
    min_dist_depot_to_req: float
    mst_odd_weight: float
    num_req_edges: int
    sammon_error: float
    sammon_layout_error: float
    req_edges_mean: float
    req_edges_median: float
    req_edges_std: float
    avg_dist_all_edge_centroids: float
    avg_dist_req_edge_centroids: float
    std_dist_req_edge_centroids: float
    graphml_path: str

    def to_string(self) -> str:
        values = [str(getattr(self, f.name)) for f in fields(self)]
        return ";".join(values)


def generate_instance_1(num_vertices, required_ratio, seed_value) -> nx.Graph:
    set_seed(seed_value)
    while True:
        vertices, edges = generateGraph(num_vertices, required_ratio)
        if testGraph(vertices, edges) != "No":
            candidate = to_networkx(vertices, edges)
            if has_single_strong_component(candidate):
                return candidate

def generate_instance_2(num_vertices, required_ratio, seed_value) -> nx.Graph:
    set_seed(seed_value)
    while True:
        vertices = generate_vertices(num_vertices)
        edges = build_voronoi_adjacency_edges(vertices)
        assign_required(edges, required_ratio)
        g = to_networkx2(vertices, edges)
        if has_single_strong_component2(g):
            return g

def get_edge_centroid_distance_matrices(g):
    pos = {n: (float(d.get('x', 0)), float(d.get('y', 0))) for n, d in g.nodes(data=True)}
    
    all_centroids = []
    req_centroids = []
    
    all_edges_list = []
    req_edges_list = []

    for u, v, d in g.edges(data=True):
        mx = (pos[u][0] + pos[v][0]) / 2.0
        my = (pos[u][1] + pos[v][1]) / 2.0
        centroid = [mx, my]
        
        all_centroids.append(centroid)
        all_edges_list.append((u, v)) 
        
        if int(d.get('required', d.get('d3', 0))) == 1:
            req_centroids.append(centroid)
            req_edges_list.append((u, v))

    def calculate_distance_matrix(centroids_list):
        A = np.array(centroids_list)
        if len(A) < 2:
            return np.zeros((len(A), len(A)))
        diff = A[:, np.newaxis, :] - A[np.newaxis, :, :]
        return np.sqrt(np.sum(diff ** 2, axis=2))

    dist_matrix_all = calculate_distance_matrix(all_centroids)
    dist_matrix_req = calculate_distance_matrix(req_centroids)

    return (
        dist_matrix_all, 
        dist_matrix_req, 
        all_edges_list, 
        req_edges_list, 
        np.array(all_centroids), 
        np.array(req_centroids) 
    )

def mean_center(points: np.ndarray) -> np.ndarray:
    return points.mean(axis=0)

def median_center(points: np.ndarray) -> np.ndarray:
    return np.median(points, axis=0)

def bbox_center(points: np.ndarray) -> np.ndarray:
    return (points.min(axis=0) + points.max(axis=0)) / 2.0

def avg_distance_to_center(points: np.ndarray, center: np.ndarray) -> float:
    diffs = points - center
    distances = np.sqrt((diffs ** 2).sum(axis=1))
    return float(distances.mean())

def calculate_dist_depot_to_hottest_cell(g, points, depot_pos, grid_size):
    pos = {n: (float(d['x']), float(d['y'])) for n, d in g.nodes(data=True)}
    req_edges = [e for e in g.edges(data=True) if int(e[2].get('required', e[2].get('d3', 0))) == 1]
    
    if not req_edges:
        return 0.0

    min_x, max_x = points[:, 0].min(), points[:, 0].max()
    min_y, max_y = points[:, 1].min(), points[:, 1].max()
    rx, ry = (max_x - min_x) or 1.0, (max_y - min_y) or 1.0

    grid_counts = np.zeros((grid_size, grid_size))
    
    for u, v, _ in req_edges:
        mx, my = (pos[u][0] + pos[v][0]) / 2, (pos[u][1] + pos[v][1]) / 2
        col = min(int((mx - min_x) / rx * grid_size), grid_size - 1)
        row = min(int((my - min_y) / ry * grid_size), grid_size - 1)
        grid_counts[row, col] += 1

    max_idx = np.unravel_index(np.argmax(grid_counts), grid_counts.shape)
    hottest_center_x = min_x + (max_idx[1] + 0.5) * (rx / grid_size)
    hottest_center_y = min_y + (max_idx[0] + 0.5) * (ry / grid_size)
    
    return float(np.linalg.norm(depot_pos - np.array([hottest_center_x, hottest_center_y])))

def calculate_avg_dist_depot_to_active_cells(g, points, depot_pos, grid_size):
    pos = {n: (float(d['x']), float(d['y'])) for n, d in g.nodes(data=True)}
    req_edges = [e for e in g.edges(data=True) if int(e[2].get('required', e[2].get('d3', 0))) == 1]
    
    if not req_edges:
        return 0.0

    min_x, max_x = points[:, 0].min(), points[:, 0].max()
    min_y, max_y = points[:, 1].min(), points[:, 1].max()
    rx, ry = (max_x - min_x) or 1.0, (max_y - min_y) or 1.0

    active_cells = set()
    for u, v, _ in req_edges:
        mx, my = (pos[u][0] + pos[v][0]) / 2, (pos[u][1] + pos[v][1]) / 2
        col = min(int((mx - min_x) / rx * grid_size), grid_size - 1)
        row = min(int((my - min_y) / ry * grid_size), grid_size - 1)
        active_cells.add((row, col))

    distances = []
    for row, col in active_cells:
        cx = min_x + (col + 0.5) * (rx / grid_size)
        cy = min_y + (row + 0.5) * (ry / grid_size)
        dist = np.linalg.norm(depot_pos - np.array([cx, cy]))
        distances.append(dist)

    return float(np.mean(distances)) if distances else 0.0

def calculate_avg_dist_between_centroids(points, grid_size):
    if len(points) < 2:
        return 0.0

    min_x, max_x = points[:, 0].min(), points[:, 0].max()
    min_y, max_y = points[:, 1].min(), points[:, 1].max()
    rx, ry = (max_x - min_x) or 1.0, (max_y - min_y) or 1.0
    norm_x = (points[:, 0] - min_x) / rx
    norm_y = (points[:, 1] - min_y) / ry
    norm_points = np.column_stack((norm_x, norm_y))
    grid = np.floor(np.clip(norm_points, 0, 0.9999) * grid_size).astype(int)
    active_cells = np.unique(grid, axis=0)

    if len(active_cells) < 2:
        return 0.0
    centers_norm = (active_cells + 0.5) / float(grid_size)
    centers_real_x = min_x + centers_norm[:, 0] * rx
    centers_real_y = min_y + centers_norm[:, 1] * ry
    centers_real = np.column_stack((centers_real_x, centers_real_y))

    return float(pdist(centers_real).mean())

def calculate_avg_internal_dist_cells(points, grid_size):
    if len(points) < 2:
        return 0.0

    min_x, max_x = points[:, 0].min(), points[:, 0].max()
    min_y, max_y = points[:, 1].min(), points[:, 1].max()
    rx, ry = (max_x - min_x) or 1.0, (max_y - min_y) or 1.0

    df = pd.DataFrame(points, columns=['x', 'y'])
    
    df['cell_x'] = np.clip(((df['x'] - min_x) / rx * grid_size).astype(int), 0, grid_size - 1)
    df['cell_y'] = np.clip(((df['y'] - min_y) / ry * grid_size).astype(int), 0, grid_size - 1)
    df['cell'] = list(zip(df['cell_x'], df['cell_y']))

    internal_means = df.groupby('cell').apply(
        lambda g: pdist(g[['x', 'y']]).mean() if len(g) > 1 else 0.0,
        include_groups=False
    )

    return float(internal_means.mean()) if not internal_means.empty else 0.0

def calculate_node_density(points):
    n_total = len(points)
    if n_total == 0:
        return 0.0

    min_x, max_x = points[:, 0].min(), points[:, 0].max()
    min_y, max_y = points[:, 1].min(), points[:, 1].max()
    
    area = (max_x - min_x) * (max_y - min_y)

    return float(n_total / area) if area > 0 else 0.0

def calculate_req_node_degrees(g):
    req_edges = [(u, v, d) for u, v, d in g.edges(data=True) if int(d.get('required', d.get('d3', 0))) == 1]
    
    if not req_edges:
        return 0, 0, 0

    sub_req = nx.Graph()
    sub_req.add_edges_from(req_edges)

    total_count = sub_req.number_of_nodes()

    even_count = sum(1 for n, d in sub_req.degree() if d % 2 == 0)
    odd_count = sum(1 for n, d in sub_req.degree() if d % 2 != 0)

    return total_count, even_count, odd_count

def calculate_sammon_error(g, points):
    n = len(points)
    if n < 2:
        return 0.0

    try:
        D_matrix = nx.floyd_warshall_numpy(g, nodelist=range(n), weight='length')
        D = np.asarray(D_matrix) 
    except Exception:
        return 0.0

    diff = points[:, np.newaxis, :] - points[np.newaxis, :, :]
    d = np.sqrt(np.sum(diff ** 2, axis=2))

    mask = D > 0
    E = np.sum(np.where(mask, (D - d) ** 2 / np.where(mask, D, 1.0), 0.0))
    den = np.sum(D) / 2.0

    return float(E / den) if den > 0 else 0.0

def calculate_req_edge_length_stats(g):
    req_lengths = []
    for u, v, d in g.edges(data=True):
        if int(d.get('required', d.get('d3', 0))) == 1:
            dist = float(d.get('length', d.get('weight', 0.0)))
            req_lengths.append(dist)
    
    if not req_lengths:
        return 0.0, 0.0, 0.0

    mean_val = float(np.mean(req_lengths))
    median_val = float(np.median(req_lengths))
    std_val = float(np.std(req_lengths))

    return mean_val, median_val, std_val

def extract_metrics(graphml_path, vertices_param, seed_value, required_ratio) -> Metrics:
    g = nx.read_graphml(graphml_path)

    positions = []
    for _, data in g.nodes(data=True):
        positions.append((float(data["x"]), float(data["y"])))
    points = np.array(positions)

    node_ids = list(g.nodes())
    node_to_idx = {n: i for i, n in enumerate(node_ids)}
    req_indices = set()
    for u, v, data in g.edges(data=True):
        if int(data.get("required", 0)) == 1:
            req_indices.add(node_to_idx[u])
            req_indices.add(node_to_idx[v])

    nodes_in_50p, nodes_in_75p = calculate_nodes_within_radius_fractions(points, depot, [0.5, 0.75])

    avg_pairwise_req, avg_depot_req, req_in_50p, req_in_75p = calculate_req_spatial_stats(
        points, req_indices, depot
    )

    c_mean = mean_center(points)
    c_median = median_center(points)

    dist_depot_bbox = float(np.linalg.norm(depot - c_bbox))

    nodes_in_50p_bbox, nodes_in_75p_bbox = calculate_nodes_within_radius_fractions(
        points, c_bbox, [0.5, 0.75]
    )

    D = build_distance_matrix(g, node_list)
    circuity_avg = compute_circuity(D, points)

    prop_dead_ends = calculate_prop_dead_ends(g)

    # total graph distance stats (upper triangle of D)
    dist_min, dist_max, dist_mean_val, dist_total_median, dist_total_std = calculate_distance_matrix_stats(
        D, node_list
    )

    # required subset distance stats
    req_indices_sorted = sorted(req_indices)
    dist_req_min, dist_req_max, dist_req_mean, dist_req_median, dist_req_std = calculate_req_distance_stats(
        D, req_indices_sorted
    )

    # depot proximity: minimum network distance from node 0 to any required node
    min_dist_depot_to_req = calculate_min_depot_to_req(D, req_indices_sorted)

    # MST on odd-degree nodes using full shortest-path distances
    req_edge_list = [(u, v) for u, v, d in g.edges(data=True) if int(d.get('required', 0)) == 1]
    if req_edge_list:
        sub_req = g.edge_subgraph(req_edge_list)
        odd_nodes = [n for n, deg in sub_req.degree() if deg % 2 != 0]
        odd_indices = [node_to_idx[n] for n in odd_nodes]
    else:
        odd_indices = []
    mst_odd = compute_mst_odd_weight(D, odd_indices)

    num_req_edges = len(req_edge_list)

    _, sammon_error_val = sammon_mapping(D)
    sammon_layout_error_val = calculate_sammon_error(g, points)

    depot_pos = points[0]
 
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
    req_mean_val, req_median_val, req_std_val = calculate_req_edge_length_stats(g)
    avg_all_centroids, avg_req_centroids, std_req_centroids = calculate_edge_centroid_stats(g)

    return Metrics(
        num_nodes=g.number_of_nodes(),
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
        avg_pairwise_dist_req_edges=avg_pairwise_dist_req_edges_val,   # Asignación
        avg_dist_depot_req_edges=avg_dist_depot_req_edges_val,         # Asignación
        avg_dist_req_edges_bbox_center=avg_dist_req_edges_bbox_center_val, # Asignación
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
        dist_total_median=dist_total_median,
        dist_total_std=dist_total_std,
        dist_req_min=dist_req_min,
        dist_req_max=dist_req_max,
        dist_req_mean=dist_req_mean,
        dist_req_median=dist_req_median,
        dist_req_std=dist_req_std,
        min_dist_depot_to_req=min_dist_depot_to_req,
        mst_odd_weight=mst_odd,
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