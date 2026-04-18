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
from scipy.spatial.distance import pdist
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

def get_edge_centroid_distance_matrices(g):
    """
    Función de INSUMO. Extrae la micro-geometría espacial de las aristas.
    Devuelve 6 elementos:
    1. dist_matrix_all: Matriz cuadrada de distancias entre todas las aristas.
    2. dist_matrix_req: Matriz cuadrada de distancias entre aristas obligatorias.
    3. all_edges_list: Lista con el orden (u, v) de todas las aristas.
    4. req_edges_list: Lista con el orden (u, v) de las aristas obligatorias.
    5. centroids_all: Matriz (N x 2) con coordenadas [X, Y] de todos los centroides.
    6. centroids_req: Matriz (M x 2) con coordenadas [X, Y] de centroides obligatorios.
    """
    import numpy as np

    # Obtenemos las coordenadas de los nodos
    pos = {n: (float(d.get('x', 0)), float(d.get('y', 0))) for n, d in g.nodes(data=True)}
    
    all_centroids = []
    req_centroids = []
    
    all_edges_list = []
    req_edges_list = []

    # Recorremos el grafo en el orden original
    for u, v, d in g.edges(data=True):
        # Calculamos el punto medio (centroide) de la arista
        mx = (pos[u][0] + pos[v][0]) / 2.0
        my = (pos[u][1] + pos[v][1]) / 2.0
        centroid = [mx, my]
        
        all_centroids.append(centroid)
        all_edges_list.append((u, v)) 
        
        if int(d.get('required', d.get('d3', 0))) == 1:
            req_centroids.append(centroid)
            req_edges_list.append((u, v))

    # Generador de matrices de distancia euclídea
    def calculate_distance_matrix(centroids_list):
        A = np.array(centroids_list)
        if len(A) < 2:
            return np.zeros((len(A), len(A)))
        diff = A[:, np.newaxis, :] - A[np.newaxis, :, :]
        return np.sqrt(np.sum(diff ** 2, axis=2))

    dist_matrix_all = calculate_distance_matrix(all_centroids)
    dist_matrix_req = calculate_distance_matrix(req_centroids)

    return (
        dist_matrix_all, #Matriz con las distancias entre centroides de aristas globales
        dist_matrix_req, #Matriz con las distancias entre centroides de aristas requeridas
        all_edges_list, #Lista para identificación de aristas globales (Nodo A, Nodo B)
        req_edges_list, #Lista para identificación de aristas requeridas (Nodo A, Nodo B)
        np.array(all_centroids), #Coodenadas de cada centroide  de aristas global
        np.array(req_centroids) #Coodenadas de cada centroide  de aristas requerida
    )

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

def calculate_dist_depot_to_hottest_cell(g, points, depot_pos, grid_size):
    """
    Divición del área en una rejilla dinámica y calculo de la distancia del depósito 
    al centro de la celda con más aristas obligatorias.
    """
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
    """
    Calculo de la distancia promedio desde el depósito a los centros de todas 
    las celdas que contienen demanda, usando una resolución de rejilla dinámica.
    """
    import numpy as np

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
        # Multiplicamos por grid_size y limitamos a grid_size - 1
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
    """
    Identifica celdas activas (con nodos) en una rejilla dinámica y calcula 
    la separación promedio entre sus centros en la escala real del grafo.
    """

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
    """
    Calculo de la distancia promedio entre pares de nodos dentro de cada celda,
    adaptándose dinámicamente a cualquier resolución (grid_size).
    """

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
    """
    Cálculo de la densidad global de nodos (cantidad de nodos por unidad de área) 
    basada en la caja envolvente (bounding box) del escenario.
    """
    n_total = len(points)
    if n_total == 0:
        return 0.0

    min_x, max_x = points[:, 0].min(), points[:, 0].max()
    min_y, max_y = points[:, 1].min(), points[:, 1].max()
    
    area = (max_x - min_x) * (max_y - min_y)

    return float(n_total / area) if area > 0 else 0.0

def calculate_req_node_degrees(g):
    """
    Construye el subgrafo de demanda una sola vez y extrae la topología completa:
    total de nodos requeridos, nodos con grado par y nodos con grado impar.
    Devuelve una tupla: (total, pares, impares).
    """

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
    """
    Calcula el Error de Sammon. Compara la matriz de distancias más cortas 
    de la red (D) con la matriz de distancias euclídeas (d) de las coordenadas.
    """
    n = len(points)
    if n < 2:
        return 0.0

    try:
        D_matrix = nx.floyd_warshall_numpy(g, nodelist=range(n), weight='length')
        D = np.asarray(D_matrix) # Convertir a array de NumPy estándar
    except Exception:
        return 0.0

    diff = points[:, np.newaxis, :] - points[np.newaxis, :, :]
    d = np.sqrt(np.sum(diff ** 2, axis=2))

    mask = D > 0
    E = np.sum(np.where(mask, (D - d) ** 2 / np.where(mask, D, 1.0), 0.0))
    den = np.sum(D) / 2.0

    return float(E / den) if den > 0 else 0.0

def calculate_req_edge_length_stats(g):
    """
    Extrae las longitudes únicamente de las aristas con demanda (obligatorias) 
    y calcula su media, mediana y desviación estándar en una sola pasada.
    Devuelve una tupla: (media, mediana, std).
    """
    req_lengths = []
    for u, v, d in g.edges(data=True):
        if int(d.get('required', d.get('d3', 0))) == 1:
            # Extraemos la distancia (usamos 'length' o 'weight' como fallback)
            dist = float(d.get('length', d.get('weight', 0.0)))
            req_lengths.append(dist)
    
    if not req_lengths:
        return 0.0, 0.0, 0.0

    # Calcular estadísticas usando NumPy
    mean_val = float(np.mean(req_lengths))
    median_val = float(np.median(req_lengths))
    std_val = float(np.std(req_lengths))

    return mean_val, median_val, std_val

def extract_metrics(graphml_path, vertices_param, seed_value, required_ratio) -> Metrics:
    """Extract geometric metrics from a generated graph instance."""
    g = nx.read_graphml(graphml_path)

    positions = []
    for _, data in g.nodes(data=True):
        positions.append((float(data["x"]), float(data["y"])))
    points = np.array(positions)

    # Variables de longitud
    centroid_matrix_all,centroid_matrix_req, edge_order_all, edge_order_req, centroids_all, centroids_req = get_edge_centroid_distance_matrices(g)

    hull = ConvexHull(points)

    min_x, max_x = points[:, 0].min(), points[:, 0].max()
    min_y, max_y = points[:, 1].min(), points[:, 1].max()

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
    sammon_val = calculate_sammon_error(g, points)
    req_mean_val, req_median_val, req_std_val = calculate_req_edge_length_stats(g)

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
