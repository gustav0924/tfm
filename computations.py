from itertools import combinations

import networkx as nx
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist


def mean_center(points: np.ndarray) -> np.ndarray:
    """Centro como la media de todas las posiciones de los nodos."""
    return points.mean(axis=0)


def median_center(points: np.ndarray) -> np.ndarray:
    """Centro como la mediana de todas las posiciones de los nodos."""
    return np.median(points, axis=0)


def bbox_center(points: np.ndarray) -> np.ndarray:
    """Centro como el punto medio del bounding box."""
    return (points.min(axis=0) + points.max(axis=0)) / 2.0


def avg_distance_to_center(points: np.ndarray, center: np.ndarray) -> float:
    """Distancia euclidiana promedio desde todos los nodos hasta un punto central dado."""
    diffs = points - center
    distances = np.sqrt((diffs ** 2).sum(axis=1))
    return float(distances.mean())


def build_distance_matrix(g: nx.Graph, node_list: list) -> np.ndarray:
    """Matriz de distancias de caminos mínimos entre todos los pares de nodos via Dijkstra (weight='length'), ordenada según node_list."""
    n = len(node_list)
    node_idx = {node: i for i, node in enumerate(node_list)}
    D = np.zeros((n, n))
    for source, lengths in nx.all_pairs_dijkstra_path_length(g, weight='length'):
        i = node_idx[source]
        for target, dist in lengths.items():
            j = node_idx[target]
            D[i, j] = dist
    return D


def compute_circuity(D: np.ndarray, points: np.ndarray) -> float:
    """Ratio promedio entre la distancia en red y la distancia euclidiana sobre todos los pares de nodos."""
    i_idx, j_idx = np.triu_indices(len(points), k=1)
    euc = np.linalg.norm(points[i_idx] - points[j_idx], axis=1)
    net = D[i_idx, j_idx]
    mask = (euc > 0) & (net > 0)
    return float(np.mean(net[mask] / euc[mask])) if mask.any() else 0.0


def compute_mst_odd_weight(D: np.ndarray, odd_indices: list) -> float:
    """Peso del MST sobre el grafo completo de nodos de grado impar usando distancias de caminos mínimos de D."""
    if len(odd_indices) < 2:
        return 0.0
    h = nx.Graph()
    for a, b in combinations(range(len(odd_indices)), 2):
        i, j = odd_indices[a], odd_indices[b]
        h.add_edge(a, b, length=float(D[i, j]))
    mst = nx.minimum_spanning_tree(h, weight='length')
    return float(mst.size(weight='length'))


def sammon_mapping(D: np.ndarray, n_iter: int = 300, lr: float = 0.3) -> tuple[np.ndarray, float]:
    """Mapeo de Sammon: retorna coordenadas 2D y el error de estrés de proyección normalizado.

    El error de estrés es la suma normalizada de diferencias cuadráticas entre D y
    las distancias euclidianas en la proyección 2D — valores menores indican mejor ajuste.
    """
    n = D.shape[0]
    rng = np.random.default_rng(42)
    Y = rng.random((n, 2))
    suma_D = np.sum(D)
    if suma_D == 0:
        return Y, 0.0
    D_safe = np.where(D == 0, 1e-10, D)
    for _ in range(n_iter):
        diff = Y[:, np.newaxis, :] - Y[np.newaxis, :, :]
        d = np.sqrt(np.sum(diff ** 2, axis=2))
        d_safe = np.where(d == 0, 1e-10, d)
        factor = (D - d) / (D_safe * d_safe)
        grad = -2.0 / suma_D * np.sum(factor[:, :, np.newaxis] * diff, axis=1)
        Y -= lr * grad
    # estrés normalizado: sum((D - d)^2) / sum(D^2) sobre triángulo superior
    i_idx, j_idx = np.triu_indices(n, k=1)
    d_final = np.sqrt(np.sum((Y[i_idx] - Y[j_idx]) ** 2, axis=1))
    D_upper = D[i_idx, j_idx]
    denom = np.sum(D_upper ** 2)
    error = float(np.sum((D_upper - d_final) ** 2) / denom) if denom > 0 else 0.0
    return Y, error


def calculate_dist_depot_to_hottest_cell(g, points, depot_pos, grid_size):
    """Distancia desde el depósito hasta el centro de la celda con más aristas requeridas."""
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
    """Distancia promedio desde el depósito hasta los centros de todas las celdas que contienen aristas requeridas."""
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
        distances.append(np.linalg.norm(depot_pos - np.array([cx, cy])))

    return float(np.mean(distances)) if distances else 0.0


def calculate_avg_dist_between_centroids(points, grid_size):
    """Distancia euclidiana promedio entre centros de celdas activas de la grilla."""
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
    """Distancia promedio entre pares de nodos dentro de cada celda de la grilla."""
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
    """Número de nodos por unidad de área basado en el bounding box."""
    n_total = len(points)
    if n_total == 0:
        return 0.0

    min_x, max_x = points[:, 0].min(), points[:, 0].max()
    min_y, max_y = points[:, 1].min(), points[:, 1].max()
    area = (max_x - min_x) * (max_y - min_y)

    return float(n_total / area) if area > 0 else 0.0


def calculate_req_node_degrees(g):
    """Total de nodos requeridos, conteo de grado par e impar en el subgrafo de aristas requeridas."""
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
    """Error de Sammon: distorsión entre las distancias de caminos mínimos en red y las distancias euclidianas 2D."""
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
    """Media, mediana y desviación estándar de las longitudes de aristas requeridas."""
    req_lengths = []
    for u, v, d in g.edges(data=True):
        if int(d.get('required', d.get('d3', 0))) == 1:
            req_lengths.append(float(d.get('length', d.get('weight', 0.0))))

    if not req_lengths:
        return 0.0, 0.0, 0.0

    return float(np.mean(req_lengths)), float(np.median(req_lengths)), float(np.std(req_lengths))


def calculate_distance_matrix_stats(D: np.ndarray, node_list: list) -> tuple[float, float, float, float, float]:
    """Mínimo, máximo, media, mediana y desviación estándar de las entradas del triángulo superior de una matriz de distancias."""
    upper = D[np.triu_indices(len(node_list), k=1)]
    if len(upper) == 0:
        return 0.0, 0.0, 0.0, 0.0, 0.0
    return (
        float(upper.min()),
        float(upper.max()),
        float(upper.mean()),
        float(np.median(upper)),
        float(upper.std()),
    )


def calculate_req_distance_stats(D: np.ndarray, req_indices_sorted: list) -> tuple[float, float, float, float, float]:
    """Mínimo, máximo, media, mediana y desviación estándar de las distancias de caminos mínimos entre nodos requeridos."""
    if len(req_indices_sorted) < 2:
        return 0.0, 0.0, 0.0, 0.0, 0.0
    D_req = D[np.ix_(req_indices_sorted, req_indices_sorted)]
    upper_req = D_req[np.triu_indices(len(req_indices_sorted), k=1)]
    return (
        float(upper_req.min()),
        float(upper_req.max()),
        float(upper_req.mean()),
        float(np.median(upper_req)),
        float(upper_req.std()),
    )


def calculate_min_depot_to_req(D: np.ndarray, req_indices_sorted: list, depot_idx: int = 0) -> float:
    """Distancia mínima en red desde el nodo depósito hasta cualquier nodo requerido."""
    if not req_indices_sorted:
        return 0.0
    return float(D[depot_idx, req_indices_sorted].min())


def calculate_nodes_within_radius_fractions(
    points: np.ndarray, center: np.ndarray, fracs: list[float]
) -> list[int]:
    """Conteo de nodos dentro de cada fracción de la distancia máxima desde el centro.

    Retorna una lista de conteos, uno por fracción en fracs.
    Excluye el nodo central mismo (omite distancia exactamente cero).
    """
    dists = np.linalg.norm(points - center, axis=1)
    # excluir distancia exactamente cero (el nodo central mismo)
    non_zero = dists[dists > 0]
    if len(non_zero) == 0:
        return [0] * len(fracs)
    max_r = float(non_zero.max())
    return [int(np.sum(non_zero <= f * max_r)) for f in fracs]


def calculate_prop_dead_ends(g: nx.Graph) -> float:
    """Proporción de nodos con grado 1 (calles sin salida) en el grafo."""
    n = g.number_of_nodes()
    if n == 0:
        return 0.0
    return len([node for node, deg in g.degree() if deg == 1]) / n


def calculate_req_spatial_stats(
    points: np.ndarray, req_indices: set, depot: np.ndarray
) -> tuple[float, float, int, int]:
    """Estadísticas espaciales de nodos requeridos: distancia promedio entre pares, distancia promedio al depósito, nodos en radio 50%/75%.

    Retorna (avg_pairwise_req, avg_depot_req, req_in_50p, req_in_75p).
    """
    if not req_indices:
        return 0.0, 0.0, 0, 0
    req_points = points[list(req_indices)]
    avg_pairwise = float(pdist(req_points).mean()) if len(req_indices) > 1 else 0.0
    avg_depot = float(np.linalg.norm(req_points - depot, axis=1).mean())
    c_req = req_points.mean(axis=0)
    dist_to_c = np.linalg.norm(req_points - c_req, axis=1)
    non_zero = dist_to_c[dist_to_c > 0]
    if len(non_zero) == 0:
        return avg_pairwise, avg_depot, 0, 0
    max_r = float(non_zero.max())
    in_50p = int(np.sum(non_zero <= 0.5 * max_r))
    in_75p = int(np.sum(non_zero <= 0.75 * max_r))
    return avg_pairwise, avg_depot, in_50p, in_75p


def get_edge_centroid_distance_matrices(g: nx.Graph) -> tuple:
    """Extrae la micro-geometría espacial de las aristas calculando matrices de distancia entre centroides.

    Devuelve 6 elementos:
    - dist_matrix_all: matriz (N×N) de distancias euclidianas entre centroides de todas las aristas.
    - dist_matrix_req: matriz (M×M) de distancias euclidianas entre centroides de aristas requeridas.
    - all_edges_list: lista de tuplas (u, v) con el orden de todas las aristas.
    - req_edges_list: lista de tuplas (u, v) con el orden de las aristas requeridas.
    - centroids_all: array (N×2) con coordenadas [X, Y] de todos los centroides.
    - centroids_req: array (M×2) con coordenadas [X, Y] de centroides de aristas requeridas.
    """
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

    def _dist_matrix(centroids_list):
        A = np.array(centroids_list)
        if len(A) < 2:
            return np.zeros((len(A), len(A)))
        diff = A[:, np.newaxis, :] - A[np.newaxis, :, :]
        return np.sqrt(np.sum(diff ** 2, axis=2))

    return (
        _dist_matrix(all_centroids),
        _dist_matrix(req_centroids),
        all_edges_list,
        req_edges_list,
        np.array(all_centroids) if all_centroids else np.empty((0, 2)),
        np.array(req_centroids) if req_centroids else np.empty((0, 2)),
    )


def calculate_edge_centroid_stats(g: nx.Graph) -> tuple[float, float, float]:
    """Estadísticas de distancia entre centroides de aristas: promedio global, promedio requeridas, std requeridas.

    Retorna (avg_dist_all_centroids, avg_dist_req_centroids, std_dist_req_centroids).
    """
    _, dist_req, _, _, centroids_all, centroids_req = get_edge_centroid_distance_matrices(g)

    if len(centroids_all) >= 2:
        diff = centroids_all[:, np.newaxis, :] - centroids_all[np.newaxis, :, :]
        d_all = np.sqrt(np.sum(diff ** 2, axis=2))
        upper_all = d_all[np.triu_indices(len(centroids_all), k=1)]
        avg_all = float(upper_all.mean())
    else:
        avg_all = 0.0

    if len(centroids_req) >= 2:
        upper_req = dist_req[np.triu_indices(len(centroids_req), k=1)]
        avg_req = float(upper_req.mean())
        std_req = float(upper_req.std())
    else:
        avg_req = 0.0
        std_req = 0.0

    return avg_all, avg_req, std_req
