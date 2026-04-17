from scipy.spatial.distance import pdist
import networkx as nx
import numpy as np
import pandas as pd


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


def build_distance_matrix(g: nx.Graph, node_list: list) -> np.ndarray:
    """All-pairs shortest path distance matrix via Dijkstra (weight='length'), ordered by node_list."""
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
    """Average ratio of network distance to Euclidean distance over all node pairs."""
    i_idx, j_idx = np.triu_indices(len(points), k=1)
    euc = np.linalg.norm(points[i_idx] - points[j_idx], axis=1)
    net = D[i_idx, j_idx]
    mask = (euc > 0) & (net > 0)
    return float(np.mean(net[mask] / euc[mask])) if mask.any() else 0.0


def compute_mst_odd_weight(g: nx.Graph) -> float:
    """MST weight on the subgraph induced by odd-degree nodes of the required-edges subgraph."""
    req_edges = [(u, v) for u, v, d in g.edges(data=True) if int(d.get('required', 0)) == 1]
    if not req_edges:
        return 0.0
    sub_req = g.edge_subgraph(req_edges)
    odd_nodes = [n for n, deg in sub_req.degree() if deg % 2 != 0]
    if len(odd_nodes) < 2:
        return 0.0
    odd_subgraph = g.subgraph(odd_nodes)
    if not nx.is_connected(odd_subgraph):
        return 0.0
    mst = nx.minimum_spanning_tree(odd_subgraph, weight='length')
    return float(mst.size(weight='length'))


def sammon_mapping(D: np.ndarray, n_iter: int = 300, lr: float = 0.3) -> np.ndarray:
    """Sammon mapping: returns a (n, 2) array of 2D coordinates that preserve pairwise distances."""
    n = D.shape[0]
    rng = np.random.default_rng(42)
    Y = rng.random((n, 2))
    suma_D = np.sum(D)
    if suma_D == 0:
        return Y
    D_safe = np.where(D == 0, 1e-10, D)
    for _ in range(n_iter):
        diff = Y[:, np.newaxis, :] - Y[np.newaxis, :, :]
        d = np.sqrt(np.sum(diff ** 2, axis=2))
        d_safe = np.where(d == 0, 1e-10, d)
        factor = (D - d) / (D_safe * d_safe)
        grad = -2.0 / suma_D * np.sum(factor[:, :, np.newaxis] * diff, axis=1)
        Y -= lr * grad
    return Y


def calculate_dist_depot_to_hottest_cell(g, points, depot_pos, grid_size):
    """Distance from depot to center of cell with most required edges."""
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
    """Average distance from depot to centers of all cells containing required edges."""
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
    """Average Euclidean distance between centers of active grid cells."""
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
    """Average pairwise distance between nodes within each grid cell."""
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
    """Node count per unit area based on bounding box."""
    n_total = len(points)
    if n_total == 0:
        return 0.0

    min_x, max_x = points[:, 0].min(), points[:, 0].max()
    min_y, max_y = points[:, 1].min(), points[:, 1].max()
    area = (max_x - min_x) * (max_y - min_y)

    return float(n_total / area) if area > 0 else 0.0


def calculate_req_node_degrees(g):
    """Total required nodes, even-degree count, odd-degree count from required-edge subgraph."""
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
    """Sammon error: distortion between network shortest-path distances and 2D Euclidean distances."""
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
    """Mean, median, std of required edge lengths."""
    req_lengths = []
    for u, v, d in g.edges(data=True):
        if int(d.get('required', d.get('d3', 0))) == 1:
            req_lengths.append(float(d.get('length', d.get('weight', 0.0))))

    if not req_lengths:
        return 0.0, 0.0, 0.0

    return float(np.mean(req_lengths)), float(np.median(req_lengths)), float(np.std(req_lengths))
