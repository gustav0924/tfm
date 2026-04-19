# TFM

Scripts para generar instancias y extracción de atributos.

## Requerimientos

```
networkx
numpy
scipy
matplotlib
gurobi
```

## Archivos

| Archivo | Descripción |
|---|---|
| `generateInstance.py` | Generador de grafos sin borde de cruce |
| `generateInstance2.py` | Generador de grafos con Delaunay |
| `generateAndAnalyze.py` | Generador y extractor |

## Uso

### generateAndAnalyze.py

Genera múltiples instancias de grafos y extrae métricas geométricas.

```bash
python generateAndAnalyze.py \
  --vertices 10 20 50 \
  --seeds 0 1 2 3 4 \
  --required 30 50 70 \
  --generator 1 \
  --output-dir instances/ \
  --csv results.csv
```

### Argumentos

| Argumento | Requerido | Default | Descripción |
|---|---|---|---|
| `--vertices` | Yes | - | Lista cantidad de vértices (un espacio por valor) |
| `--seeds` | No | `0` | Lista con random seeds (un espacio por valor) |
| `--required` | Yes | - | Lista con ejes requeridos porcentualmente 0-100 (un espacio por valor) |
| `--generator` | No | `1` | `1` = proximity/planar, `2` = Delaunay |
| `--output-dir` | No | `instances` | Carpeta para generar archivos GraphML |
| `--csv` | No | - | Directorio de salida para CSV (`;` -separación) |

El script genera todas las combinaciones de los `vertices x seeds x required`, siendo `--vertices 10 20 --seeds 0 1 2 --required 30 50` produce 2 x 3 x 2 = 12 instancias.

### Output

#### Archivos GraphML

Un `.graphml` por instancia, Se guarda en `--output-dir` con patrón:

```
graph_<generator>_v<vertices>_s<seed>_r<required>.graphml
```

Ejemplo: `graph_planar_v20_s3_r50.graphml`

#### Output

Métricas impresas separadas por ";", seguido de una tabla:

Ejemplo:
```
num_nodes;num_edges;vertices_param;seed;required_ratio;convex_hull_area;convex_hull_perimeter;width;height;graphml_path
```

#### Archivo CSV (Opcional)

Mismas columnas que output de consola separadas por `;`. Guardado en el directorio especificado con `--csv`.

### Metricas

| Campo | Descripción |
|---|---|
| `num_nodes` | Cantidad de nodos |
| `sqrt_num_nodes` | Raíz cuadrada del nº de nodos |
| `num_edges` | Cantidad de ejes |
| `vertices_param` | Cantidad de vértices (parámetro de entrada) |
| `seed` | Semilla aleatoria |
| `required_ratio` | Porcentaje de ejes requeridos |
| `convex_hull_area` | Área del casco convexo |
| `bbox_area` | Área de la caja envolvente |
| `bbox_perimeter` | Perímetro de la caja envolvente |
| `convex_hull_perimeter` | Perímetro del casco convexo |
| `width` | Ancho de todas las posiciones de los nodos |
| `height` | Largo de todas las posiciones de los nodos |
| `avg_pairwise_dist_req` | Distancia media entre pares de nodos obligatorios |
| `avg_dist_depot_req` | Distancia media depot→nodos obligatorios |
| `dist_depot_bbox_center` | Distancia depot→centro de la caja envolvente |
| `nodes_within_50p_radius` | Nodos dentro del 50% del radio máximo al depot |
| `nodes_within_75p_radius` | Nodos dentro del 75% del radio máximo al depot |
| `req_nodes_within_50p_radius` | Nodos obligatorios dentro del 50% de su radio |
| `req_nodes_within_75p_radius` | Nodos obligatorios dentro del 75% de su radio |
| `nodes_within_50p_bbox` | Nodos dentro del 50% del radio al centro bbox |
| `nodes_within_75p_bbox` | Nodos dentro del 75% del radio al centro bbox |
| `avg_dist_mean` | Distancia media de todos los nodos al centroide (media) |
| `avg_dist_median` | Distancia media de todos los nodos al centroide (mediana) |
| `avg_dist_bbox` | Distancia media de todos los nodos al centro de la caja envolvente |
| `dist_depot_to_hottest_cell_10` | Distancia depot→celda más densa (rejilla 10×10) |
| `dist_depot_to_hottest_cell_15` | Distancia depot→celda más densa (rejilla 15×15) |
| `avg_dist_depot_to_active_cells_10` | Distancia media depot→celdas activas (rejilla 10×10) |
| `avg_dist_depot_to_active_cells_15` | Distancia media depot→celdas activas (rejilla 15×15) |
| `avg_internal_dist_cells_10` | Distancia interna media entre nodos por celda (10×10) |
| `avg_internal_dist_cells_15` | Distancia interna media entre nodos por celda (15×15) |
| `avg_dist_between_centroids_10` | Distancia media entre centroides de celdas activas (10×10) |
| `avg_dist_between_centroids_15` | Distancia media entre centroides de celdas activas (15×15) |
| `node_density` | Densidad de nodos por unidad de área (caja envolvente) |
| `num_req_nodes` | Nodos que tocan al menos una arista obligatoria |
| `num_even_req_nodes` | Nodos obligatorios con grado par |
| `num_odd_req_nodes` | Nodos obligatorios con grado impar |
| `circuity_avg` | Circuidad promedio (distancia red / distancia euclídea) |
| `prop_dead_ends` | Proporción de nodos con grado 1 |
| `dist_min` | Distancia mínima entre pares (caminos más cortos) |
| `dist_max` | Distancia máxima entre pares (caminos más cortos) |
| `dist_mean` | Distancia media entre pares (caminos más cortos) |
| `mst_odd_weight` | Peso del MST sobre nodos de grado impar del subgrafo obligatorio |
| `num_req_edges` | Cantidad de aristas obligatorias |
| `sammon_path` | Ruta del archivo de coordenadas Sammon (`.npy`) |
| `sammon_error` | Error de distorsión Sammon (red vs. plano 2D) |
| `req_edges_mean` | Media de la longitud de aristas obligatorias |
| `req_edges_median` | Mediana de la longitud de aristas obligatorias |
| `req_edges_std` | Desviación estándar de la longitud de aristas obligatorias |
| `graphml_path` | Nombre del archivo GraphML generado |

