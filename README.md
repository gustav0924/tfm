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
| `--mode` | No | `generate` | `generate` = crea instancias nuevas; `analyze` = analiza GraphMLs existentes |
| `--input-dir` | No | `mapas` | Carpeta con GraphMLs existentes (solo modo `analyze`) |
| `--vertices` | Solo `generate` | - | Lista cantidad de vértices (un espacio por valor) |
| `--seeds` | No | `0` | Lista con random seeds (un espacio por valor) |
| `--required` | Solo `generate` | - | Lista con ejes requeridos porcentualmente 0-100 (un espacio por valor) |
| `--generator` | No | `1` | `1` = proximity/planar, `2` = Delaunay |
| `--output-dir` | No | `instances` | Carpeta para generar/copiar archivos GraphML |
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
num_nodes;sqrt_num_nodes;num_edges;...;avg_dist_all_edge_centroids;avg_dist_req_edge_centroids;std_dist_req_edge_centroids;graphml_path
```

#### Archivo CSV (Opcional)

Mismas columnas que output de consola separadas por `;`. Guardado en el directorio especificado con `--csv`.

### Metricas

| Campo | Descripción |
|---|---|
| `num_nodes` | Cantidad de nodos |
| `sqrt_num_nodes` | Raíz cuadrada del número de nodos |
| `num_edges` | Cantidad de aristas |
| `vertices_param` | Cantidad de vértices (parámetro de entrada) |
| `seed` | Semilla aleatoria |
| `required_ratio` | Proporción de aristas requeridas (0–1) |
| `convex_hull_area` | Área del casco convexo de todos los nodos |
| `bbox_area` | Área de la caja envolvente de todos los nodos |
| `req_bbox_area` | Área de la caja envolvente de los nodos requeridos |
| `req_convex_hull_area` | Área del casco convexo de los nodos requeridos |
| `bbox_perimeter` | Perímetro de la caja envolvente |
| `convex_hull_perimeter` | Perímetro del casco convexo |
| `width` | Ancho (extensión en X) de las posiciones de los nodos |
| `height` | Alto (extensión en Y) de las posiciones de los nodos |
| `avg_pairwise_dist_req` | Distancia media entre pares de nodos requeridos |
| `avg_dist_depot_req` | Distancia media depot→nodos requeridos |
| `avg_pairwise_dist_req_edges` | Distancia media entre pares de centroides de aristas requeridas |
| `avg_dist_depot_req_edges` | Distancia media depot→centroides de aristas requeridas |
| `avg_dist_req_edges_bbox_center` | Distancia media centro bbox→centroides de aristas requeridas |
| `dist_depot_bbox_center` | Distancia depot→centro de la caja envolvente |
| `dist_depot_to_node_centroid` | Distancia depot→centroide medio de los nodos |
| `avg_dist_bbox_center_to_req_nodes` | Distancia media centro bbox→nodos requeridos |
| `avg_dist_node_centroid_to_req_nodes` | Distancia media centroide medio→nodos requeridos |
| `nodes_within_50p_radius` | Nodos dentro del 50% del radio máximo al depot |
| `nodes_within_75p_radius` | Nodos dentro del 75% del radio máximo al depot |
| `req_nodes_within_50p_radius` | Nodos requeridos dentro del 50% de su radio desde su centroide |
| `req_nodes_within_75p_radius` | Nodos requeridos dentro del 75% de su radio desde su centroide |
| `nodes_within_50p_bbox` | Nodos dentro del 50% del radio al centro bbox |
| `nodes_within_75p_bbox` | Nodos dentro del 75% del radio al centro bbox |
| `avg_dist_mean` | Distancia media de todos los nodos al centroide medio |
| `avg_dist_median` | Distancia media de todos los nodos al centroide mediana |
| `avg_dist_bbox` | Distancia media de todos los nodos al centro de la caja envolvente |
| `std_bearing_depot` | Desviación estándar de los ángulos de orientación desde el depot |
| `std_bearing_node_centroid` | Desviación estándar de los ángulos de orientación desde el centroide medio |
| `std_bearing_bbox_centroid` | Desviación estándar de los ángulos de orientación desde el centro bbox |
| `avg_std_xy` | Promedio de las desviaciones estándar en X e Y de los nodos |
| `prod_std_xy` | Producto de las desviaciones estándar en X e Y de los nodos |
| `std_dist_depot` | Desviación estándar de distancias de nodos al depot |
| `std_dist_node_centroid` | Desviación estándar de distancias de nodos al centroide medio |
| `std_dist_bbox_centroid` | Desviación estándar de distancias de nodos al centro bbox |
| `dist_depot_to_hottest_cell_10` | Distancia depot→celda más densa (rejilla 10×10) |
| `dist_depot_to_hottest_cell_15` | Distancia depot→celda más densa (rejilla 15×15) |
| `avg_dist_depot_to_active_cells_10` | Distancia media depot→celdas activas (rejilla 10×10) |
| `avg_dist_depot_to_active_cells_15` | Distancia media depot→celdas activas (rejilla 15×15) |
| `avg_internal_dist_cells_10` | Distancia interna media entre nodos por celda (10×10) |
| `avg_internal_dist_cells_15` | Distancia interna media entre nodos por celda (15×15) |
| `avg_dist_between_centroids_10` | Distancia media entre centroides de celdas activas (10×10) |
| `avg_dist_between_centroids_15` | Distancia media entre centroides de celdas activas (15×15) |
| `node_density` | Densidad de nodos por unidad de área (caja envolvente) |
| `edge_density` | Densidad de aristas por unidad de área (caja envolvente) |
| `num_req_nodes` | Nodos que tocan al menos una arista requerida |
| `num_even_req_nodes` | Nodos requeridos con grado par en el subgrafo requerido |
| `num_odd_req_nodes` | Nodos requeridos con grado impar en el subgrafo requerido |
| `circuity_avg` | Circuidad promedio (distancia en red / distancia euclídea) |
| `prop_dead_ends` | Proporción de nodos con grado 1 (calles sin salida) |
| `prop_odd_req_nodes` | Proporción de nodos de grado impar en el subgrafo requerido |
| `dist_min` | Distancia mínima entre pares (caminos más cortos, red completa) |
| `dist_max` | Distancia máxima entre pares (caminos más cortos, red completa) |
| `dist_mean` | Distancia media entre pares (caminos más cortos, red completa) |
| `dist_total_median` | Mediana de distancias entre pares (red completa) |
| `dist_total_std` | Desviación estándar de distancias entre pares (red completa) |
| `dist_req_min` | Distancia mínima entre pares de nodos requeridos |
| `dist_req_max` | Distancia máxima entre pares de nodos requeridos |
| `dist_req_mean` | Distancia media entre pares de nodos requeridos |
| `dist_req_median` | Mediana de distancias entre pares de nodos requeridos |
| `dist_req_std` | Desviación estándar de distancias entre pares de nodos requeridos |
| `min_dist_depot_to_req` | Distancia mínima en red desde el depot hasta cualquier nodo requerido |
| `mst_odd_weight` | Peso del MST sobre nodos de grado impar del subgrafo requerido |
| `std_dist_odd_nodes` | Desviación estándar de distancias entre nodos de grado impar |
| `num_req_edges` | Cantidad de aristas requeridas |
| `sammon_error` | Error de distorsión Sammon (proyección 2D vs distancias en red) |
| `sammon_layout_error` | Error de Sammon del layout original (posiciones 2D del grafo vs distancias en red) |
| `req_edges_mean` | Media de la longitud de aristas requeridas |
| `req_edges_median` | Mediana de la longitud de aristas requeridas |
| `req_edges_std` | Desviación estándar de la longitud de aristas requeridas |
| `avg_dist_all_edge_centroids` | Distancia media entre centroides de todas las aristas |
| `avg_dist_req_edge_centroids` | Distancia media entre centroides de aristas requeridas |
| `std_dist_req_edge_centroids` | Desviación estándar de distancias entre centroides de aristas requeridas |
| `graphml_path` | Nombre del archivo GraphML generado |

