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

### Arguments

| Argumento | Requerido | Default | Descripción |
|---|---|---|---|
| `--vertices` | Yes | - | Lista cantidad de vértices (un espacio por valor) |
| `--seeds` | No | `0` | Lista con random seeds (un espacio por valor) |
| `--required` | Yes | - | Lista con ejes requeridos porcentualmente 0-100 (un espacio por valor) |
| `--generator` | No | `1` | `1` = proximity/planar, `2` = Delaunay |
| `--output-dir` | No | `instances` | Carpeta para generar archivos GraphML |
| `--csv` | No | - | Directorio de salida para CSV (`;` -separación) |

El script genenera todas las combinaciones de los `vertices x seeds x required`, siendo `--vertices 10 20 --seeds 0 1 2 --required 30 50` produce 2 x 3 x 2 = 12 instancias.

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

En proceso.

