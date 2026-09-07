# src/viz_overlay_centrality.py
# Superpone centralidad (PageRank por defecto) sobre el ÚLTIMO snapshot (ALL).
# - Color base de nodos: tiempo de infección (gris = jamás/todavía).
# - Overlay: TOP_K por centralidad con halo y borde resaltado, tamaño mayor y etiqueta.
# - Exporta además un CSV con el ranking de centralidad.

import pandas as pd, networkx as nx, matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
from pathlib import Path

# ---------- PARAMS ----------
METRIC   = "pagerank"   # 'pagerank' | 'betweenness' | 'out_degree'
TOP_K    = 10           # cuántos nodos resaltar
K_FACTOR = 3.8          # layout base (más alto = más separación)
SCALE    = 8.0
ITER     = 1000
MIN_SEP      = 0.10     # separación mínima post-proceso
SEP_STEPS    = 500
SEP_STEPSIZE = 0.02

ROOT = Path(__file__).resolve().parents[1]
P    = ROOT/"data/processed"
R    = ROOT/"results"
F    = ROOT/"figs"; F.mkdir(exist_ok=True)

# ---------- helpers ----------
def enforce_min_separation(pos, min_sep=0.08, steps=400, step_size=0.02):
    keys = list(pos.keys())
    P = np.array([pos[k] for k in keys], dtype=float)
    P = (P - P.min(0)) / np.maximum(np.ptp(P, axis=0), 1e-9)
    for _ in range(steps):
        moved = False
        for i in range(len(P)):
            delta = P[i] - P
            dist2 = (delta**2).sum(axis=1)
            mask = (dist2 > 0) & (dist2 < (min_sep**2))
            if not np.any(mask): continue
            d = np.sqrt(dist2[mask])[:, None]
            push = (delta[mask] / np.maximum(d,1e-9) * (min_sep - d)).sum(axis=0)
            P[i] += step_size * push
            moved = True
        if not moved: break
        P = (P - P.min(0)) / np.maximum(np.ptp(P, axis=0), 1e-9)
    for k, xy in zip(keys, P):
        pos[k] = xy
    return pos

# ---------- carga ----------
nodes    = pd.read_csv(P/"nodes.csv")
contacts = pd.read_csv(P/"contacts.csv")
inf      = pd.read_csv(R/"infection_edges.csv")
series   = pd.read_csv(R/"states_by_t.csv")

nodes["id"]        = nodes["id"].astype(str)
contacts["source"] = contacts["source"].astype(str)
contacts["target"] = contacts["target"].astype(str)
inf["source"]      = inf["source"].astype(str)
inf["target"]      = inf["target"].astype(str)
inf["t_infect"]    = inf["t_infect"].astype(int)

ids  = nodes["id"].tolist()
TALL = int(inf["t_infect"].max()) if len(inf) else 0

# Colores base: por tiempo de infección
first_inf = inf.groupby("target")["t_infect"].min().to_dict()
if len(first_inf)>0:
    TMIN, TMAX = int(min(first_inf.values())), int(max(first_inf.values()))
else:
    TMIN, TMAX = 0, 1
cmap_time = plt.colormaps["plasma"]
norm_time = colors.Normalize(vmin=TMIN, vmax=TMAX)
def color_by_time(v):
    ti = first_inf.get(v, None)
    return "#d9d9d9" if (ti is None) else cmap_time(norm_time(ti))

# Layout según contactos + separación
G_layout = nx.from_pandas_edgelist(contacts, "source", "target", create_using=nx.Graph())
G_layout.add_nodes_from(ids)
N = max(len(ids), 1)
k = K_FACTOR / np.sqrt(N)
pos = nx.spring_layout(G_layout, seed=7, k=k, iterations=ITER, scale=SCALE)
pos = enforce_min_separation(pos, MIN_SEP, SEP_STEPS, SEP_STEPSIZE)

# Red de contagios completa (hasta TALL)
Ginf = nx.DiGraph(); Ginf.add_nodes_from(ids)
for r in inf.itertuples():
    Ginf.add_edge(r.source, r.target, t_infect=int(r.t_infect))

# ---------- centralidad ----------
cent = {v:0.0 for v in ids}
if METRIC == "pagerank":
    try:
        cent = nx.pagerank(Ginf, alpha=0.85)
    except Exception:
        cent = nx.pagerank_numpy(Ginf, alpha=0.85)
elif METRIC == "betweenness":
    k_sample = min(200, Ginf.number_of_nodes())
    cent = nx.betweenness_centrality(Ginf, k=k_sample, seed=42) if k_sample>0 else cent
elif METRIC == "out_degree":
    for v, d in Ginf.out_degree(): cent[v] = float(d)
else:
    raise SystemExit(f"Métrica desconocida: {METRIC}")

# ranking CSV
pd.DataFrame({"id": ids, METRIC: [cent.get(v,0.0) for v in ids]})\
  .sort_values(METRIC, ascending=False)\
  .to_csv(R/f"rank_{METRIC}.csv", index=False)

# normalización para tamaño del overlay
vals = np.array([cent.get(v,0.0) for v in ids], dtype=float)
vmin, vmax = float(vals.min()), float(vals.max())
norm_vals = np.zeros_like(vals) if vmax-vmin<1e-12 else (vals - vmin)/(vmax - vmin)

# top-K para halo y etiqueta
top_idx   = np.argsort(-norm_vals)[:TOP_K]
top_nodes = [ids[i] for i in top_idx]

# ---------- dibujo base (color por tiempo) ----------
fig, ax = plt.subplots(figsize=(14,12))
nodelist = ids
base_colors = [color_by_time(v) for v in ids]

# 1) fantasma (halo) para top-K — nodos más grandes semi-transparentes
nx.draw_networkx_nodes(
    Ginf, pos, nodelist=top_nodes,
    node_color="#ffcc00", alpha=0.25, node_size=360, linewidths=0, ax=ax
)

# 2) todos los nodos con color base por tiempo
nx.draw_networkx_nodes(
    Ginf, pos, nodelist=nodelist,
    node_color=base_colors, node_size=68, linewidths=0.6, edgecolors="#444", ax=ax
)

# 3) bordes resaltados + tamaño mayor para top-K (según centralidad)
sizes_top = [100 + 320*norm_vals[ids.index(v)] for v in top_nodes]
nx.draw_networkx_nodes(
    Ginf, pos, nodelist=top_nodes,
    node_color="none", edgecolors="#d62728", linewidths=2.2,
    node_size=sizes_top, ax=ax
)

# 4) aristas (ancho por t_infect)
if len(inf) > 0:
    tmin, tmax = int(inf["t_infect"].min()), int(inf["t_infect"].max())
    def w(tt): return 2.0 if tmax==tmin else 0.8 + 2.7*(tt - tmin)/(tmax - tmin)
    widths = [w(d.get("t_infect", tmin)) for _,_,d in Ginf.edges(data=True)]
else:
    widths = []
nx.draw_networkx_edges(
    Ginf, pos, ax=ax,
    arrows=True, arrowstyle="-|>", arrowsize=22,
    width=widths, alpha=0.70, connectionstyle="arc3,rad=0.10"
)

# 5) etiquetas SOLO para top-K
labels = {v: v for v in top_nodes}
nx.draw_networkx_labels(Ginf, pos, labels=labels, font_size=8, font_weight="bold", ax=ax)

# barra de color del tiempo (para base)
sm = plt.cm.ScalarMappable(norm=norm_time, cmap=plt.colormaps["plasma"]); sm.set_array([])
fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.03).set_label("Tiempo de infección (t)")

ax.set_title(f"Último snapshot (ALL) — overlay centralidad: {METRIC} (top {TOP_K})")
ax.axis("off"); fig.tight_layout()
out = F/f"infections_ALL_overlay_{METRIC}.png"
fig.savefig(out, dpi=180); plt.close(fig)
print(f"Listo → {out}")
print(f"Ranking → results/rank_{METRIC}.csv")
