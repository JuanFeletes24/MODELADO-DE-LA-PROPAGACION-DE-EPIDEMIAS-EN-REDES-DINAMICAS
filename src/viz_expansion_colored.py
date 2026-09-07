# src/viz_expansion_colored.py
# Snapshots de CONTAGIOS con flechas y PESO=t_infect, más separación y color por tiempo/grupo.
# CORREGIDO: los colores de nodos se asignan en el MISMO orden de `ids`.

import pandas as pd, networkx as nx, matplotlib.pyplot as plt
from matplotlib import colors
import numpy as np
from pathlib import Path

# ---------- CONFIG ----------
COLOR_BY      = "time"   # "time" | "group"
NUM_SNAPSHOTS = 10       # cantidad de snapshots automáticos
MAX_LABELS    = 120      # etiquetas de arista (t_infect) como máx.
# Layout (más espacio)
K_FACTOR   = 3.8
SCALE      = 8.0
ITER       = 1000
# Separación mínima extra entre nodos (post-proceso)
MIN_SEP      = 0.10
SEP_STEPS    = 500
SEP_STEPSIZE = 0.02

ROOT = Path(__file__).resolve().parents[1]
P    = ROOT/"data/processed"
R    = ROOT/"results"
F    = ROOT/"figs"; F.mkdir(exist_ok=True)

# ---------- Utilidad: separación mínima (compat. NumPy 2.x) ----------
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

# ---------- Datos ----------
nodes    = pd.read_csv(P/"nodes.csv")
contacts = pd.read_csv(P/"contacts.csv")
inf      = pd.read_csv(R/"infection_edges.csv")
series   = pd.read_csv(R/"states_by_t.csv")

nodes["id"]        = nodes["id"].astype(str)
if "group" not in nodes.columns: nodes["group"] = ""
contacts["source"] = contacts["source"].astype(str)
contacts["target"] = contacts["target"].astype(str)
inf["source"]      = inf["source"].astype(str)
inf["target"]      = inf["target"].astype(str)
inf["t_infect"]    = inf["t_infect"].astype(int)

ids  = nodes["id"].tolist()
TALL = int(inf["t_infect"].max()) if len(inf) else 0

# Primer/último contagio (para espaciar 10 snapshots)
first_inf = inf.groupby("target")["t_infect"].min().to_dict()
if len(first_inf)>0:
    TMIN, TMAX = int(min(first_inf.values())), int(max(first_inf.values()))
else:
    TMIN, TMAX = 1, max(1, TALL)

# Calcula 10 tiempos equiespaciados + ALL
if TMAX > 0:
    snaps = list(np.linspace(TMIN, TMAX, NUM_SNAPSHOTS, dtype=int))
    SNAPSHOTS = sorted(set(snaps + [TMAX]))
else:
    SNAPSHOTS = [0]

# ---------- Layout base (según CONTACTOS) + separación ----------
G_layout = nx.from_pandas_edgelist(contacts, "source", "target", create_using=nx.Graph())
G_layout.add_nodes_from(ids)
N = max(len(ids), 1)
k = K_FACTOR / np.sqrt(N)
pos = nx.spring_layout(G_layout, seed=7, k=k, iterations=ITER, scale=SCALE)
pos = enforce_min_separation(pos, min_sep=MIN_SEP, steps=SEP_STEPS, step_size=SEP_STEPSIZE)

# ---------- Coloreo ----------
cmap_time = plt.colormaps["plasma"]
cmap_tab  = plt.colormaps["tab20"]
norm_time = colors.Normalize(vmin=TMIN, vmax=TMAX)

def color_by_time(v, t_cut):
    ti = first_inf.get(v, None)
    if (ti is None) or (ti > t_cut):  # no infectado aún
        return "#d9d9d9"
    return cmap_time(norm_time(ti))

groups = nodes["group"].astype(str).tolist()
uniq_g = sorted(set([g for g in groups if g != ""]))
g2c = {g: cmap_tab(i/ max(1,len(uniq_g)-1)) for i,g in enumerate(uniq_g)}
def color_by_group(v):
    g = nodes.loc[nodes["id"]==v, "group"].astype(str).values[0]
    return g2c.get(g, "#d9d9d9")

# ---------- Dibujo ----------
def draw_upto(t_cut, fname):
    sub = inf[inf["t_infect"] <= t_cut].copy()
    H = nx.DiGraph(); H.add_nodes_from(ids)
    for r in sub.itertuples():
        H.add_edge(r.source, r.target, t_infect=int(r.t_infect))

    # Colores de nodos (en el MISMO orden de `ids`)
    if COLOR_BY == "time":
        base_colors = [color_by_time(v, t_cut) for v in ids]
    elif COLOR_BY == "group":
        base_colors = [color_by_group(v) for v in ids]
    else:
        base_colors = ["#d9d9d9"]*len(ids)

    # --- PARCHE DE ORDEN ---
    nodelist = ids                              # usamos exactamente este orden
    color_map = dict(zip(ids, base_colors))     # id -> color
    node_colors_ordered = [color_map[v] for v in nodelist]  # colores alineados

    # Anchos de aristas según t_infect
    if len(sub) > 0:
        tmin, tmax = int(sub["t_infect"].min()), int(sub["t_infect"].max())
        def w(tt): return 2.0 if tmax==tmin else 0.8 + 2.7*(tt - tmin)/(tmax - tmin)
        widths = [w(d.get("t_infect", tmin)) for _,_,d in H.edges(data=True)]
    else:
        widths = []

    fig, ax = plt.subplots(figsize=(14,12))

    nx.draw_networkx_nodes(
        H, pos, nodelist=nodelist, node_color=node_colors_ordered,
        node_size=58, linewidths=0.4, edgecolors="#444444", ax=ax
    )
    nx.draw_networkx_edges(
        H, pos, ax=ax,
        arrows=True, arrowstyle="-|>", arrowsize=22,
        width=widths, alpha=0.70, connectionstyle="arc3,rad=0.10"
    )

    # Etiquetas de arista (limitadas)
    if len(sub) <= MAX_LABELS:
        labels = {(r.source, r.target): int(r.t_infect) for r in sub.itertuples()}
    else:
        labels = {(r.source, r.target): int(r.t_infect)
                  for r in sub.sort_values("t_infect").head(MAX_LABELS).itertuples()}
    if labels:
        nx.draw_networkx_edge_labels(
            H, pos, edge_labels=labels, font_size=7, rotate=False,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.6, pad=1.0),
            label_pos=0.5, ax=ax
        )

    if COLOR_BY == "time":
        sm = plt.cm.ScalarMappable(norm=norm_time, cmap=cmap_time); sm.set_array([])
        fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.03).set_label("Tiempo de infección (t)")

    ax.set_title(f"Contagios dirigidos hasta t={t_cut} — color por {COLOR_BY}")
    ax.axis("off"); fig.tight_layout()
    fig.savefig(F/fname, dpi=180); plt.close(fig)

# Curva I(t)
fig, ax = plt.subplots(figsize=(10,5))
ax.plot(series["t"], series["I"])
ax.set_xlabel("t (20 s)"); ax.set_ylabel("I"); ax.set_title("SIS — evolución de casos")
fig.tight_layout(); fig.savefig(F/"sis_I_curve.png", dpi=170); plt.close(fig)

# Snapshots automáticos (NUM_SNAPSHOTS) + ALL
for tt in SNAPSHOTS:
    draw_upto(tt, f"infections_upto_t{tt}_{COLOR_BY}_spaced.png")
draw_upto(TALL, f"infections_ALL_{COLOR_BY}_spaced.png")

print("Listo → figs/: *_{COLOR_BY}_spaced.png (", len(SNAPSHOTS), "snapshots )")
