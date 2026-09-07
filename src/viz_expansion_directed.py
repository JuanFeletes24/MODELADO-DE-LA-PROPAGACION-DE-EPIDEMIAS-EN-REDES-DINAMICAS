# src/viz_expansion_directed.py
# Grafo de CONTAGIOS dirigido con flechas y PESO=t_infect (robusto a "Node X has no position")

import pandas as pd, networkx as nx, matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
P    = ROOT/"data/processed"
RES  = ROOT/"results"
FIGS = ROOT/"figs"; FIGS.mkdir(exist_ok=True)

# --- Cargar datos y normalizar tipos ---
nodes    = pd.read_csv(P/"nodes.csv")
contacts = pd.read_csv(P/"contacts.csv")
inf      = pd.read_csv(RES/"infection_edges.csv")
series   = pd.read_csv(RES/"states_by_t.csv")
TMAX = int(inf["t_infect"].max()) if len(inf) else 0


# Asegurar que TODO es string (coherente)
nodes["id"]    = nodes["id"].astype(str)
contacts["source"] = contacts["source"].astype(str)
contacts["target"] = contacts["target"].astype(str)
inf["source"]  = inf["source"].astype(str)
inf["target"]  = inf["target"].astype(str)
if "t_infect" in inf.columns:
    inf["t_infect"] = inf["t_infect"].astype(int)

ids = nodes["id"].tolist()

# --- Grafo de CONTACTOS (para layout): garantiza posición para todos los nodos ---
G_layout = nx.from_pandas_edgelist(contacts, "source", "target", create_using=nx.Graph())
G_layout.add_nodes_from(ids)   # por si hay nodos sin contactos en alguna ventana
pos = nx.spring_layout(G_layout, seed=7)   # <-- POS para TODOS los nodos

# --- Grafo maestro de CONTAGIOS con TODOS los nodos (para snapshots) ---
G_all = nx.DiGraph(); G_all.add_nodes_from(ids)
for r in inf.itertuples():
    G_all.add_edge(r.source, r.target, t_infect=int(r.t_infect))

# --- Curva I(t) ---
plt.figure(figsize=(8,4))
plt.plot(series["t"], series["I"])
plt.xlabel("t (ventanas de 20 s)"); plt.ylabel("Casos activos (I)")
plt.title("SIS — evolución de casos")
plt.tight_layout(); plt.savefig(FIGS/"sis_I_curve.png", dpi=160); plt.close()

# --- Función de dibujo (dirigido + peso visible) ---
def draw_upto(t, fname, max_labels=150):
    sub = inf[inf["t_infect"] <= t].copy()

    H = nx.DiGraph(); H.add_nodes_from(ids)
    for r in sub.itertuples():
        H.add_edge(r.source, r.target, t_infect=int(r.t_infect))

    # Nodelist robusto: SOLO los que tienen posición
    nodelist = list(pos.keys())

    # Anchos proporcionales a t_infect (más tarde = más ancho)
    if len(sub) > 0:
        tmin, tmax = int(sub["t_infect"].min()), int(sub["t_infect"].max())
        def w(tt): return 2.0 if tmax==tmin else 0.8 + 2.7*(tt - tmin)/(tmax - tmin)
        widths = [w(d.get("t_infect", tmin)) for _,_,d in H.edges(data=True)]
    else:
        widths = []

    plt.figure(figsize=(9,7))
    nx.draw_networkx_nodes(H, pos, nodelist=nodelist, node_size=60, node_color="#d9d9d9", linewidths=0)
    nx.draw_networkx_edges(
        H, pos,
        arrows=True, arrowstyle="-|>", arrowsize=20,
        width=widths, alpha=0.65, connectionstyle="arc3,rad=0.06"
    )

    # Etiquetas de PESO = t_infect (limita para no saturar)
    if len(sub) <= max_labels:
        labels = {(r.source, r.target): int(r.t_infect) for r in sub.itertuples()}
    else:
        sub2 = sub.sort_values("t_infect").head(max_labels)
        labels = {(r.source, r.target): int(r.t_infect) for r in sub2.itertuples()}

    if labels:
        nx.draw_networkx_edge_labels(
            H, pos, edge_labels=labels, font_size=7, rotate=False,
            bbox=dict(facecolor="white", edgecolor="none", alpha=0.6, pad=1.0),
            label_pos=0.5
        )

    plt.title(f"Red de CONTAGIOS dirigida hasta t={t} (peso = tiempo de infección)")
    plt.axis("off"); plt.tight_layout()
    plt.savefig(FIGS/fname, dpi=160); plt.close()

# Snapshots en tiempos fijos + uno con "todos"
for tt in [40, 100, 200, 400, TMAX]:
    draw_upto(tt, f"infections_upto_t{tt}_directed_weighted.png")

# Figura con TODOS los contagios acumulados
draw_upto(TMAX, "infections_ALL_directed_weighted.png")


print("Listo → figs/: sis_I_curve.png + infections_upto_t*_directed_weighted.png")
