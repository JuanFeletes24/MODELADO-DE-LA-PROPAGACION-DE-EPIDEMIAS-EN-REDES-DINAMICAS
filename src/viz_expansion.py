import pandas as pd, networkx as nx, matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
P    = ROOT/"data/processed"
RES  = ROOT/"results"
FIGS = ROOT/"figs"; FIGS.mkdir(exist_ok=True)

series = pd.read_csv(RES/"states_by_t.csv")
inf    = pd.read_csv(RES/"infection_edges.csv")
nodes  = pd.read_csv(P/"nodes.csv")

# Curva I(t)
plt.figure(); plt.plot(series["t"], series["I"])
plt.xlabel("t (ventanas de 20 s)"); plt.ylabel("Casos activos (I)")
plt.title("SIS — evolución de casos"); plt.tight_layout()
plt.savefig(FIGS/"sis_I_curve.png", dpi=150); plt.close()

# Red de contagios acumulada
def draw_upto(t, fname):
    sub = inf[inf["t_infect"]<=t]
    G = nx.DiGraph(); G.add_nodes_from(nodes["id"])
    for r in sub.itertuples(): G.add_edge(r.source, r.target)
    pos = nx.spring_layout(G, seed=5)
    plt.figure(figsize=(7,5))
    nx.draw_networkx_nodes(G, pos, node_size=50)
    nx.draw_networkx_edges(G, pos, alpha=0.35, arrows=True)
    plt.title(f"Contagios acumulados hasta t={t}")
    plt.axis("off"); plt.tight_layout()
    plt.savefig(FIGS/fname, dpi=150); plt.close()

for tt in [40, 100, 200, 400]:
    draw_upto(tt, f"infections_upto_t{tt}.png")

print("Listo → figs/: sis_I_curve.png, infections_upto_t*.png")
