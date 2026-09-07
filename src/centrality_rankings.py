# src/centrality_rankings.py
import pandas as pd, networkx as nx
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
P    = ROOT/"data/processed"
R    = ROOT/"results"; R.mkdir(exist_ok=True)

# --- Red de contagios (dirigida) ---
inf = pd.read_csv(R/"infection_edges.csv").astype({"source":str,"target":str})
Ginf = nx.DiGraph()
Ginf.add_edges_from(inf[["source","target"]].itertuples(index=False, name=None))

outdeg = dict(Ginf.out_degree())
pr     = nx.pagerank(Ginf, alpha=0.85) if Ginf.number_of_nodes()>0 else {}
btw    = nx.betweenness_centrality(Ginf, k=min(200, Ginf.number_of_nodes()), seed=42) if Ginf.number_of_nodes()>0 else {}

rank_inf = pd.DataFrame({
    "id": list(Ginf.nodes()),
    "out_degree": [outdeg.get(n,0) for n in Ginf.nodes()],
    "pagerank":   [pr.get(n,0.0) for n in Ginf.nodes()],
    "betw_approx":[btw.get(n,0.0) for n in Ginf.nodes()],
}).sort_values(["out_degree","pagerank"], ascending=False)
rank_inf.to_csv(R/"rank_infection_centrality.csv", index=False)

# --- Grafo de contactos agregado (no temporal, no dirigido) ---
contacts = pd.read_csv(P/"contacts.csv").astype({"source":str,"target":str})
Gcont = nx.from_pandas_edgelist(contacts, "source","target", create_using=nx.Graph())
deg_cont = dict(Gcont.degree())
rank_cont = pd.DataFrame({"id": list(Gcont.nodes()), "degree_contact": [deg_cont[n] for n in Gcont.nodes()]})\
             .sort_values("degree_contact", ascending=False)
rank_cont.to_csv(R/"rank_contact_degree.csv", index=False)

print("OK -> results/rank_infection_centrality.csv")
print("OK -> results/rank_contact_degree.csv")
