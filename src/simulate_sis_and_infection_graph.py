import pandas as pd, networkx as nx, random
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
P    = ROOT/"data/processed"
RES  = ROOT/"results"; RES.mkdir(exist_ok=True)

nodes = pd.read_csv(P/"nodes.csv")
contacts = pd.read_csv(P/"contacts.csv")
T = int(contacts["t_end"].max()) if len(contacts)>0 else 0

beta  = 0.5   # contagio por ventana
delta = 0.02   # recuperación por ventana
seed_schedule = {0: 10, 40: 50, 120: 25}   # puedes editar libremente

by_t = {t:g for t,g in contacts.groupby("t_start")}
state = {v:0 for v in nodes["id"]}  # 0=S, 1=I
rng = random.Random(11)

def reseed(k):
    S = [v for v,s in state.items() if s==0]
    for v in rng.sample(S, min(k, len(S))):
        state[v]=1

if 0 in seed_schedule: reseed(seed_schedule[0])

infection_edges, series = [], []

for t in range(1, T+1):
    if t in seed_schedule and t!=0:
        reseed(seed_schedule[t])

    # recuperaciones
    for v in list(state.keys()):
        if state[v]==1 and rng.random()<delta:
            state[v]=0

    g = by_t.get(t)
    if g is not None and len(g)>0:
        Gt = nx.from_pandas_edgelist(g, "source","target")
        newly = []
        for v in Gt.nodes():
            if state[v]==0:
                inf_nbrs = [u for u in Gt.neighbors(v) if state[u]==1]
                if not inf_nbrs: continue
                # prob. ≈ 1 - (1-beta)^(w_total/20)
                w_total = g.query("source==@v or target==@v")["weight"].sum()
                blocks = max(1, int(round(w_total/20)))
                p = 1 - (1-beta)**blocks
                if rng.random() < p:
                    u = rng.choice(inf_nbrs)
                    newly.append((u,v))
        for (u,v) in newly:
            if state[v]==0:
                state[v]=1
                infection_edges.append((u,v,t))

    I = sum(1 for s in state.values() if s==1)
    S = len(state)-I
    series.append((t,S,I))

pd.DataFrame(infection_edges, columns=["source","target","t_infect"])\
  .assign(weight=lambda d: d["t_infect"], directed=True)\
  .to_csv(RES/"infection_edges.csv", index=False)

pd.DataFrame(series, columns=["t","S","I"]).to_csv(RES/"states_by_t.csv", index=False)

print("OK -> results/infection_edges.csv")
print("OK -> results/states_by_t.csv")
