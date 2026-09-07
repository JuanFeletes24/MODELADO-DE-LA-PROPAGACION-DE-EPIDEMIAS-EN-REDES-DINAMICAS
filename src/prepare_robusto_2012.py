# src/prepare_robusto_2012.py
# Normaliza data/raw/{nodes.csv, edges.csv} (Netzschleuder/SocioPatterns 2012)
# a:
#   data/processed/nodes.csv     (id,label,group)
#   data/processed/contacts.csv  (source,target,t_start,t_end,weight,directed)
#
# No asume cabeceras: toma la 1ª columna de nodes como ID;
# detecta automáticamente en edges:
#   - dos columnas de nodos (extremos)
#   - una columna de tiempo 't' (númerica, con rango amplio)
# Si no hay duración, usa 20 s.

import pandas as pd
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RAW  = ROOT/"data/raw"
OUT  = ROOT/"data/processed"
OUT.mkdir(parents=True, exist_ok=True)

def read_any(path):
    # Intenta autodetección de separador; si falla, coma
    try:
        return pd.read_csv(path, sep=None, engine="python")
    except Exception:
        return pd.read_csv(path)

# ---------- NODES ----------
nodes_path = RAW/"nodes.csv"
nodes_raw = read_any(nodes_path)

# Si no hay cabeceras o no coinciden, usa header=None
if nodes_raw.columns.size == 1 and nodes_raw.columns[0] == 0:
    nodes_df = nodes_raw.copy()
else:
    # Asegurar que la primera columna sea el ID, ignorando nombres
    nodes_df = nodes_raw.copy()

# Usar SIEMPRE la 1ª columna como ID (blindado)
first_col = nodes_df.columns[0]
id_series = nodes_df[first_col].astype(str)

# Intentar detectar una columna de grupo si existe (clase/sexo/etc.)
group_candidates = [c for c in nodes_df.columns if str(c).lower() in
                    ["class","group","curso","classe","gender","sexo"]]
group_series = nodes_df[group_candidates[0]].astype(str) if group_candidates else ""

nodes_out = pd.DataFrame({"id": id_series, "label": id_series, "group": group_series})
nodes_out.to_csv(OUT/"nodes.csv", index=False)

print(f"[NODES] id_col='{first_col}'  group_col='{group_candidates[0] if group_candidates else '(none)'}'  n={len(nodes_out)}")

# ---------- EDGES ----------
edges_path = RAW/"edges.csv"
edges_raw = read_any(edges_path)

# Si no hay cabeceras fiables, tratemos todo como sin cabecera (posición)
edges_df = edges_raw.copy()

# Identificar columnas numéricas candidatas a tiempo (rango amplio y valores grandes)
num_cols = [c for c in edges_df.columns if pd.api.types.is_numeric_dtype(edges_df[c])]

# Si ninguna es numérica, intentar convertir todo a numérico "suave"
if not num_cols:
    for c in edges_df.columns:
        try:
            edges_df[c] = pd.to_numeric(edges_df[c], errors="coerce")
        except Exception:
            pass
    num_cols = [c for c in edges_df.columns if pd.api.types.is_numeric_dtype(edges_df[c])]

if not num_cols:
    raise SystemExit("edges.csv: no hay columnas numéricas; revisa separadores o formato.")

# Heurística: la columna de tiempo tendrá
#  - muchos valores distintos y
#  - un valor máximo (p.ej., epoch) más alto que el de los ID
uniq_counts = {c: edges_df[c].nunique(dropna=True) for c in num_cols}
max_vals    = {c: edges_df[c].max(skipna=True)    for c in num_cols}

# Elegir como tiempo la columna con mayor 'uniq' y mayor 'max' relativo
time_col = sorted(num_cols, key=lambda c: (uniq_counts[c], max_vals[c]))[-1]

# Los extremos (nodos) serán las dos columnas más "intermedias" entre las restantes
rest = [c for c in edges_df.columns if c != time_col]
# Tomar las dos primeras de 'rest' como extremos por simplicidad robusta
if len(rest) < 2:
    raise SystemExit("edges.csv: no encuentro al menos 2 columnas además del tiempo para extremos.")
u_col, v_col = rest[0], rest[1]

df = edges_df[[u_col, v_col, time_col]].copy()
df.columns = ["i","j","t"]

# Quitar filas con t NaN y asegurar numérico
df["t"] = pd.to_numeric(df["t"], errors="coerce")
df = df.dropna(subset=["t"]).reset_index(drop=True)

# Duración: si hay una 4ª/5ª columna tipo 'dur/weight', úsala; si no, 20 s
dur_col = None
for cand in edges_df.columns:
    if str(cand).lower() in ["duration","dur","dt","delta","weight","w"]:
        dur_col = cand; break

if dur_col is not None:
    dur = pd.to_numeric(edges_df.loc[df.index, dur_col], errors="coerce").fillna(20.0).values
else:
    dur = np.full(len(df), 20.0)

# Discretizar a ventanas de 20 s
DELTA = 20
t0 = int(df["t"].min())
df["t_start"] = ((df["t"] - t0)//DELTA + 1).astype(int)
df["t_end"]   = ((df["t"] + dur - t0)//DELTA + 1).astype(int)
df["t_end"]   = df[["t_end","t_start"]].max(axis=1)
df["weight"]  = dur
df["directed"]= False

# Filtrar a nodos válidos
valid = set(nodes_out["id"])
before = len(df)
df = df[df["i"].astype(str).isin(valid) & df["j"].astype(str).isin(valid)]
after = len(df)

contacts_out = df.rename(columns={"i":"source","j":"target"})
contacts_out = contacts_out[["source","target","t_start","t_end","weight","directed"]]
contacts_out.to_csv(OUT/"contacts.csv", index=False)

print(f"[EDGES] time_col='{time_col}'  u_col='{u_col}'  v_col='{v_col}'  dur_col='{dur_col if dur_col else '(default 20s)'}'")
if after < before:
    print(f"[WARN] Se descartaron {before-after} filas de edges con IDs no presentes en nodes.")
print("OK -> data/processed/nodes.csv")
print("OK -> data/processed/contacts.csv")
print("T =", int(contacts_out["t_end"].max()) if len(contacts_out)>0 else 0)
