import pandas as pd
import networkx as nx
from sqlalchemy import create_engine
import json
import os

db_path = snakemake.input[0]
out_csv = snakemake.output[0]
# also write a JSON summary next to the CSV
out_json = os.path.splitext(out_csv)[0] + "_summary.json"
engine = create_engine(f"sqlite:///{db_path}")
edges = pd.read_sql("SELECT * FROM edges;", engine)

# Build weighted graph (normalize combined_score to [0,1] if present)
G = nx.Graph()
for _, r in edges.iterrows():
    u = r["protein1"]
    v = r["protein2"]
    w = float(r.get("combined_score", 0)) / 1000.0 if "combined_score" in r and pd.notna(r["combined_score"]) else 1.0
    G.add_edge(u, v, weight=w)
# Seed: DprE1 (use a best-effort match)
seed_candidates = ["dprE1", "DprE1", "dpre1", "RvXXX"]  # adjust if your IDs differ
seed_node = None
for s in seed_candidates:
    if s in G:
        seed_node = s
        break
# fallback to the first node if none of the candidates match
if seed_node is None:
    seed_node = next(iter(G.nodes()))
# initialize scores
scores = {n: 0.0 for n in G.nodes()}
scores[seed_node] = 1.0

# diffusion parameters
alpha = 0.5
steps = 30
for step in range(steps):
    new_scores = scores.copy()
    for n in G.nodes():
        neigh = list(G[n])
        if not neigh:
            continue
        s = 0.0
        total_w = 0.0
        for nbr in neigh:
            w = G[n][nbr].get("weight", 1.0)
            s += scores[nbr] * w
            total_w += w
        if total_w > 0:
            new_scores[n] = (1 - alpha) * scores[n] + alpha * (s / total_w)
    scores = new_scores
# Normalize scores to [0,1]
min_s = min(scores.values())
max_s = max(scores.values())
if max_s - min_s > 0:
    norm_scores = {k: (v - min_s) / (max_s - min_s) for k, v in scores.items()}
else:
    norm_scores = {k: 0.0 for k in scores.keys()}

# Write CSV (sorted descending)
df = pd.DataFrame.from_dict(norm_scores, orient="index", columns=["score"]).reset_index().rename(columns={"index":"protein"})
df = df.sort_values("score", ascending=False)
df.to_csv(out_csv, index=False)
# Prepare JSON summary: top-10 proteins and basic stats
topk = df.head(10).to_dict(orient="records")
summary = {
    "seed_node": seed_node,
    "steps": steps,
    "alpha": alpha,
    "n_nodes": len(G.nodes()),
    "n_edges": len(G.edges()),
    "top_influencers": topk
}
with open(out_json, "w") as fh:
    json.dump(summary, fh, indent=2)
# Print top influencers to console
print("Diffusion done. Seed node:", seed_node)
print("Top influencers (protein : score):")
for r in topk:
    print(f'  {r["protein"]} : {r["score"]:.4f}')

print("Wrote:", out_csv)
print("Wrote:", out_json)
