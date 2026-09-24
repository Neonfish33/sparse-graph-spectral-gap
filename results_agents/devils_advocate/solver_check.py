"""Cross-check eigsh(giant) vs dense eigh on saved n=3000 graphs (dense>threshold)."""
import os, sys, glob
import numpy as np
import networkx as nx

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
G_DIR = os.path.join(ROOT, "results_v2", "graphs")

def dense_lambda2(G):
    nodes = list(G.nodes())
    A = nx.to_numpy_array(G, nodelist=nodes, dtype=float)
    deg = A.sum(1)
    inv = np.zeros_like(deg); nz = deg > 0; inv[nz] = 1/np.sqrt(deg[nz])
    L = np.eye(len(nodes)) - (inv[:, None] * A * inv[None, :])
    vals = np.linalg.eigvalsh(L)
    return vals[:4]

for rep in range(3):
    p = os.path.join(G_DIR, f"pl_dmin2_g2.5_n3000_d2_r{rep}.npz")
    d = np.load(p)
    n = int(d["n"][0]); edges = d["edges"]
    G = nx.Graph(); G.add_nodes_from(range(n)); G.add_edges_from(map(tuple, edges.tolist()))
    giant_nodes = max(nx.connected_components(G), key=len)
    giant = G.subgraph(giant_nodes).copy()
    vals = dense_lambda2(giant)
    print(f"rep={rep} n_full={n} giant_n={giant.number_of_nodes()} dense_lam2={vals[1]:.8f} "
          f"lam3={vals[2]:.8f} lam4={vals[3]:.8f}")
