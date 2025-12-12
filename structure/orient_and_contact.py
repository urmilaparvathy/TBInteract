import numpy as np
import pandas as pd
from sqlalchemy import create_engine
from Bio.PDB import PDBParser
import itertools

def load_ca(pdb_file):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("s", pdb_file)
    coords = []
    model = next(structure.get_models())
    for chain in model:
        for res in chain:
            if "CA" in res:
                coords.append(res["CA"].coord)
    return np.array(coords)

def center(coords):
    return coords - coords.mean(axis=0)

def rot_x(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[1,0,0],[0,c,-s],[0,s,c]])

def rot_y(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c,0,s],[0,1,0],[-s,0,c]])

def rot_z(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c,-s,0],[s,c,0],[0,0,1]])

def pairwise(a, b):
    return np.sqrt(((a[:,None,:] - b[None,:,:])**2).sum(axis=2))

db = snakemake.input.db
p1 = snakemake.input.pdb1
p2 = snakemake.input.pdb2
out = snakemake.output[0]

engine = create_engine(f"sqlite:///{db}")
edges = pd.read_sql("SELECT * FROM edges;", engine)

A = load_ca(p1)
B = load_ca(p2)
A0 = center(A)
B0 = center(B)

angles = np.deg2rad([0, 45, 90, 135, 180])

results = []

for _, row in edges.iterrows():
    u, v = row["protein1"], row["protein2"]
    score = float(row["combined_score"]) / 1000.0

    best_frac = -1
    best_mean = None

    for ax, ay, az in itertools.product(angles, angles, angles):
        R = rot_z(az) @ rot_y(ay) @ rot_x(ax)
        B_rot = B0 @ R.T

        D = pairwise(A0, B_rot)
        mean_d = float(D.mean())
        frac = float((D <= 8.0).sum() / D.size)

        if frac > best_frac:
            best_frac = frac
            best_mean = mean_d

    results.append({
        "protein1": u,
        "protein2": v,
        "contact_score": score,
        "best_interface_fraction": best_frac,
        "best_mean_distance": best_mean
    })

pd.DataFrame(results).to_csv(out, index=False)
print("Docking-like orientation completed:", out)
