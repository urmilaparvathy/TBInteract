rule all:
    input:
        "database/tb_protein_network.db",
        "results/diffusion_simulation.csv",
        "results/contact_scores.csv"
rule build_db:
    input:
        "data/processed/string_edges.csv"
    output:
        "database/tb_protein_network.db"
    script:
        "database/build_db.py"

rule diffusion:
    input:
        "database/tb_protein_network.db"
    output:
        "results/diffusion_simulation.csv"
    script:
        "dynamics/diffusion_model.py"
rule structure_contacts:
    input:
        db="database/tb_protein_network.db",
        pdb1="data/raw/alphafold_dprE1.pdb",
        pdb2="data/raw/alphafold_dprE2.pdb"
    output:
        "results/contact_scores.csv"
    script:
        "structure/contact_map.py"
rule oriented_contacts:
    input:
        db="database/tb_protein_network.db",
        pdb1="data/raw/alphafold_dprE1.pdb",
        pdb2="data/raw/alphafold_dprE2.pdb"
    output:
        "results/oriented_contact_scores.csv"
    script:
        "structure/orient_and_contact.py"
