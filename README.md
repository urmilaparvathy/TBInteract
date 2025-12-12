# 🧬 TBInteract  
*A reproducible network + structural modeling pipeline for Mycobacterium tuberculosis protein interactions*

TBInteract is a research-grade computational pipeline that integrates:

- **STRING PPI networks** for *M. tuberculosis*  
- **Diffusion-based influence propagation** (DprE1 as seed)  
- **AlphaFold structural analysis**  
- **CA–CA contact scoring**  
- **A quasi-docking rotational search** to approximate protein–protein interfaces  
- **Snakemake workflow** for full reproducibility  

This project demonstrates reproducible computational biology practices useful for structural bioinformatics, network biology, antibiotic target analysis, and systems-level TB research.

---

## 🚀 Features

### **1. Network dynamics**
- Loads STRING network into SQLite  
- Performs deterministic **heat diffusion** from a seed node (DprE1)  
- Outputs ranked influence scores + JSON summary  

### **2. Structural Contact Analysis**
- Parses AlphaFold models  
- Calculates CA–CA distances and interface fractions  
- Generates interpretable structural metrics  

### **3. Quasi-Docking Orientation Algorithm**
- Rotates protein B in 3D across predefined angles  
- Detects best possible interaction orientation  
- Reports best interface fraction, mean distance & contact scoring  

---

## 📁 Workflow Overview (Snakemake)
```
TBInteract/
├── data/
│   ├── raw/ (AlphaFold models)
│   └── processed/
├── database/ (SQLite DB)
├── dynamics/ (diffusion model)
├── structure/ (contact + quasi-docking)
├── results/
├── Snakefile
└── requirements.txt
```
