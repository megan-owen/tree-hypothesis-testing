# Tree Hypothesis Testing Project (`summer2025`)

This project implements a complete statistical pipeline for comparing two sets of **gene trees** using **Fréchet means**, **difference-in-means permutation tests**, and **cross-match tests** under both **weighted Robinson–Foulds (wRF)** and **BHV (Billera–Holmes–Vogtmann)** distance metrics.

---

## 🧩 Overview of Components

### Java Tools
- **`SturmMean_201102.jar`** – Computes Fréchet (Sturm) means of a set of trees.  
- **`gtp.jar`** – Computes BHV tree-space geodesics (used via JPype or command line).

### Python Scripts
- **`trees_to_edge_lengths.py`** – Simulates gene trees using DendroPy, computes edge lengths via Fréchet mean (SturmMean), and writes outputs.  
- **`crossmatch_functions.py`** – Implements the **Cross-Match test** on a set of items given a distance matrix.   
- **`crossmatch_bhv.py`** – Runs Cross-Match test using **BHV distances** computed via `gtp.jar` through JPype.  
- **`crossmatch_sample.py`** – Generalized script for performing Cross-Match tests on custom datasets.  
- **`eric_functions.py`** – Utility functions for additional analyses.  

### Directories
- **`BHV/`** – Practice scripts for JPype calls to Java classes from Python.  
- **`Locus Analysis/`** – Visualization suite that aggregates and compares test results across loci, producing summary graphs and heatmaps.  
- **`batch_output/`** – output directory replacing the old `mean_test_results/`.  
- **`crossmatch_results_summary/`** – Contains per-locus Cross-Match outputs and summary CSVs.  
- **`rootedtrees/`** – Input dataset (Dryad) with 304 files, each containing 1000 rooted trees.  

---

## 📊 Test Overview

Both **wRF** and **BHV** Cross-Match tests were benchmarked on locus trees.  
The **Robinson–Foulds** implementation was ultimately chosen for efficiency and comparable runtime to the BHV approach.

---

## 📁 Output Summary

### Mean Test Results
- The folder `mean_test_results/` (now replaced by `batch_output/`) contains permutation test outputs across **304 loci**.
- **Graphs:** Histograms available for loci **1–304**.  
- **CSV Summary:** `new_permutation_summary.csv` includes all loci (1–304), summarizing **p-values** and **test statistics** comparing the two gene tree sets.

---

## ⚙️ How to Run the Pipeline

### 1️⃣ Generate Gene Trees

Use `trees_to_edge_lengths.py` (accepts `argparse` parameters):

```bash
python trees_to_edge_lengths.py --tau1 10000 --tau2 10500 --tau3 11000 \
--N1 10000 --N2 10000 --N3 10000 --out output/

python trees_to_edge_lengths.py --tau1 10000 --tau2 10100 --tau3 11100 \
--N1 10000 --N2 10000 --N3 10000 --out output/
```

**Example Output Files:**
```
output/gts_dendropy_CAT_tauAB-10000_tauABC-10500_tauRoot-11000_pAB-10000_pABC-10000_pRoot-10000
output/gts_dendropy_CAT_tauAB-10000_tauABC-10100_tauRoot-11100_pAB-10000_pABC-10000_pRoot-10000
```

---

## 🧠 Notes

- The project’s current data and analysis results are located in:
  - `batch_output/` → mean test results  
  - `crossmatch_results_summary/` → per-locus cross-match statistics  
  - `Locus Analysis/` → scripts for comparing test outcomes across loci



