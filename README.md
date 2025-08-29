# Referee-Fair Comparison of Galactic Rotation Curve Models (SPARC)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16996009.svg)](https://doi.org/10.5281/zenodo.16996009)

This repository contains the full code and analysis pipeline to reproduce the results in the paper:

> Tupay, J. A. M. (2025). *An Entropy-Inspired Phenomenological Relation Competes with NFW on 175 SPARC Rotation Curves under Referee-Fair, Cross-Validated Tests*. Zenodo. [https://doi.org/10.5281/zenodo.16996009](https://doi.org/10.5281/zenodo.16996009)

The analysis performs a head-to-head, referee-fair comparison of four models for describing galactic rotation curves on the SPARC dataset (175 galaxies).

---

## Models Compared

The comparison includes four models with strictly controlled parameter budgets:

* **M0: EBC (2 params)**: A simple entropy-inspired power law, $V(r) = Ar^{\alpha}$.
* **M1: Baryons + EBC (2 params)**: A composite model where the EBC component accounts for the missing velocity after considering baryons, $V^2(r) = V_{\text{bar}}^2(r) + (Ar^{\alpha})^2$.
* **M2: Baryons + NFW (2 params)**: The standard dark matter model, combining the baryonic component with an NFW halo.
* **M3: MOND (0 params)**: Modified Newtonian Dynamics, using the "simple" interpolation function with a fixed standard acceleration parameter $a_0$.

The comparison is run on two tracks:
* **Track A**: A purely phenomenological test where the NFW model has **no** cosmology prior.
* **Track B**: An ΛCDM-anchored test where the NFW model is constrained by the **Dutton–Macciò (2014) concentration-mass relation prior**.

---

## Key Results

Out of 175 galaxies, 165 passed quality control (≥6 usable data points per curve). The winning model for each galaxy was determined using the Bayesian Information Criterion (BIC).

| Track                   | M0 (EBC) Wins | M1 (Bar+EBC) Wins | M2 (Bar+NFW) Wins | M3 (MOND) Wins |
| ----------------------- | :-----------: | :---------------: | :---------------: | :------------: |
| **A (no NFW prior)** |      62       |        44         |        55         |       4        |
| **B (DM14 NFW prior)** |      62       |        66         |        33         |       4        |

*For a breakdown of decisive wins ($\Delta \text{BIC} > 10$), see the paper.*

---

## ⚙️ Reproduction Guide

### 1. Installation

Clone this repository and install the required Python packages.

```bash
git clone [https://github.com/your-username/NeoAmorfic_EBC_IC.git](https://github.com/your-username/NeoAmorfic_EBC_IC.git)
cd NeoAmorfic_EBC_IC
pip install -r requirements.txt
```

2. Data Setup

The analysis relies on the SPARC galaxy data.

Download the Data: Download the 175 _rotmod.csv files from the official SPARC website.

Place the Data: Move all 175 CSV files into the rc_parsed/ directory in this repository.

Verify Integrity (Optional): Run the following command to ensure all files are correct and match the versions used in the paper.

sha256sum -c rc_parsed/sha256s.txt

3. Running the Analysis

A single script is provided to run both tracks and generate all figures and tables.

3. Running the Analysis

A single script is provided to run both tracks and generate all figures and tables.

This script will execute the main analysis script (scripts/ebc_ic_referee.py) for both Track A and Track B, followed by the plotting scripts to generate all figures. It will take some time to complete.

(For advanced use, the individual commands can be found inside the run_all.sh script.)

📂 Repository Structure & Outputs
paper/: A copy of the manuscript PDF.

rc_parsed/: Directory for the input SPARC data CSVs.

scripts/: Contains all Python scripts for analysis and plotting.

results_trackA_dense/ & results_trackB_dense/: Output directories containing result tables, winner lists, and quality-control exclusion lists (exclusions.csv).

figs_trackA_dense/ & figs_trackB_dense/: Output directories for all figures, including the CV-WRMS boxplots.

figs_m0_examples/: Contains the 12-panel galaxy fit overlay plots shown in the appendix.

📜 Citation
If you use this code or the results in your research, please cite both the paper and this repository's Zenodo DOI.

@article{Tupay2025,
  author       = {Tupay, Johann Anton Michael},
  title        = {An Entropy-Inspired Phenomenological Relation Competes with NFW on 175 SPARC Rotation Curves under Referee-Fair, Cross-Validated Tests},
  year         = {2025},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.16996009},
  url          = {[https://doi.org/10.5281/zenod.16996009](https://doi.org/10.5281/zenod.16996009)}
}

⚖️ License
This project is licensed under the MIT License.
