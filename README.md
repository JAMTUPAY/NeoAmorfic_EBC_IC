# NeoAmorfic_EBC_IC — Referee-fair replication (175 SPARC rotation curves)

This repository reproduces the two-track, referee-fair comparison in the paper:

> **An Entropy-Inspired Phenomenological Relation Competes with NFW on 175 SPARC Rotation Curves under Referee-Fair, Cross-Validated Tests**

- **Track A**: NFW with **no** cosmology prior  
- **Track B**: NFW with the **Dutton–Macciò (2014) c(M) prior** (DM14)

**Sample & results.** 175 inputs; 165 pass QC (≥6 usable points).  
- Track A winners: M0 62, M1 44, M2 55, M3 4; decisive: M0 45, M1 13, M2 10  
- Track B winners: M0 62, M1 66, M2 33, M3 4; decisive: M0 46, M1 22, M2 6

**Figure 1** (CV-WRMS boxplots) shows predictive accuracy (5-fold, radius-blocked CV; lower is better).  
Figures and tables are in `figs_trackA_dense/` and `figs_trackB_dense/`. The 12-panel overlays appear in `figs_m0_examples/`.

---

## Reproduction (dense grids, CV=5, σ-floor=3 km/s)

```bash
# Track A (no NFW prior)
python -u scripts/ebc_ic_referee.py \
  --data rc_parsed \
  --out  results_trackA_dense \
  --figs figs_trackA_dense \
  --sigma-floor 3.0 \
  --H0 70 \
  --cv-folds 5 \
  --nfw-prior none

# Track B (DM14 c(M) prior)
python -u scripts/ebc_ic_referee.py \
  --data rc_parsed \
  --out  results_trackB_dense \
  --figs figs_trackB_dense \
  --sigma-floor 3.0 \
  --H0 70 \
  --cv-folds 5 \
  --nfw-prior dm14
CV boxplots (Figure 1)
python -u scripts/make_cv_boxplots.py \
  --summary results_trackA_dense/ic_summary.csv \
  --out     figs_trackA_dense \
  --label   "Track A (dense)"

python -u scripts/make_cv_boxplots.py \
  --summary results_trackB_dense/ic_summary.csv \
  --out     figs_trackB_dense \
  --label   "Track B (dense, DM14 prior)"
Overlay montages (Appendix)
# Track A
python -u scripts/plot_examples_overlay.py \
  --data    rc_parsed \
  --summary results_trackA_dense/ic_summary.csv \
  --out     figs_m0_examples/trackA_dense_all \
  --k 12 --sigma-floor 3.0 --curves all --track "Track A (dense)"

# Track B (DM14 prior)
python -u scripts/plot_examples_overlay.py \
  --data    rc_parsed \
  --summary results_trackB_dense/ic_summary.csv \
  --out     figs_m0_examples/trackB_dense_all \
  --k 12 --sigma-floor 3.0 --curves all --dm14-prior --H0 70 --track "Track B (dense)"
QC lists: results_trackA_dense/exclusions.csv and results_trackB_dense/exclusions.csv.

⸻

Data

If SPARC redistribution is restricted in your jurisdiction, do not commit raw CSVs. Instead keep:
	•	rc_parsed/sha256s.txt (checksums) and
	•	a short “fetch” note in this README linking to SPARC.

Environment

Python ≥3.10; packages in requirements.txt.

License / Citation

MIT. Please cite the Zenodo DOI (badge to be added after release) and the paper PDF in paper/.
