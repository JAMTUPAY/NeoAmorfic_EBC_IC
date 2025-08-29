#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Referee-safe IC/AICc/BIC recompute for 175 SPARC _rotmod galaxies.

Models:
  M0: EBC-only (V = A r^alpha)               [2 params]
  M1: Baryons + EBC                          [2 params]
  M2: Baryons + NFW (with optional c(M) prior) [2 params]
  M3: MOND "simple" nu, strict a0            [0 params]

Fairness options:
  - Symmetric coarse grid for M0, M1, M2 (default)
  - Optional 5-fold within-galaxy CV (default off, enable via --cv-folds 5)
  - NFW c(M) prior (Dutton & Macciò 2014) ON by default (--nfw-prior dm14)

Outputs:
  results/ic_summary.csv
  figs/winner_counts_bar.png
  figs/dBIC_M2_minus_M0.png
  figs/dBIC_M3_minus_M0.png
  figs/margin_of_victory_hist.png
  figs/ic_sanity_check.png
"""
import os, sys, math, csv, argparse, random
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ------------------------- Physical constants / helpers -------------------------
G_KPC = 4.30091e-6  # Newton's G in kpc * (km/s)^2 / Msun
KPC_TO_M = 3.0856775814913673e19
KM_TO_M  = 1.0e3

def dm14_log10c(M200_Msun, h):
    """
    Dutton & Macciò (2014) Planck z=0 c200(M200) for NFW fits:
      log10(c200) = a + b * log10( (M200 / (1e12 h^-1 Msun)) )
                  = a + b * ( log10(h*M200/Msun) - 12 )
    with a = 0.905, b = -0.101 (Table 3; NFW, 200c, z=0).
    Intrinsic scatter ~ 0.11 dex.
    Source: Dutton & Macciò 2014, MNRAS 441, 3359. 
    """
    a, b = 0.905, -0.101
    return a + b * (np.log10(max(h * M200_Msun, 1e-30)) - 12.0)

def safe_mkdir(path):
    os.makedirs(path, exist_ok=True)

def is_comment_or_blank(line):
    s = line.strip()
    return (not s) or s.startswith("#")

def parse_rotmod_row(line):
    """
    Robust parser for SPARC _rotmod lines:
    Expected common forms (numbers may be separated by spaces or commas):
       R  Vobs  eV  Vgas  Vdisk  Vbul
       R  Vobs  eV  Vgas  Vdisk
       R  Vobs  eV
    We treat commas as whitespace and split.
    Returns tuple of floats: (R, Vobs, eV, Vgas, Vdisk, Vbul)
    Missing components are set to 0.
    Raises ValueError for non-numeric rows.
    """
    line = line.replace(",", " ").replace("\t", " ")
    toks = [t for t in line.strip().split() if t]
    vals = [float(t) for t in toks]
    if len(vals) >= 6:
        R, Vobs, eV, Vgas, Vdisk, Vbul = vals[:6]
    elif len(vals) == 5:
        R, Vobs, eV, Vgas, Vdisk = vals
        Vbul = 0.0
    elif len(vals) == 4:
        # (rare) assume R, Vobs, Vgas, Vdisk
        R, Vobs, Vgas, Vdisk = vals
        eV = float("nan")
        Vbul = 0.0
    elif len(vals) == 3:
        R, Vobs, eV = vals
        Vgas = Vdisk = Vbul = 0.0
    else:
        raise ValueError("Too few numeric columns")
    return R, Vobs, eV, Vgas, Vdisk, Vbul

def load_galaxy_csv(path, sigma_floor):
    """
    Loads a single SPARC _rotmod CSV or whitespace file.
    Returns dict with arrays; filters bad rows; enforces sigma floor.
    """
    R, Vobs, eV, Vgas, Vdisk, Vbul = [], [], [], [], [], []
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if is_comment_or_blank(line):
                continue
            try:
                r, vo, ev, vg, vd, vb = parse_rotmod_row(line)
            except Exception:
                # skip non-data lines quietly
                continue
            R.append(r); Vobs.append(vo); eV.append(ev)
            Vgas.append(vg); Vdisk.append(vd); Vbul.append(vb)
    R = np.array(R, dtype=float)
    Vobs = np.array(Vobs, dtype=float)
    eV = np.array(eV, dtype=float)
    Vgas = np.array(Vgas, dtype=float)
    Vdisk= np.array(Vdisk, dtype=float)
    Vbul = np.array(Vbul, dtype=float)

    # Basic cleaning
    good = np.isfinite(R) & np.isfinite(Vobs)
    R, Vobs, eV, Vgas, Vdisk, Vbul = R[good], Vobs[good], eV[good], Vgas[good], Vdisk[good], Vbul[good]
    good = (R > 0) & np.isfinite(Vgas) & np.isfinite(Vdisk) & np.isfinite(Vbul)
    R, Vobs, eV, Vgas, Vdisk, Vbul = R[good], Vobs[good], eV[good], Vgas[good], Vdisk[good], Vbul[good]
    if R.size == 0:
        return None

    # sigma handling
    sig = np.where(np.isfinite(eV) & (eV > 0), eV, sigma_floor)
    sig = np.maximum(sig, sigma_floor)

    Vbar = np.sqrt(np.maximum(0.0, Vgas**2 + Vdisk**2 + Vbul**2))
    return dict(R=R, Vobs=Vobs, sig=sig, Vbar=Vbar)

# ------------------------- Models -------------------------
def model_M0(R, A, alpha):
    return A * np.power(R, alpha)

def model_M1(R, Vbar, A, alpha):
    return np.sqrt(np.maximum(0.0, Vbar**2 + (A*np.power(R, alpha))**2))

def v_nfw_kms(R_kpc, V200, c, H0):
    """
    Standard NFW circular velocity (km/s).
    x = r/R200; R200 = V200 / (10 H0) [Mpc] = V200/(10 H0)*1000 [kpc]
    V(r) = V200 * sqrt( [ln(1+c x) - c x/(1+c x)] / [ x (ln(1+c) - c/(1+c)) ] )
    """
    R200_kpc = (V200 / (10.0 * H0)) * 1000.0
    x = np.maximum(1e-12, R_kpc / max(R200_kpc, 1e-30))
    num = np.log(1.0 + c*x) - (c*x)/(1.0 + c*x)
    den = x * (np.log(1.0 + c) - c/(1.0 + c))
    f = np.maximum(1e-30, num/np.maximum(1e-30, den))
    return V200 * np.sqrt(f)

def model_M2(R, Vbar, V200, c, H0):
    Vh = v_nfw_kms(R, V200, c, H0)
    return np.sqrt(np.maximum(0.0, Vbar**2 + Vh**2))

def model_M3_strict_mond(R_kpc, Vbar_kms, a0=1.2e-10):
    """
    MOND 'simple' nu-function with strict a0 (in m s^-2).
      a_N = Vbar^2 / r  (SI units)
      nu(y) = 0.5 + sqrt(0.25 + 1/y)
      V = sqrt(nu) * Vbar
    """
    r_m = R_kpc * KPC_TO_M
    v_mps = Vbar_kms * KM_TO_M
    aN = np.maximum(1e-40, (v_mps**2) / np.maximum(1e-30, r_m))
    y = aN / a0
    nu = 0.5 + np.sqrt(0.25 + 1.0/np.maximum(1e-30, y))
    return np.sqrt(np.maximum(0.0, nu)) * Vbar_kms

# ------------------------- Fit / metrics -------------------------
def wrms(resid, sig):
    w = 1.0 / np.maximum(1e-30, sig**2)
    return math.sqrt(float(np.sum(w * resid**2) / np.sum(w)))

def chi2(resid, sig):
    return float(np.sum((resid / np.maximum(1e-30, sig))**2))

def aicc_from_chi2(ch2, k, n):
    if n <= (k + 1):
        return ch2 + 2*k
    return ch2 + 2*k + 2*k*(k+1)/(n - k - 1.0)

def bic_from_chi2(ch2, k, n):
    return ch2 + k * math.log(max(n,1))

def grid_best_M0(R, Vobs, sig, A_grid, alpha_grid):
    best = None
    for alpha in alpha_grid:
        Ra = np.power(R, alpha)
        for A in A_grid:
            V = A * Ra
            res = V - Vobs
            ch2 = chi2(res, sig)
            if (best is None) or (ch2 < best[0]):
                best = (ch2, A, alpha, wrms(res, sig))
    ch2, A, alpha, w = best
    n, k = len(R), 2
    return dict(wrms=w, chi2=ch2, aicc=aicc_from_chi2(ch2,k,n), bic=bic_from_chi2(ch2,k,n),
                pars=dict(A=A, alpha=alpha))

def grid_best_M1(R, Vobs, Vbar, sig, A_grid, alpha_grid):
    best = None
    for alpha in alpha_grid:
        Ra = np.power(R, alpha)
        for A in A_grid:
            V = np.sqrt(np.maximum(0.0, Vbar**2 + (A*Ra)**2))
            res = V - Vobs
            ch2 = chi2(res, sig)
            if (best is None) or (ch2 < best[0]):
                best = (ch2, A, alpha, wrms(res, sig))
    ch2, A, alpha, w = best
    n, k = len(R), 2
    return dict(wrms=w, chi2=ch2, aicc=aicc_from_chi2(ch2,k,n), bic=bic_from_chi2(ch2,k,n),
                pars=dict(A=A, alpha=alpha))

def mass200_from_V200(V200, H0):
    # M200 = V200^2 * R200 / G, with R200 = V200/(10 H0) * 1000 kpc
    R200_kpc = (V200 / (10.0 * H0)) * 1000.0
    return (V200*V200) * R200_kpc / G_KPC  # Msun

def grid_best_M2(R, Vobs, Vbar, sig, V200_grid, c_grid, H0, nfw_prior, sigma_logc):
    best = None
    h = H0/100.0
    ln_n = math.log(len(R))
    denom_c = (np.log(1.0 + c_grid) - c_grid/(1.0 + c_grid))  # cached term for speed (not strictly needed)
    for V200 in V200_grid:
        M200 = mass200_from_V200(V200, H0)   # Msun
        if nfw_prior == "dm14":
            log10c_pred = dm14_log10c(M200, h)
        for c in c_grid:
            V = model_M2(R, Vbar, V200, c, H0)
            res = V - Vobs
            ch2 = chi2(res, sig)
            if nfw_prior == "dm14":
                prior = ((math.log10(max(c,1e-30)) - log10c_pred)/sigma_logc)**2
                ch2_tot = ch2 + prior
            else:
                ch2_tot = ch2
            if (best is None) or (ch2_tot < best[0]):
                best = (ch2_tot, V200, c, wrms(res, sig))
    ch2, V200, c, w = best
    n, k = len(R), 2
    return dict(wrms=w, chi2=ch2, aicc=aicc_from_chi2(ch2,k,n), bic=bic_from_chi2(ch2,k,n),
                pars=dict(V200=V200, c=c))

def fit_M3(R, Vobs, Vbar, sig):
    V = model_M3_strict_mond(R, Vbar)
    res = V - Vobs
    ch2 = chi2(res, sig)
    n, k = len(R), 0
    return dict(wrms=wrms(res, sig), chi2=ch2, aicc=aicc_from_chi2(ch2,k,n), bic=bic_from_chi2(ch2,k,n),
                pars={})

# ------------------------- CV (optional) -------------------------
def kfold_indices(n, k):
    idx = np.arange(n)
    # deterministic fold assignment for reproducibility
    return [idx[i::k] for i in range(k)]

def cv_wrms_for_model(gal, model_name, grids, H0, nfw_prior, sigma_logc, cv_folds):
    if cv_folds <= 1: return float("nan")
    R, Vobs, sig, Vbar = gal['R'], gal['Vobs'], gal['sig'], gal['Vbar']
    folds = kfold_indices(len(R), cv_folds)
    wrms_list = []
    for i in range(cv_folds):
        test = folds[i]
        train = np.setdiff1d(np.arange(len(R)), test)
        Rtr, Vtr, Str, Btr = R[train], Vobs[train], sig[train], Vbar[train]
        Rte, Vte, Ste, Bte = R[test],  Vobs[test],  sig[test],  Vbar[test]
        if model_name == "M0":
            fit = grid_best_M0(Rtr, Vtr, Str, grids['A'], grids['alpha'])
            Vhat = model_M0(Rte, fit['pars']['A'], fit['pars']['alpha'])
        elif model_name == "M1":
            fit = grid_best_M1(Rtr, Vtr, Btr, Str, grids['A'], grids['alpha'])
            Vhat = model_M1(Rte, Bte, fit['pars']['A'], fit['pars']['alpha'])
        elif model_name == "M2":
            fit = grid_best_M2(Rtr, Vtr, Btr, Str, grids['V200'], grids['c'], H0, nfw_prior, sigma_logc)
            Vhat = model_M2(Rte, Bte, fit['pars']['V200'], fit['pars']['c'], H0)
        else: # M3
            Vhat = model_M3_strict_mond(Rtr, Btr)  # parameters free? none. We can just use full-data; use train here
            # Recompute on test using same formula (no parameters)
            Vhat = model_M3_strict_mond(Rte, Bte)
        wrms_list.append(wrms(Vhat - Vte, Ste))
    return float(np.nanmean(wrms_list)) if wrms_list else float("nan")

# ------------------------- Figures -------------------------
def winner_bar(winners, outpng, outsvg=None):
    lbls = ["M0","M1","M2","M3"]
    vals = [winners.get(k,0) for k in lbls]
    plt.figure(figsize=(6,3.8))
    plt.bar(lbls, vals)
    plt.ylabel("Winners (by BIC)")
    plt.title("Per-galaxy winners")
    for i,v in enumerate(vals):
        plt.text(i, v+0.5, str(v), ha='center', va='bottom', fontsize=9)
    plt.tight_layout()
    plt.savefig(outpng, dpi=150)
    if outsvg: plt.savefig(outsvg)
    plt.close()

def plot_hist(data, outpng, xlabel, title, bins=40, xlim=None):
    x = np.array([d for d in data if np.isfinite(d)])
    plt.figure(figsize=(6.4,4.0))
    plt.hist(x, bins=bins)
    if xlim: plt.xlim(*xlim)
    plt.xlabel(xlabel); plt.ylabel("Count"); plt.title(title)
    plt.tight_layout(); plt.savefig(outpng, dpi=150); plt.close()

def sanity_check_plot(outpng):
    plt.figure(figsize=(6.4,3.0))
    plt.axis('off')
    plt.text(0.02,0.7,"IC sanity check:", fontsize=12, weight='bold')
    plt.text(0.02,0.5,"BIC = χ² + k ln n;  AICc = χ² + 2k + 2k(k+1)/(n-k-1)", fontsize=11)
    plt.text(0.02,0.3,"NFW prior: Dutton & Macciò 2014 (Planck, z=0):  log10 c = 0.905 - 0.101 log10(M200/(1e12 h^{-1} Msun)), σ=0.11 dex", fontsize=10)
    plt.tight_layout(); plt.savefig(outpng, dpi=150); plt.close()

# ------------------------- Main -------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", required=True, help="Folder with *_rotmod*.csv files")
    ap.add_argument("--out",  required=True, help="Output folder for CSV")
    ap.add_argument("--figs", required=True, help="Output folder for figures")
    ap.add_argument("--sigma-floor", type=float, default=3.0)
    ap.add_argument("--H0", type=float, default=70.0)
    ap.add_argument("--cv-folds", type=int, default=1)
    ap.add_argument("--min-points", type=int, default=6)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--bins", type=int, default=40)
    ap.add_argument("--nfw-prior", choices=["dm14","none"], default="dm14")
    ap.add_argument("--sigma-logc", type=float, default=0.11, help="scatter (dex) for c(M) prior")
    args = ap.parse_args()

    safe_mkdir(args.out); safe_mkdir(args.figs)

    # Grids (symmetric across models)
    A_grid     = np.linspace(1.0, 300.0, 60)     # EBC amplitude
    alpha_grid = np.linspace(0.1, 1.20, 40)      # EBC slope
    V200_grid  = np.linspace(30.0, 350.0, 40)    # NFW scale vel (km/s)
    c_grid     = np.linspace(2.0, 30.0, 32)      # NFW concentration

    grids = dict(A=A_grid, alpha=alpha_grid, V200=V200_grid, c=c_grid)

    # Collect galaxy files (rotmod only)
    candidates = []
    for fn in sorted(os.listdir(args.data)):
        if not fn.lower().endswith(".csv"): continue
        if "_rotmod" not in fn.lower():     continue
        candidates.append(os.path.join(args.data, fn))
    if not candidates:
        print("No *_rotmod*.csv files found")
        sys.exit(1)
    if args.limit and args.limit > 0:
        candidates = candidates[:args.limit]

    rows = []
    winners = {}
    dBIC_M2M0, dBIC_M3M0, dBIC2 = [], [], []

    for i, path in enumerate(candidates, 1):
        gname = os.path.splitext(os.path.basename(path))[0]
        gal = load_galaxy_csv(path, args.sigma_floor)
        if gal is None or len(gal['R']) < args.min_points:
            print(f"[WARN] {gname}.csv skipped: too few valid points")
            continue

        R, Vobs, Vbar, sig = gal['R'], gal['Vobs'], gal['Vbar'], gal['sig']

        # Fit all models
        fit0 = grid_best_M0(R, Vobs, sig, A_grid, alpha_grid)
        fit1 = grid_best_M1(R, Vobs, Vbar, sig, A_grid, alpha_grid)
        fit2 = grid_best_M2(R, Vobs, Vbar, sig, V200_grid, c_grid, args.H0, args.nfw_prior, args.sigma_logc)
        fit3 = fit_M3(R, Vobs, Vbar, sig)

        # Optional CV WRMS
        cv0 = cv_wrms_for_model(gal, "M0", grids, args.H0, args.nfw_prior, args.sigma_logc, args.cv_folds)
        cv1 = cv_wrms_for_model(gal, "M1", grids, args.H0, args.nfw_prior, args.sigma_logc, args.cv_folds)
        cv2 = cv_wrms_for_model(gal, "M2", grids, args.H0, args.nfw_prior, args.sigma_logc, args.cv_folds)
        cv3 = cv_wrms_for_model(gal, "M3", grids, args.H0, args.nfw_prior, args.sigma_logc, args.cv_folds)

        # Winner by BIC
        bics = dict(M0=fit0['bic'], M1=fit1['bic'], M2=fit2['bic'], M3=fit3['bic'])
        winner = min(bics, key=bics.get)
        winners[winner] = winners.get(winner, 0) + 1
        bic_sorted = sorted(bics.items(), key=lambda kv: kv[1])
        margin = bic_sorted[1][1] - bic_sorted[0][1]  # second-best minus best

        dBIC_M2M0.append(fit2['bic'] - fit0['bic'])
        dBIC_M3M0.append(fit3['bic'] - fit0['bic'])
        dBIC2.append(margin)

        rows.append([
            gname, len(R),
            fit0['wrms'], fit0['aicc'], fit0['bic'],
            fit1['wrms'], fit1['aicc'], fit1['bic'],
            fit2['wrms'], fit2['aicc'], fit2['bic'],
            fit3['wrms'], fit3['aicc'], fit3['bic'],
            winner, margin, fit2['bic']-fit0['bic'], fit3['bic']-fit0['bic'],
            cv0, cv1, cv2, cv3
        ])

        print(f"[{i}/{len(candidates)}] {gname} → {winner}")

    # Write CSV
    out_csv = os.path.join(args.out, "ic_summary.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["galaxy","n",
                    "wrms_M0","aicc_M0","bic_M0",
                    "wrms_M1","aicc_M1","bic_M1",
                    "wrms_M2","aicc_M2","bic_M2",
                    "wrms_M3","aicc_M3","bic_M3",
                    "winner","dBIC_second_minus_best","dBIC_M2_minus_M0","dBIC_M3_minus_M0",
                    "cv_wrms_M0","cv_wrms_M1","cv_wrms_M2","cv_wrms_M3"])
        for r in rows: w.writerow(r)

    # Figures
    winner_bar(winners,
        os.path.join(args.figs, "winner_counts_bar.png"),
        outsvg=os.path.join(args.figs, "winner_counts_bar.svg"))

    a = np.array(dBIC_M2M0); b = np.array(dBIC_M3M0); m = np.array(dBIC2)
    def safemin(x, default):
        x = x[np.isfinite(x)]
        return float(np.min(x)) if x.size else default
    def safemax(x, default):
        x = x[np.isfinite(x)]
        return float(np.max(x)) if x.size else default

    plot_hist(a, os.path.join(args.figs, 'dBIC_M2_minus_M0.png'),
              xlabel='ΔBIC (NFW − EBC)', title='ΔBIC: NFW vs EBC', bins=40,
              xlim=(min(-10.0, safemin(a, -10.0)), max(50.0, safemax(a, 50.0))))
    plot_hist(b, os.path.join(args.figs, 'dBIC_M3_minus_M0.png'),
              xlabel='ΔBIC (MOND(strict) − EBC)', title='ΔBIC: MOND(strict) vs EBC', bins=40,
              xlim=(min(-10.0, safemin(b, -10.0)), max(50.0, safemax(b, 50.0))))
    vmax = np.nanpercentile(m, 99) if len(m)>0 else 10.0
    plot_hist(m, os.path.join(args.figs, 'margin_of_victory_hist.png'),
              xlabel='ΔBIC (second-best − best)', title='Margin of victory (Jeffreys scale)', bins=40,
              xlim=(0.0, max(10.0, float(vmax))))
    sanity_check_plot(os.path.join(args.figs, 'ic_sanity_check.png'))

    print("Winner counts:", winners)
    print("Output CSV:", out_csv)
    print("Figures:", os.path.join(args.figs, 'winner_counts_bar.png'))

if __name__ == '__main__':
    main()
