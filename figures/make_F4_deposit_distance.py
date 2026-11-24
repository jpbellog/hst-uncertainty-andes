# -*- coding: utf-8 -*-
"""Figure 4 -- Deposit-track distance against the uncertainty of the track.

Shows the distribution of the deposits as a function of their distance to the
track, compared with the uncertainty band of the track position (in the style
of Fig. 10 of Bello-Gonzalez et al. 2026).

(A) Empirical CDF of the restored distance of each deposit to the track of its
    hotspot (porphyries, epithermal deposits and igneous magmatic background),
    with the MEDIAN 68% and 95% half-width of the collision band (in km) marked
    as vertical lines: the fraction of the curve to the left of each line is the
    fraction of cases whose N/S assignment is not resolved above the
    uncertainty of the track itself.
(B) Distance of each deposit to the track against its age, over the 68% (fill)
    and 95% (line) half-width of the collision band of the corresponding track,
    in km.

Data: the same files as step 07 of the pipeline (deposit database, v7
porphyries and v6 epithermal deposits; geochemistry database v12_geodyn2; age
quality filter v3) and the Monte Carlo bands of this paper
(models/Q22_WK08-A/mc_collision_bands_total.csv, f_PAC = 0.44).

Output: figures/F4_deposit_distance.pdf/.png
Usage:  py make_F4_deposit_distance.py
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
BANDS = os.path.join(os.path.dirname(HERE), "Q22_WK08-A",
                     "mc_collision_bands_total.csv")
# The two circuit versions are overlaid. The deposit distances are the SAME in both
# models (they come from the database); only the bands change, so the revised
# version is drawn as a dashed line without fill.
BANDS_WD = os.path.join(os.path.dirname(HERE), "WD25_WK08-A",
                        "mc_collision_bands_total.csv")
DEP = os.path.join(__import__("os").environ["PAPER_ROOT"], "inputs", "deposits")  # copy of the deposit database (see rotation_models/README.md)
PORF = os.path.join(DEP, "Porphyry_Plate_ID_geodyn_v7.csv")
EPIT = os.path.join(DEP, "Epithermal_Plate_Id_geodyn_v6_thesis.csv")
GEO = os.path.join(__import__("os").environ["PAPER_ROOT"], "inputs", "geochemistry",
      "Full_Data_v11_tt_with_Clusters_v12_geodyn2.xlsx")
KM = 111.19492664455873
PF = "0.44"
HS2RIDGE = {13: "JFR", 14: "Easter", 15: "TT", 16: "Foundation"}
RIDGE_LABEL = {"Easter": "Easter", "TT": "San Félix", "JFR": "Juan Fernández",
               "Foundation": "Foundation"}
RIDGE_COLOR = {"Easter": "#d95f02", "TT": "#7570b3", "JFR": "#1b9e77", "Foundation": "#8c510a"}

# With PAPER_FIG_TAG (tag of the model outputs, e.g. _hsb_cc_g13) the figure reads the step 07
# outputs of PAPER_PIPE (deposit_tracks_*.csv: nearest of the four chains, HS4, and its
# distance, Dist4_km) with the background of step 12 (method filter, 0-63.5 Ma), and the bands
# with that tag from PAPER_PIPE and PAPER_PIPE_WD. Without it, the figure uses the HS column
# of the deposit database (three chains).
TAG = os.environ.get("PAPER_FIG_TAG", "")
ROOT = os.path.dirname(HERE)
PIPE = os.path.join(ROOT, "models", os.environ.get("PAPER_PIPE", "Q22_WK08-A"))
FIG_OUT = os.environ.get("PAPER_FIG_OUT", HERE)
if TAG:
    BANDS = os.path.join(PIPE, f"mc_collision_bands_total{TAG}.csv")
    BANDS_WD = os.path.join(ROOT, "models", os.environ.get("PAPER_PIPE_WD", "WD25_WK08-A"),
                            f"mc_collision_bands_total{TAG}.csv")

# age quality filter v3 (identical to step 07)
EXCLUDED = {m.strip().lower() for m in [
    "K/Ar plagioclase", "K/Ar amphibole", "inference from K/Ar plagioclase",
    "Inference", "inference from K/Ar whole rock", "K/Ar whole rock",
    "Inference from K/Ar (biotite plagioclase whole rock) data",
    "Inference from K/Ar (biotite amphibole) data",
    "Inference from K/Ar biotite data", "K/Ar biotite", "K/Ar feldespar",
    "K/Ar glass", "K/Ar", "N.D.", "Inference from K/Ar",
    "The age of sample S-16 is assigned", "?",
    "The age of sample S-59 is assigned",
    "The age of sample S-H22 (176m) is assigned",
    "The age of sample S-T2 is assigned",
    "The age of sample S-50 is assigned",
    "The age of sample S-31 is assigned",
    "K-Ar feldspar", "K-Ar", "inference with dated units", "K-Ar biotite",
    "K-Ar whole-rock", "K-Ar amphibole", "K-Ar hornblende", "K-Ar plagioclase",
    "radiometric", "K-Ar whole rock",
    "inference from dated samples of same unit with K-Ar"]}

plt.rcParams.update({"font.family": "serif", "mathtext.fontset": "cm",
                     "font.size": 9.5, "axes.linewidth": 0.8})


def load_bands(path=BANDS):
    d = pd.read_csv(path, dtype={"pf": str})
    d = d[d.pf == PF].sort_values(["ridge", "age"])
    out = {}
    # Foundation does not collide at all ages: without those rows the interpolation gives no NaN
    d = d.dropna(subset=["lat_p2_5", "lat_p16", "lat_p84", "lat_p97_5"])
    for ridge, t in d.groupby("ridge"):
        out[ridge] = dict(age=t.age.values,
                          h68=((t.lat_p84 - t.lat_p16) / 2).values * KM,
                          h95=((t.lat_p97_5 - t.lat_p2_5) / 2).values * KM)
    return out


def prep(df, age_col):
    df = df.copy()
    df["age"] = pd.to_numeric(df[age_col], errors="coerce")
    hs, dcol = ("HS4", "Dist4_km") if TAG else ("HS", "Dist_Sample_Track")
    df["dist"] = pd.to_numeric(df[dcol], errors="coerce")
    df["ridge"] = df[hs].map(HS2RIDGE)
    df = df.dropna(subset=["age", "dist", "ridge"])
    return df[(df.age >= 0) & (df.age <= 66)]


def ecdf(x):
    x = np.sort(np.asarray(x, float))
    return x, np.arange(1, len(x) + 1) / len(x)


def main():
    bands = load_bands()
    bands_wd = load_bands(BANDS_WD)
    if TAG:
        por = prep(pd.read_csv(os.path.join(PIPE, "deposit_tracks_porphyries.csv")), "Age MY")
        epi = prep(pd.read_csv(os.path.join(PIPE, "deposit_tracks_epithermal.csv")), "Age MY")
        g = pd.read_csv(os.path.join(PIPE, "deposit_tracks_samples.csv"), low_memory=False)
    else:
        por = prep(pd.read_csv(PORF), "Age MY")
        epi = prep(pd.read_csv(EPIT), "Age MY")
        g = pd.read_excel(GEO)
    meth = g["Method"].astype(str).str.strip().str.lower()
    age = pd.to_numeric(g["Age (Ma)"], errors="coerce")
    bad = meth.isin(EXCLUDED) | g["Method"].isna() | (meth == "") | (meth == "nan")
    g = g[(~bad) | (age < 1)]
    bg = prep(g, "Age (Ma)")
    if TAG:      # background of step 12
        bg = bg[bg.age <= 63.5]
    else:
        bg = bg[(bg["Latitude_restored"] <= -2.9) & (bg["Latitude_restored"] >= -34.5)
                & (bg.age <= 63.5)]

    def halfwidth(df, what, bb=None):
        bb = bands if bb is None else bb
        return np.array([np.interp(a, bb[r]["age"], bb[r][what])
                         for r, a in zip(df.ridge, df.age)])

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.3, 3.5))

    # ---- (A) CDF of the distance to the track ----
    series = [(por, "Porphyry deposits", "#b2182b", "-", 1.9),
              (epi, "Epithermal deposits", "#2166ac", (0, (5, 2)), 1.6),
              (bg, "Igneous samples", "#555555", (0, (1, 1.4)), 1.5)]
    for df, label, color, ls, lw in series:
        x, y = ecdf(df.dist)
        axA.plot(x, y, color=color, ls=ls, lw=lw, zorder=3,
                 label=f"{label} (n = {len(x)})")
    h68m = float(np.median(halfwidth(por, "h68")))
    h95m = float(np.median(halfwidth(por, "h95")))
    axA.axvspan(0, h68m, color="#bdbdbd", alpha=0.45, zorder=0)
    axA.axvline(h95m, color="#777777", lw=0.9, ls=":", zorder=1)
    # revised circuit version: same medians, marked with a vertical line
    h68w = float(np.median(halfwidth(por, "h68", bands_wd)))
    h95w = float(np.median(halfwidth(por, "h95", bands_wd)))
    for xv in (h68w, h95w):
        axA.axvline(xv, color="#333333", lw=1.0, ls=(0, (2, 1.6)), zorder=2)
    axA.text(h68m - 10, 0.985, "median 68%\nhalf-width", fontsize=7.2,
             ha="right", va="top", color="#444444")
    axA.text(h95m + 10, 0.985, "median 95%\nhalf-width", fontsize=7.2,
             ha="left", va="top", color="#666666")
    axA.set_xlim(0, 900)
    axA.set_ylim(0, 1)
    axA.set_xlabel("Restored distance to nearest hotspot track (km)")
    axA.set_ylabel("Cumulative fraction")
    axA.grid(True, color="#dddddd", lw=0.5, zorder=0)
    axA.set_axisbelow(True)
    h, l = axA.get_legend_handles_labels()
    h.append(plt.Line2D([0], [0], color="#333333", lw=1.0, ls=(0, (2, 1.6))))
    l.append("Same half-widths, revised circuit")
    axA.legend(h, l, loc="lower right", frameon=False, fontsize=7.8)
    axA.set_title("A", loc="left", fontweight="bold", fontsize=10)

    # ---- (B) distance vs age, over the bands ----
    for ridge in [r for r in ("Easter", "TT", "JFR", "Foundation") if r in bands]:
        b = bands[ridge]
        axB.fill_between(b["age"], 0, b["h68"], color=RIDGE_COLOR[ridge],
                         alpha=0.16, lw=0, zorder=1)
        axB.plot(b["age"], b["h95"], color=RIDGE_COLOR[ridge], lw=1.0,
                 ls=(0, (4, 2)), zorder=2, label=f"{RIDGE_LABEL[ridge]} 95%")
        axB.plot(b["age"], b["h68"], color=RIDGE_COLOR[ridge], lw=1.2,
                 zorder=2)
        bw = bands_wd[ridge]
        axB.plot(bw["age"], bw["h68"], color=RIDGE_COLOR[ridge], lw=1.0,
                 ls=(0, (2, 1.6)), zorder=3)
    axB.scatter(bg.age, bg.dist, s=4, c="#bbbbbb", lw=0, zorder=3,
                label="Igneous samples")
    axB.scatter(epi.age, epi.dist, s=18, marker="x", c="#2166ac", lw=0.9,
                zorder=4, label="Epithermal")
    axB.scatter(por.age, por.dist, s=22, marker="*", c="#b2182b", lw=0,
                zorder=5, label="Porphyry")
    axB.set_xlim(0, 66)
    axB.set_ylim(0, 900)
    axB.set_xlabel("Age (Ma)")
    axB.set_ylabel("Distance to its hotspot track (km)")
    axB.grid(True, color="#dddddd", lw=0.5, zorder=0)
    axB.set_axisbelow(True)
    axB.legend(loc="upper left", fontsize=6.4, ncol=2, framealpha=0.9,
               edgecolor="none", facecolor="white", borderpad=0.4,
               columnspacing=1.0, handlelength=1.8)
    axB.set_title("B", loc="left", fontweight="bold", fontsize=10)

    fig.tight_layout()
    for ext in ("pdf", "png"):
        fn = os.path.join(FIG_OUT, f"F4_deposit_distance.{ext}")
        fig.savefig(fn, dpi=300 if ext == "png" else None)
        print("written:", fn)

    # numbers for the text (cited in the Discussion)
    for name, df in (("porphyry", por), ("epithermal", epi), ("background", bg)):
        h68 = halfwidth(df, "h68")
        h95 = halfwidth(df, "h95")
        print(f"{name:11s} n={len(df):4d}  median dist {df.dist.median():5.0f} km"
              f" | within the 68% half-width: {np.mean(df.dist.values <= h68)*100:4.1f}%"
              f" | 95%: {np.mean(df.dist.values <= h95)*100:4.1f}%"
              f" | median 68% half-width {np.median(h68):3.0f} km")


def main_modelos():
    """All models of the set (paper_models.py), shown equally.
    (A) CDF of the distance of each class in each model (thin lines); shading = median over
        the models of the median 68% half-width; dotted = same for 95%.
    (B) distance of each deposit: median over the models (marker) and range between models
        (bar); 68% (fill) and 95% (line) half-widths of each chain, median over the models.
        The rows of deposit_tracks_*.csv come in the same order in all models."""
    import paper_models as MP
    clases = {"por": ("deposit_tracks_porphyries.csv", "Age MY"),
              "epi": ("deposit_tracks_epithermal.csv", "Age MY"),
              "bg": ("deposit_tracks_samples.csv", "Age (Ma)")}
    D = {k: [] for k in clases}
    H = {k: [] for k in clases}
    B68, B95 = [], []
    rows = []
    for c in MP.CARPETAS:
        P_ = os.path.join(ROOT, "models", c)
        bb = load_bands(os.path.join(P_, f"mc_collision_bands_total{TAG}.csv"))
        for k, (fn, acol) in clases.items():
            df = pd.read_csv(os.path.join(P_, fn), low_memory=False)
            if k == "bg":
                meth = df["Method"].astype(str).str.strip().str.lower()
                age = pd.to_numeric(df["Age (Ma)"], errors="coerce")
                bad = meth.isin(EXCLUDED) | df["Method"].isna() | (meth == "") | (meth == "nan")
                df = df[(~bad) | (age < 1)]
            df = prep(df, acol)
            if k == "bg":
                df = df[df.age <= 63.5]
            h68 = np.array([np.interp(a, bb[r]["age"], bb[r]["h68"]) for r, a in zip(df.ridge, df.age)])
            h95 = np.array([np.interp(a, bb[r]["age"], bb[r]["h95"]) for r, a in zip(df.ridge, df.age)])
            D[k].append(df.assign(modelo=c, h68=h68, h95=h95))
            rows.append(dict(modelo=c, clase=k, n=len(df), dist_mediana_km=df.dist.median(),
                             frac_dentro_68=np.mean(df.dist.values <= h68),
                             frac_dentro_95=np.mean(df.dist.values <= h95),
                             semiancho68_mediano_km=np.median(h68)))
        for r, b in bb.items():
            B68.append(pd.Series(b["h68"], index=b["age"], name=(c, r)))
            B95.append(pd.Series(b["h95"], index=b["age"], name=(c, r)))
    R = pd.DataFrame(rows)

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(7.3, 3.7))
    sty = {"por": ("Porphyry deposits", "#b2182b", "-", 1.0),
           "epi": ("Epithermal deposits", "#2166ac", (0, (5, 2)), 0.9),
           "bg": ("Igneous samples", "#555555", (0, (1, 1.4)), 0.9)}
    for k, (label, color, ls, lw) in sty.items():
        for df in D[k]:
            x, y = ecdf(df.dist)
            axA.plot(x, y, color=color, ls=ls, lw=lw, alpha=0.75, zorder=3)
        n = int(R[R.clase == k].n.median())
        axA.plot([], [], color=color, ls=ls, lw=1.6, label=f"{label} (n = {n})")
    rp = R[R.clase == "por"]
    h68m = float(np.median([np.median(df.h68) for df in D["por"]]))
    h95m = float(np.median([np.median(df.h95) for df in D["por"]]))
    axA.axvspan(0, h68m, color="#bdbdbd", alpha=0.45, zorder=0)
    axA.axvline(h95m, color="#777777", lw=0.9, ls=":", zorder=1)
    # vertical inside the band: on the left it overlaps the axis
    axA.text(h68m / 2, 0.97, "median 68% half-width", fontsize=7.0, rotation=90,
             ha="center", va="top", color="#444444")
    axA.text(h95m + 10, 0.985, "median 95%\nhalf-width", fontsize=7.2,
             ha="left", va="top", color="#666666")
    axA.set_xlim(0, 900); axA.set_ylim(0, 1)
    axA.set_xlabel("Restored distance to nearest hotspot track (km)")
    axA.set_ylabel("Cumulative fraction")
    axA.grid(True, color="#dddddd", lw=0.5, zorder=0); axA.set_axisbelow(True)
    h, l = axA.get_legend_handles_labels()
    h.append(plt.Line2D([0], [0], color="0.4", lw=0.8))
    l.append(f"one line per model ({len(MP.CARPETAS)} models)")
    axA.legend(h, l, loc="lower right", frameon=False, fontsize=7.4)
    axA.set_title("A", loc="left", fontweight="bold", fontsize=10)

    b68 = pd.concat(B68, axis=1); b95 = pd.concat(B95, axis=1)
    for ridge in [r for r in ("Easter", "TT", "JFR", "Foundation")
                  if r in b68.columns.get_level_values(1)]:
        m68 = b68.xs(ridge, axis=1, level=1).median(axis=1).dropna()
        m95 = b95.xs(ridge, axis=1, level=1).median(axis=1).dropna()
        col = RIDGE_COLOR[ridge]
        axB.fill_between(m68.index, 0, m68.values, color=col, alpha=0.16, lw=0, zorder=1)
        axB.plot(m95.index, m95.values, color=col, lw=1.0, ls=(0, (4, 2)), zorder=2,
                 label=f"{RIDGE_LABEL[ridge]} 95%")
        axB.plot(m68.index, m68.values, color=col, lw=1.2, zorder=2)
    for k, mk, s, col, lab, z in (("bg", "o", 4, "#bbbbbb", "Igneous samples", 3),
                                  ("epi", "x", 18, "#2166ac", "Epithermal", 4),
                                  ("por", "*", 22, "#b2182b", "Porphyry", 5)):
        M = np.column_stack([df.dist.values for df in D[k]])
        age = D[k][0].age.values
        med = np.median(M, axis=1)
        if k != "bg":
            axB.vlines(age, M.min(axis=1), M.max(axis=1), color=col, lw=0.5, alpha=0.6,
                       zorder=z - 0.5)
        axB.scatter(age, med, s=s, marker=mk, c=col, lw=0.9 if mk == "x" else 0, zorder=z,
                    label=lab + (" (median, range\nbetween models)" if k == "por" else ""))
    axB.set_xlim(0, 66); axB.set_ylim(0, 900)
    axB.set_xlabel("Age (Ma)")
    axB.set_ylabel("Distance to its hotspot track (km)")
    axB.grid(True, color="#dddddd", lw=0.5, zorder=0); axB.set_axisbelow(True)
    axB.legend(loc="upper left", fontsize=6.2, ncol=2, framealpha=0.9, edgecolor="none",
               facecolor="white", borderpad=0.4, columnspacing=1.0, handlelength=1.8)
    axB.set_title("B", loc="left", fontweight="bold", fontsize=10)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fn = os.path.join(FIG_OUT, f"F4_deposit_distance{MP.SUFIJO_SALIDA}.{ext}")
        fig.savefig(fn, dpi=300 if ext == "png" else None)
        print("written:", fn)
    R.to_csv(os.path.join(FIG_OUT, f"F4_deposit_distance{MP.SUFIJO_SALIDA}_numbers.csv"),
             index=False, float_format="%.3f")
    print(R.groupby("clase")[["dist_mediana_km", "frac_dentro_68", "frac_dentro_95",
                              "semiancho68_mediano_km"]].agg(["min", "median", "max"])
          .round(2).to_string())


if __name__ == "__main__":
    if TAG:
        main_modelos()
    else:
        main()
