# -*- coding: utf-8 -*-
"""Builds the .rot files of the article models from the starting rotation model
(inputs/rotations/Global_Model_R_BG_2025.rot, with the sources of each row in its comment and
in its original timescale) and checks the result against the .rot files of the runs.

Order per model (circuit x frame):
  a. circuit on the starting .rot: Q22 as is; WD25 with models/WD25_WK08-A/build_rot_wd25.py
     (tables S2.2-S2.4 of Wilson and DeMets 2025 in models/WD25_WK08-A/wd25/);
  b. frame (plate 002): WK08-A as is; GMHRF with rotation_models/frames/gmhrf/build_rot_gmhrf.py
     (Table S3 of Doubrovine et al. 2012) and its CK95 -> GTS2012 ages (remap_gmhrf_gts2012.py,
     chrons from analysis/outputs/chrons_CK95_vs_GTS2012_full.csv); GWG20 with
     rotation_models/frames/gwg20/build_rot_gwg20.py (Table S6 of Gaastra et al. 2022);
  c. common steps 1-4 below;
  d. SG12: analysis/build_rot_sg12.py on the finished Q22 of the same frame (Somoza and
     Ghidella 2012, Table 1, written in the script).
GMHRF-A' and GMHRF-B share the .rot; only the prior of step 03b distinguishes them.

Common steps:
  1. West-East Antarctica (804 rel 802) from Granot et al. (2013, Table 1), inverted, between
     27.859 and 40.13 Ma (analysis/rot_sens_granot2013.py --base).
  2. Node 804 at 43.8 Ma from Cande et al. (2000) in place of the constant node at 45.03 Ma
     added by step 1.
  3. Notes on the Farallon -> Nazca change at 21.767 Ma (rows 902 and 911). The tag in these
     notes makes rot_gts2020.py treat row 902 FAR-PAC as GTS2012, like its partner NAZ-PAC
     (C6Bn.1 of Quiero et al. 2022).
  4. Ages to GTS2020 (analysis/rot_gts2020.py; it is not idempotent: it is applied only once,
     on the freshly built copy).
The hotspot rows (plates 13-16) are written later by step 03b; the covariances added to the
circuit are applied when the .rot is read (analysis/circuit_covariance.py).

Verification: the 12 .rot files match those of the runs outside plates 13-16 (rows, poles,
angles, fixed plates and covariances); in the three GWG20 models, 76 covariances of plate 002
differ in the fifth significant figure (those of the runs were rotated after rounding, with
analysis/gwg20_rotate_covariance.py).

Usage:  py build_rotation_models.py [output_folder]     (default outputs/rot/; compares with models/*/rot)
"""
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
AN = os.path.join(ROOT, "analysis")
BASE = os.path.join(ROOT, "inputs", "rotations", "Global_Model_R_BG_2025.rot")
ENC = "latin-1"                      # the starting .rot mixes encodings in its comments

NOTA_902_FAR = ("[this work: FAR=NZ until the breakup, propagated cov null; "
                "see row 911-902 at this age]")
NOTA_902_NAZ = ("[this work: copy of row 911 NAZ-PAC 6Bn.1 (same pole and cov); "
                "in the MC it shares its draw (SAME_DRAW in 04_mc_collisions.py)]")
NOTA_911 = (" NZ=FAR before the breakup. [this work] Propagated cov null: NZ-FAR = "
            "inv(FAR-PAC) x NZ-PAC, and the two rows at this age are the same rotation and the same "
            "measurement (6Bn.1), with correlation one. It is therefore not drawn in the MC; row 902 "
            "FAR-PAC at this age shares the draw of 911 NAZ-PAC (SAME_DRAW in 04_mc_collisions.py)")
CANDE_804 = "804 43.8  -18.15  -17.85    1.7  802 !CEED MBL-ANT Cande et al, Nature 2000"


def modulo(ruta, **globales):
    """Loads a project script as a module, sets its input/output paths and runs main()."""
    import importlib.util
    d = os.path.dirname(ruta)
    sys.path.insert(0, d)
    try:
        spec = importlib.util.spec_from_file_location(os.path.splitext(os.path.basename(ruta))[0] + "_rep", ruta)
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        for k, v in globales.items():
            setattr(m, k, v)
        import contextlib
        import io
        with contextlib.redirect_stdout(io.StringIO()):
            m.main()
    finally:
        sys.path.remove(d)


def py(*args):
    subprocess.run([sys.executable, *args], cwd=AN, check=True, stdout=subprocess.DEVNULL)


def editar(fn, f):
    lineas = open(fn, encoding=ENC).read().split("\n")
    open(fn, "w", encoding=ENC).write("\n".join(f(lineas)))


def una(lineas, patron):
    i = [k for k, l in enumerate(lineas) if re.match(patron, l)]
    if len(i) != 1:
        raise RuntimeError(f"{patron!r}: {len(i)} rows")
    return i[0]


def paso_cande(lineas):
    lineas[una(lineas, r"804 45\.03 ")] = CANDE_804
    return lineas


def paso_notas_21767(lineas):
    k = una(lineas, r"902 21\.767 +90\.0 +0\.0 +0\.0 +911 ")
    lineas[k] = lineas[k].rstrip() + " " + NOTA_902_FAR
    k = una(lineas, r"902 +21\.767 .* 901 +!NAZ-PAC ")
    lineas[k] = lineas[k].replace("!NAZ-PAC ", "!NAZ-PAC " + NOTA_902_NAZ + " ", 1)
    k = una(lineas, r"911 21\.767 +0 +90 +0 +902!")
    lineas[k] = lineas[k].rstrip() + NOTA_911
    return lineas


def comunes(out, notas=True):
    """Steps 1-4 on a .rot already built with its circuit and its frame."""
    py("rot_sens_granot2013.py", out, out, "--base")
    editar(out, paso_cande)
    if notas:                       # WD25 has no FAR -> NZ change at 21.767 Ma (it moves it to 33.705)
        editar(out, paso_notas_21767)
    py("rot_gts2020.py", out)
    shutil.rmtree(os.path.join(os.path.dirname(out), "backup_pre_gts2020"), ignore_errors=True)


def cargar(ruta):
    """Imports a project script as a module (without running its main)."""
    import importlib.util
    d = os.path.dirname(ruta)
    if d not in sys.path:
        sys.path.insert(0, d)
    spec = importlib.util.spec_from_file_location(os.path.splitext(os.path.basename(ruta))[0] + "_rep", ruta)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


PRE = os.path.join(ROOT, "rotation_models", "frames")


def circuito(nombre, out):
    """Circuit on the starting .rot, before the common steps."""
    if nombre == "Q22":
        shutil.copyfile(BASE, out)
    elif nombre == "WD25":
        # models/WD25_WK08-A/build_rot_wd25.py: plates 911 and 901 from Wilson and DeMets (2025), tables S2.2-S2.4
        # (models/WD25_WK08-A/wd25/), and Farallon/Nazca crossing moved to 33.705 Ma
        modulo(os.path.join(ROOT, "models", "WD25_WK08-A", "build_rot_wd25.py"),
               SRC=BASE, OUT=out, OUTDIR=os.path.dirname(out))


def marco_gmhrf(out):
    """Plate 002 <- inverse of Table S3 "Pacific" of Doubrovine et al. (2012) (build_rot_gmhrf.py,
    version A) and its CK95 -> GTS2012 ages (remap_gmhrf_gts2012.py); the common steps take it to
    GTS2020. The hotspot drift (GMHRF-B) is applied by step 03b with PAPER_HS_PRIOR=gmhrf."""
    g = cargar(os.path.join(PRE, "gmhrf", "build_rot_gmhrf.py"))
    g.BASE, g.OUT_A = out, out + ".tmp"
    g.write_version_a()
    r = cargar(os.path.join(PRE, "gmhrf", "remap_gmhrf_gts2012.py"))
    f, _ = r.remap_function()
    lineas = []
    with open(out + ".tmp", encoding="utf-8", errors="replace") as fh:
        for ln in fh:
            tok = ln.split()
            if tok and tok[0] in ("002", "2") and r.KEY in ln:
                i = ln.index(tok[1])
                ln = ln[:i] + f"{float(f(float(tok[1]))):.3f}" + ln[i + len(tok[1]):]
                ln = ln.replace(r.KEY, "Table S3 Pacific, age CK95->GTS2012 (", 1)
            lineas.append(ln)
    with open(out, "w", encoding="utf-8", newline="\n") as fh:
        fh.writelines(lineas)
    os.remove(out + ".tmp")


def marco_gwg20(out):
    """Plate 002 <- inverse of Table S6 of Gaastra et al. (2022), cov/kappa and C' = B C B^T
    (build_rot_gwg20.py, write_version)."""
    w = cargar(os.path.join(PRE, "gwg20", "build_rot_gwg20.py"))
    w.write_version(out, out + ".tmp", w.read_s6())
    os.replace(out + ".tmp", out)


def modelo(circ, marco, out, q22=None):
    """A complete .rot (without the hotspot rows of 03b). SG12 is built on the finished Q22 of the
    same frame (analysis/build_rot_sg12.py, Somoza and Ghidella 2012, Table 1)."""
    if circ == "SG12":
        py("build_rot_sg12.py", q22, out)
        return
    circuito(circ, out)
    if marco == "GMHRF":
        marco_gmhrf(out)
    elif marco == "GWG20":
        marco_gwg20(out)
    comunes(out, notas=circ != "WD25")


# (circuit, frame) -> (folder, .rot of the run)
MODELOS = {
    ("Q22", "WK08-A"): ("Q22_WK08-A", "Global_Model_R_BG_2025_hsb.rot"),
    ("Q22", "GMHRF"): ("Q22_GMHRF-B", "Global_Model_R_BG_2025_GMHRF_gts12_hsbgmhrf.rot"),
    ("Q22", "GWG20"): ("Q22_GWG20", "Global_Model_R_BG_2025_GWG20_hsb.rot"),
    ("WD25", "WK08-A"): ("WD25_WK08-A", "Global_Model_R_BG_2025_WD25_hsb.rot"),
    ("WD25", "GMHRF"): ("WD25_GMHRF-B", "Global_Model_R_BG_2025_WD25_GMHRF_gts12_hsbgmhrf.rot"),
    ("WD25", "GWG20"): ("WD25_GWG20", "Global_Model_R_BG_2025_WD25_GWG20_hsb.rot"),
    ("SG12", "WK08-A"): ("SG12_WK08-A", "Global_Model_R_BG_2025_SG12_hsb.rot"),
    ("SG12", "GMHRF"): ("SG12_GMHRF-B", "Global_Model_R_BG_2025_SG12_GMHRF_gts12_hsbgmhrf.rot"),
    ("SG12", "GWG20"): ("SG12_GWG20", "Global_Model_R_BG_2025_SG12_GWG20_hsb.rot"),
}
# GMHRF-A' (supplementary) uses the same .rot as GMHRF-B: only the prior of 03b changes
A_PRIMA = {"Q22": ("Q22_GMHRF-A", "Global_Model_R_BG_2025_GMHRF_gts12_hsb.rot"),
           "WD25": ("WD25_GMHRF-A", "Global_Model_R_BG_2025_WD25_GMHRF_gts12_hsb.rot"),
           "SG12": ("SG12_GMHRF-A", "Global_Model_R_BG_2025_SG12_GMHRF_gts12_hsb.rot")}


# ------------------------------------------------------------------ verification
HOTSPOTS = {13, 14, 15, 16}


def tabla(fn):
    eng = cargar(os.path.join(ROOT, "models", "Q22_WK08-A", "02_rot_engine.py"))
    return {p: r for p, r in eng.parse_rot(fn).items() if p not in HOTSPOTS}


def comparar(a, b, tol=1e-6):
    """Same rows (plate, age, pole, angle, fixed) and same covariances, without hotspots."""
    import numpy as np
    ta, tb = tabla(a), tabla(b)
    malas, redondeo = [], 0
    for p in sorted(set(ta) | set(tb)):
        ra, rb = ta.get(p, []), tb.get(p, [])
        if len(ra) != len(rb):
            malas.append((p, f"{len(ra)} rows vs {len(rb)}"))
            continue
        for x, y in zip(ra, rb):
            if x[4] != y[4] or not np.allclose(x[:4], y[:4], atol=tol, rtol=0):
                malas.append((p, x[:5], y[:5]))
            elif (x[5] is None) != (y[5] is None):
                malas.append((p, x[0], "covariance in only one"))
            elif x[5] is not None:
                d = np.abs(x[5] - y[5]).max() / np.abs(y[5]).max()
                if d > 1e-4:          # the covariances are written with 5 significant figures
                    malas.append((p, x[0], f"covariance, relative diff. {d:.1e}"))
                elif d > 0:
                    redondeo += 1
    n = sum(len(v) for v in ta.values())
    print(f"  {os.path.basename(b)}: {n} rows; different {len(malas)}"
          + (f"; {redondeo} covariances differ in the 5th significant figure" if redondeo else ""))
    for m in malas[:10]:
        print("    ", m)
    return not malas


def main():
    sal = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "outputs", "rot")
    os.makedirs(sal, exist_ok=True)
    ok, hechos = True, {}
    for circ in ("Q22", "WD25", "SG12"):
        for marco in ("WK08-A", "GMHRF", "GWG20"):
            carpeta, ref = MODELOS[(circ, marco)]
            out = os.path.join(sal, ref.replace("_hsbgmhrf", "").replace("_hsb", ""))
            modelo(circ, marco, out, q22=hechos.get(("Q22", marco)))
            hechos[(circ, marco)] = out
            print(f"{circ} + {marco} -> {os.path.basename(out)}")
            ok &= comparar(out, os.path.join(ROOT, "models", carpeta, "rot", ref))
            if marco == "GMHRF":
                c2, r2 = A_PRIMA[circ]
                ok &= comparar(out, os.path.join(ROOT, "models", c2, "rot", r2))
    print("identical (without hotspots)" if ok else "THERE ARE DIFFERENCES")


if __name__ == "__main__":
    main()
