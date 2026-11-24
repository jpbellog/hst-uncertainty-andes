# -*- coding: utf-8 -*-
"""Converts to GTS2020 (Ogg 2020) the ages of the nodes of the
hotspot -> South America chain in the .rot files of the eight models.

The model mixed timescales. The Neogene branches (Quiero et al. 2022, DeMets and
Merkouriev 2019, Granot and Dyment 2018, Wilson and DeMets 2025, Gaastra et al. 2022,
GMHRF and S&G 2012 already remapped) are in GTS2012; the Paleogene ones (Croon et al.
2008, Rowan and Rowley 2014, CEED: Cande et al. 2010, Cande and Stock 2004, Tikku and
Cande 2000, Bernard et al. 2005, Gaina et al. 1998; Mueller et al. 1999;
Wright et al. 2015; and the ages we assigned to Granot et al. 2013) are in CK95.

Rules per row (only plates of the chain; the rest of the .rot is not touched):
  - hotspots 13, 14, 15: not touched (03b rewrites them when refitting);
  - 002 of WK08-A: not touched (radiometric ages, no polarity timescale);
  - age >= 84 Ma: not touched;
  - top of C34n (age 83.0 or 83.5 labelled C34y/An34/A34y/34y, or 83.0/83.5 on
    a plate other than 002): 83.65 (GTS2020);
  - comment naming a source in GTS2012 (FUENTES_GTS12): GTS2012 -> GTS2020;
  - the rest: CK95 -> GTS2020.
Remapping functions: timescale_GTS2020.remap (piecewise linear through the boundaries of
the normal chrons). Poles, angles and covariances do not change.

A row without a recognizable source with the same age as a GTS2012 row of the
same plate (the partner row of a fixed-plate change, e.g. 802 28.278 709) follows
the rule of its partner. A CK95 row that, once remapped, falls within 0.05 Myr
of a GTS2012 row of the same plate and the same fixed plate is the same chron dated
twice (e.g. C10n.2n(o): Quiero 901@28.278 and Croon 901@28.75 -> 28.285) and is
removed; the row from the GTS2012 source is kept. If two different rows end up within
0.002 Myr of each other, the script aborts without writing.

Usage:  py rot_gts2020.py --dry <rot> ...    (report only, does not write)
        py rot_gts2020.py <rot> ...          (backs up the file in a subfolder of rot/ and rewrites it)
"""
import importlib.util
import os
import re
import shutil
import sys
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
_s = importlib.util.spec_from_file_location("esc", os.path.join(HERE, "timescale_GTS2020.py"))
esc = importlib.util.module_from_spec(_s)
_s.loader.exec_module(esc)
F = {"ck95": esc.remap("ck95"), "gts12": esc.remap("gts12")}

CADENA = {(2, 901), (201, 701), (709, 701), (801, 802), (802, 701), (802, 709),
          (804, 802), (833, 801), (901, 804), (901, 833), (902, 901), (902, 911),
          (908, 9022), (911, 901), (911, 902), (9022, 901), (90201, 901), (90201, 902),
          (90202, 901), (90202, 902), (902001, 908), (902001, 90201), (902002, 908),
          (902002, 90202), (902003, 902), (902003, 908), (201, 911), (901, 0)}
FUENTES_GTS12 = ("DeMets & Merkouriev", "REDBACK", "@REF ???", "[this work",
                 "Quiero", "Wilson", "WD25", "W&D", "Somoza", "SG12", "Gaastra",
                 "GWG", "Table S3 Pacific (", "CK95->GTS2012")
C34 = re.compile(r"C34y|An34|A34y|\b34y|\"34y")
EDAD = re.compile(r"^(\s*\d+\s+)(-?\d+(?:\.\d+)?)(\s.*)$", re.S)


def regla(pl, fx, age, com):
    if pl in (13, 14, 15) or age >= 84 or age <= 0 or "GTS2020" in com:
        return "igual"
    if pl == 2 and not any(k in com for k in FUENTES_GTS12):
        return "igual"                                   # WK08-A
    if age in (83.0, 83.5) and (C34.search(com) or pl != 2):
        return "c34"
    if any(k in com for k in FUENTES_GTS12):
        return "gts12"
    return "ck95"


def nueva_edad(r, age):
    if r == "igual":
        return age
    if r == "c34":
        return esc.C34_TECHO["gts20"]
    return round(float(F[r](age)), 3)


def fila(ln):
    """(plate, age, fixed plate, comment, rotation) or None. The fixed plate may be
    attached to the "!" of the comment (e.g. "911 21.767  0  90  0 902! NZ=FAR")."""
    pre = ln.split("!", 1)[0].split()
    if len(pre) >= 6 and pre[0].isdigit() and re.match(r"^-?\d+$", pre[5]):
        com = ln.split("!", 1)[1] if "!" in ln else ""
        return int(pre[0]), float(pre[1]), int(pre[5]), com, tuple(pre[2:5])
    return None


def procesar(fn, escribir):
    bak = os.path.join(os.path.dirname(fn), "backup_pre_gts2020", os.path.basename(fn))
    if os.path.exists(bak):
        # not idempotent: an already converted .rot would be remapped again
        raise SystemExit(f"{fn}: already converted ({bak} exists)")
    lines = open(fn, encoding="utf-8", errors="surrogateescape", newline="").readlines()
    out, informe = [], defaultdict(list)
    vistos = defaultdict(dict)

    def campos(ln):
        return fila(ln)
    # first pass: rules and GTS2012 ages per plate and per (plate, fixed plate)
    reg, g12_pl, g12_pf = {}, defaultdict(set), defaultdict(list)
    for k, ln in enumerate(lines):
        c = campos(ln)
        if c and (c[0], c[2]) in CADENA:
            reg[k] = regla(c[0], c[2], c[1], c[3])
            if reg[k] == "gts12":
                g12_pl[c[0]].add(c[1])
                g12_pf[(c[0], c[2])].append(nueva_edad("gts12", c[1]))
    for k, ln in enumerate(lines):
        c = campos(ln)
        if k in reg and reg[k] == "ck95" and c[1] in g12_pl[c[0]] and not c[3].strip():
            reg[k] = "gts12"                         # partner row of a fixed-plate change
    for k, ln in enumerate(lines):
        c = campos(ln)
        if c and k in reg:
            pl, age, fx, com, rot = c
            if True:
                r = reg[k]
                a2 = nueva_edad(r, age)
                if r == "ck95" and any(abs(a2 - b) < 0.05 for b in g12_pf[(pl, fx)]):
                    informe[(pl, fx, "BORRADA")].append((age, a2))
                    continue
                informe[(pl, fx, r)].append((age, a2))
                for a_prev, rot_prev in vistos[(pl, fx)].items():
                    if abs(a_prev - a2) < 0.002 and rot_prev != rot:
                        raise SystemExit(f"{fn}: {pl}/{fx} {age} -> {a2} collides with {a_prev}")
                vistos[(pl, fx)][a2] = rot
                if a2 != age:
                    m = EDAD.match(ln)
                    ln = m.group(1) + f"{a2:.3f}" + m.group(3)
        out.append(ln)
    for (pl, fx, r), v in sorted(informe.items()):
        ch = [f"{a:g}->{b:g}" for a, b in v if a != b]
        print(f"  {pl:>6}/{fx:<6} {r:6s} n={len(v):3d}  " + (" ".join(ch[:12]) + (" ..." if len(ch) > 12 else "")))
    # check that the age order within each (plate, fixed plate) is unchanged
    seq = defaultdict(list)
    for ln in out:
        c = fila(ln)
        if c:
            seq[c[0]].append(c[1])
    for pl in {k[0] for k in CADENA}:
        s = seq.get(pl, [])
        if any(b < a for a, b in zip(s, s[1:])):
            raise SystemExit(f"{fn}: plate {pl} ends up out of order: {s}")
    if escribir:
        os.makedirs(os.path.dirname(bak), exist_ok=True)
        shutil.copy2(fn, bak)
        open(fn, "w", encoding="utf-8", errors="surrogateescape", newline="").writelines(out)


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    dry = "--dry" in sys.argv
    for fn in [a for a in sys.argv[1:] if a != "--dry"]:
        print(("[dry] " if dry else "") + fn)
        procesar(fn, not dry)
