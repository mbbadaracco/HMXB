#!/usr/bin/env python
"""List the systems that can be drawn as a mosaic, most observations first.

    env/bin/python build_mosaic_list.py

Writes two tables, because the observations of this sample fall into two kinds
and only one of them can be ordered by orbital phase.

Output/mosaic_systems.txt
    The systems with at least one observation carrying a phase. These are the
    phase mosaics, ordered from phase 0 to ~0.99:
        ./ds9_phase_mosaic.sh <system>

Output/mosaic_systems_no_phase.txt
    The systems with observations and no phase on any of them. They are still
    worth looking at, several having many observations, but the only order
    available is time, so they are drawn with
        ./ds9_phase_mosaic.sh <system> --by-mjd
    which tiles every observation by the mid-point of its exposure and writes
    the MJD on each tile instead of the phase. The `why` column says which of
    the four reasons of X-09 and X-10 left the system without a phase.

Both are ordered by how many observations the mosaic will hold, most first,
ties broken by name so the files are stable between runs. The `system` column
is the directory name under Observations/Images/, which is exactly the argument
the script takes, so a line can be copied straight into the command.

Counts are read from the per-system tables written by build_image_tables.py:
a row carries a phase when it has four fields, and does not when it has three.
Tper_from is read from Output/hmxb_orbital_parameters.csv and says how the
epoch reached the common convention, "as published" meaning it was already a
periastron passage and everything else naming the conversion of X-09.
"""

import glob
import os
import re

import pandas as pd

base = "/home/marina/Doctorado/2026/HMXB_project"
work = os.path.join(base, "work")
images = os.path.join(base, "Observations", "Images")
out_phase = os.path.join(work, "Output", "mosaic_systems.txt")
out_nophase = os.path.join(work, "Output", "mosaic_systems_no_phase.txt")

orb = pd.read_csv(os.path.join(work, "Output", "hmxb_orbital_parameters.csv"))
orb["ID"] = orb["ID"].astype(str).str.strip()
tper_from = dict(zip(orb["ID"], orb["Tper_from"].fillna("")))
porb = dict(zip(orb["ID"], orb["Porb_d"]))
mjd_t0 = dict(zip(orb["ID"], orb["MJD_T0"]))
t0_kind = dict(zip(orb["ID"], orb["T0_kind"].fillna("")))
mjd_tper = dict(zip(orb["ID"], orb["MJD_Tper"]))


def period(sid):
    """%.10g, not str(): the period comes out of the CSV as a float and
    4U 1907+097's 8.3753 prints as 8.375299999999998, which is binary noise
    and not a digit of the value."""
    p = porb.get(sid)
    return "%.10g" % p if pd.notna(p) else ""


# The epoch kinds that bear a fixed relation to the orbit, and so can be carried
# to a periastron passage once the longitude of periastron is known (X-09).
CONVERTIBLE = ("Tpi/2", "mid-eclipse", "superior conjunction", "inferior conjunction")


def why_no_phase(sid):
    """Which step of X-09 and X-10 this system stops at."""
    if pd.isna(porb.get(sid)):
        return "no orbital period"
    if pd.isna(mjd_t0.get(sid)):
        return "period, no epoch published"
    if pd.isna(mjd_tper.get(sid)):
        kind = t0_kind.get(sid) or "epoch of unknown kind"
        # Of the two ways a conversion fails, which one it is matters: a missing
        # omega can be filled by one paper, an outburst maximum never can.
        if kind in CONVERTIBLE:
            return "%s epoch, no omega published" % kind
        return "%s epoch, no fixed relation to the orbit" % kind
    return "phase lost to the period error"


phased, unphased = [], []
for table in sorted(glob.glob(os.path.join(images, "*", "observations.txt"))):
    folder = os.path.basename(os.path.dirname(table))
    lines = open(table).read().split("\n")
    sid = re.sub(r"^#\s*", "", lines[0]).strip()
    n_phase = n_total = n_chandra = n_xmm = 0
    for line in lines:
        f = line.split()
        if not f or line.startswith("#") or f[0] == "ChandraObsID":
            continue
        n_total += 1
        if len(f) == 4:
            n_phase += 1
        if f[0] != "-":
            n_chandra += 1
        else:
            n_xmm += 1
    if not n_total:
        continue
    if n_phase:
        phased.append((n_phase, n_total, n_chandra, n_xmm, folder, sid))
    else:
        unphased.append((n_total, n_chandra, n_xmm, folder, sid))

phased.sort(key=lambda r: (-r[0], r[4]))
unphased.sort(key=lambda r: (-r[0], r[3]))

# ------------------------------------------------------------------ phase order
head = []
head.append("# Systems that can be drawn as a phase mosaic, most observations first.")
head.append("# Built by work/build_mosaic_list.py from the tables of Observations/Images/.")
head.append("# The systems with no phase at all are in mosaic_systems_no_phase.txt.")
head.append("#")
head.append("# n_phase   observations that carry an orbital phase; these are the tiles of the mosaic")
head.append("# n_obs     observations in total, Chandra and XMM-Newton together")
head.append("# Chandra   of those, Chandra; XMM, XMM-Newton. A row without a phase is skipped by the script.")
head.append("# system    the directory under Observations/Images/, and the argument of ds9_phase_mosaic.sh:")
head.append("#               ./ds9_phase_mosaic.sh <system>")
head.append("# ID        the identifier of the working sample, which carries the spaces the directory cannot")
head.append("# Porb      orbital period in days")
head.append("# Tper_from how the epoch reached the periastron convention: \"as published\" means it was")
head.append("#               already one, anything else names the conversion of X-09")
head.append("#")
head.append("%-8s %-6s %-8s %-4s %-26s %-26s %-14s %s"
            % ("n_phase", "n_obs", "Chandra", "XMM", "system", "ID", "Porb_d", "Tper_from"))

body = []
for n_phase, n_total, n_chandra, n_xmm, folder, sid in phased:
    body.append("%-8d %-6d %-8d %-4d %-26s %-26s %-14s %s"
                % (n_phase, n_total, n_chandra, n_xmm, folder, sid,
                   period(sid), tper_from.get(sid, "")))

with open(out_phase, "w") as f:
    f.write("\n".join(head + body) + "\n")

# ------------------------------------------------------------------- time order
head = []
head.append("# Systems with observations and no orbital phase on any of them, most observations first.")
head.append("# Built by work/build_mosaic_list.py from the tables of Observations/Images/.")
head.append("# The systems that do have a phase are in mosaic_systems.txt.")
head.append("#")
head.append("# These cannot be ordered by phase, so they are drawn in time order, each tile")
head.append("# labelled with its MJD instead:")
head.append("#               ./ds9_phase_mosaic.sh <system> --by-mjd")
head.append("# The script falls back to that on its own for these systems, so the flag only")
head.append("# makes the intention explicit.")
head.append("#")
head.append("# n_obs     observations in total, and the number of tiles the mosaic will hold")
head.append("# Chandra   of those, Chandra; XMM, XMM-Newton")
head.append("# system    the directory under Observations/Images/, and the argument of the script")
head.append("# ID        the identifier of the working sample, which carries the spaces the directory cannot")
head.append("# Porb      orbital period in days, empty where none is published")
head.append("# why       where this system stops on the way to a phase, in the terms of X-09 and X-10")
head.append("#")
head.append("%-6s %-8s %-4s %-26s %-26s %-14s %s"
            % ("n_obs", "Chandra", "XMM", "system", "ID", "Porb_d", "why"))

body = []
for n_total, n_chandra, n_xmm, folder, sid in unphased:
    body.append("%-6d %-8d %-4d %-26s %-26s %-14s %s"
                % (n_total, n_chandra, n_xmm, folder, sid, period(sid), why_no_phase(sid)))

with open(out_nophase, "w") as f:
    f.write("\n".join(head + body) + "\n")

print("%d systems with a phase      -> %s" % (len(phased), out_phase))
print("   observations with a phase: %d of %d"
      % (sum(r[0] for r in phased), sum(r[1] for r in phased)))
print("%d systems without one       -> %s" % (len(unphased), out_nophase))
print("   observations there       : %d" % sum(r[0] for r in unphased))
print("   of those, %d systems have more than one observation"
      % sum(1 for r in unphased if r[0] > 1))
