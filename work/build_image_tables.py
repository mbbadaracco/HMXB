#!/usr/bin/env python
"""Lay out Observations/Images/ and write one observation table per system.

    env/bin/python build_image_tables.py

For every HMXB of the working sample that has at least one observation on
disk, this creates

    Observations/Images/<system>/
        Chandra_images/          (only if the system has Chandra data)
        XMM-Newton_images/       (only if it has XMM-Newton data)
        observations.txt

The directory name is the system identifier with runs of whitespace replaced
by a single underscore; the identifier itself is written in the header of the
table, so nothing is lost.

`observations.txt` has one row per observation, with the columns the analysis
needs: the Chandra observation identifier, the XMM-Newton one, the MJD and the
orbital phase. An observation belongs to one mission, so exactly one of the two
identifier columns is filled on each row.

MJD is the mid-point of the observation, (TSTART + TSTOP) / 2 expressed through
MJDREF, read from the reprocessed event files themselves rather than from a
catalogue. For XMM-Newton an observation has several exposures and cameras, so
the span is taken over all of its imaging event files.

The orbital phase is counted from the periastron passage, phase 0 = periastron,
using MJD_Tper and Porb_d of Output/hmxb_orbital_parameters.csv, that is the
epochs brought to the common convention in X-09. Where a system has no such
epoch the column is left empty: a phase measured from an outburst maximum or
from an arbitrary folding origin is not the same quantity and is not written
here as though it were.

It is also left empty where the phase has no information left in it. Carrying
an epoch across the n cycles that separate it from an observation costs
n * sigma_Porb, and that term dominates: at 151 cycles a period known to 4.5
per cent, as SAX J0635.2+0533's is, gives an uncertainty of nearly seven whole
cycles. The phase is written when its uncertainty is below 0.5 cycles and left
empty above, since at half a cycle the value is consistent with any phase.
"""

import glob
import os
import re

import numpy as np
import pandas as pd
from astropy.io import fits

base = "/home/marina/Doctorado/2026/HMXB_project"
work = os.path.join(base, "work")
images = os.path.join(base, "Observations", "Images")

orbital = os.path.join(work, "Output", "hmxb_orbital_parameters.csv")
csc_kept = os.path.join(work, "Output", "hmxb_csc_match_5arcsec.csv")
csc_pull = os.path.join(base, "hmxb_5arcsec_csc.tsv")
xmm_ours = os.path.join(work, "Output", "hmxb_5xmm_crossmatch.csv")
xmm_cat = os.path.join(base, "5XMM-DR15", "5XMM_DR15.fits.gz")

# Above this the phase is consistent with any value and is not written.
LOST = 0.5

CHANDRA = os.path.join(base, "Observations", "Chandra_repro")
XMM = os.path.join(base, "Observations", "XMM-Newton_repro")


def read_cscview(path):
    """CSCview writes one commented description per column, then the header
    line, then the data, with the header repeated among the rows."""
    lines = open(path, encoding="utf-8", errors="replace").read().split("\n")
    start = next(i for i, l in enumerate(lines)
                 if l and not l.startswith("#") and "\t" in l)
    names = [n.strip() for n in lines[start].split("\t")]
    d = pd.read_csv(path, sep="\t", skiprows=start, header=None, names=names)
    d["name"] = d["name"].astype(str).str.strip()
    return d[d["name"] != "name"].copy()


def folder(sid):
    return re.sub(r"\s+", "_", sid.strip())


def span(paths):
    """(MJD_start, MJD_stop) over a set of event files, or None."""
    lo, hi, ref = np.inf, -np.inf, None
    for p in paths:
        try:
            h = fits.getheader(p, 1)
            t0, t1, ref = float(h["TSTART"]), float(h["TSTOP"]), float(h["MJDREF"])
        except Exception:
            continue
        lo, hi = min(lo, t0), max(hi, t1)
    if ref is None or not np.isfinite(lo):
        return None
    return ref + lo / 86400.0, ref + hi / 86400.0


# --------------------------------------------------- which observation is whose
kept = pd.read_csv(csc_kept)
kept = kept[kept["kept"]].copy()
kept["ID"] = kept["ID"].astype(str).str.strip()
kept["name"] = kept["name"].astype(str).str.strip()

pull = read_cscview(csc_pull)
pull["obsid"] = pd.to_numeric(pull["obsid"], errors="coerce").astype("Int64")
by_name = pull.dropna(subset=["obsid"]).groupby("name")["obsid"].apply(
    lambda s: sorted({int(v) for v in s}))

chandra_of = {}
for sid, group in kept.groupby("ID"):
    obs = sorted({o for n in group["name"] for o in by_name.get(n, [])})
    if obs:
        chandra_of[sid] = obs

x5 = pd.read_csv(xmm_ours)
x5["ID"] = x5["ID"].astype(str).str.strip()
srcid_of = x5.groupby("ID")["SRCID"].apply(lambda s: {int(v) for v in s}).to_dict()
wanted = {v for s in srcid_of.values() for v in s}
with fits.open(xmm_cat, memmap=True) as h:
    cat = h[1].data
    keep = np.isin(cat["SRCID"], list(wanted))
    det = pd.DataFrame({"srcid": cat["SRCID"][keep].astype("int64"),
                        "obsid": cat["OBS_ID"][keep].astype(str)})
xmm_of = {}
for sid, srcids in srcid_of.items():
    obs = sorted({o for o in det[det.srcid.isin(srcids)]["obsid"]})
    if obs:
        xmm_of[sid] = obs

# ------------------------------------------------------- the MJD of every obsid
chandra_mjd, xmm_mjd = {}, {}
for obsid in {o for v in chandra_of.values() for o in v}:
    s = span(glob.glob(os.path.join(CHANDRA, str(obsid), "*repro_evt2.fits")))
    if s:
        chandra_mjd[obsid] = s
for obsid in {o for v in xmm_of.values() for o in v}:
    s = span(glob.glob(os.path.join(XMM, obsid, "*ImagingEvts.ds")))
    if s:
        xmm_mjd[obsid] = s

# ------------------------------------------------------------- the orbital phase
orb = pd.read_csv(orbital)
orb["ID"] = orb["ID"].astype(str).str.strip()
eph = {r.ID: (r.Porb_d, r.MJD_Tper, r.e_MJD_Tper, r.e_Porb_d, r.r_T0)
       for r in orb.itertuples() if pd.notna(r.MJD_Tper) and pd.notna(r.Porb_d)}

os.makedirs(images, exist_ok=True)
rows_written = phases_written = 0
systems = sorted(set(chandra_of) | set(xmm_of))
summary = []
for sid in systems:
    rows = []
    for obsid in chandra_of.get(sid, []):
        if obsid in chandra_mjd:
            lo, hi = chandra_mjd[obsid]
            rows.append((str(obsid), "", (lo + hi) / 2))
    for obsid in xmm_of.get(sid, []):
        if obsid in xmm_mjd:
            lo, hi = xmm_mjd[obsid]
            rows.append(("", obsid, (lo + hi) / 2))
    if not rows:
        continue
    rows.sort(key=lambda r: r[2])

    root = os.path.join(images, folder(sid))
    os.makedirs(root, exist_ok=True)
    if sid in chandra_of:
        os.makedirs(os.path.join(root, "Chandra_images"), exist_ok=True)
    if sid in xmm_of:
        os.makedirs(os.path.join(root, "XMM-Newton_images"), exist_ok=True)

    def sigma_phase(mjd):
        """Uncertainty of the phase at this MJD, cycle count included."""
        P, tper, e_tper, e_P = eph[sid][:4]
        n = abs(mjd - tper) / P
        return np.hypot(e_tper if pd.notna(e_tper) else 0,
                        n * (e_P if pd.notna(e_P) else 0)) / P, n

    head = ["# %s" % sid]
    if sid in eph:
        P, tper, e_tper, e_P, ref = eph[sid]
        sig, n = sigma_phase(max(rows, key=lambda r: abs(r[2] - tper))[2])
        head.append("# orbital phase counted from periastron, phase 0 = periastron:")
        head.append("#   phase = ((MJD - %.5f) / %.6f) mod 1" % (tper, P))
        head.append("#   Porb = %.6f d, T_per = MJD %.5f, epoch from %s" % (P, tper, ref))
        head.append("#   phase uncertainty at the furthest observation: %.3f cycles"
                    " (%.0f cycles from the epoch)" % (sig, n))
        if sig >= LOST:
            head.append("#   left empty: at this distance from the epoch the period"
                        " error, %s d," % e_P)
            head.append("#   has spread the phase over more than half a cycle, so it"
                        " carries no information.")
        elif sig >= 0.1:
            head.append("#   this is a large fraction of a cycle; read the phase with"
                        " that in mind.")
        summary.append((sid, len(rows), sig, n))
    else:
        head.append("# orbital phase: not computed. This system has no epoch on the")
        head.append("#   common convention, so a phase would not be comparable with"
                    " the others.")
    head.append("#")
    head.append("%-14s %-18s %-14s %s"
                % ("ChandraObsID", "XMM-NewtonObsID", "MJD", "orbital_phase"))

    body = []
    for cid, xid, mjd in rows:
        phase = ""
        if sid in eph and sigma_phase(mjd)[0] < LOST:
            P, tper = eph[sid][0], eph[sid][1]
            phase = "%.4f" % (((mjd - tper) / P) % 1.0)
            phases_written += 1
        body.append("%-14s %-18s %-14.5f %s" % (cid or "-", xid or "-", mjd, phase))
        rows_written += 1
    with open(os.path.join(root, "observations.txt"), "w") as f:
        f.write("\n".join(head + body) + "\n")

print("%d systems -> %s" % (len(systems), images))
print("  with Chandra data   : %d" % len(chandra_of))
print("  with XMM-Newton data: %d" % len(xmm_of))
print("  observation rows    : %d" % rows_written)
print("  rows carrying a phase: %d" % phases_written)
print()
print("phase uncertainty, cycle count included, worst observation per system:")
for sid, n_rows, sig, n in sorted(summary, key=lambda s: -s[2]):
    print("  %-20s %3d obs  %6.0f cycles  sigma_phase = %.3f" % (sid, n_rows, n, sig))
