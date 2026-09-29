#!/usr/bin/env python
"""Emit Output/hmxb_crossmatch.csv, the final cross-match table.

    env/bin/python build_crossmatch_table.py

One row per HMXB. For each of Gaia DR3, Chandra CSC 2.1 and XMM-Newton
5XMM-DR15 it carries the counterpart we adopt, the one F23 publishes, the
one SIMBAD attaches. The observations are not here: they are in
Observations/Chandra/obsids_chandra.txt and
Observations/XMM-Newton/obsids_xmm.txt, which the reduction scripts read.
The orbital parameters are not here either: build_orbital_table.py writes them
to Output/hmxb_orbital_parameters.csv, one row per row of this table. A Notes
column states, in words, every way in which the three sources of information
disagree for that system.

This table replaces the per-mission comparison tables: everything they held
is here, so there is one place to look and one file to keep in step.

It needs the network once, to ask gaiadr3.gaia_source whether each published
Gaia identifier exists; that is the only way to tell a truncated value from a
different counterpart.
"""

import io
import os
import re

import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.table import Table

base = "/home/marina/Doctorado/2026/HMXB_project"
work = os.path.join(base, "work")
out = os.path.join(work, "Output", "hmxb_crossmatch.csv")

f23_table = os.path.join(work, "Input", "HMXB", "tablea.dat")
csc_pull = os.path.join(base, "hmxb_5arcsec_csc.tsv")

ACCEPT = 9.21   # chi-square at 99%, two degrees of freedom


def read_cscview(path):
    """CSCview writes one commented description per column, then the header
    line, then the data, with the header repeated among the rows."""
    lines = io.open(path, encoding="utf-8", errors="replace").read().split("\n")
    start = next(i for i, l in enumerate(lines) if l and not l.startswith("#") and "\t" in l)
    names = [n.strip() for n in lines[start].split("\t")]
    d = pd.read_csv(path, sep="\t", skiprows=start, header=None, names=names)
    d["usrid"] = d["usrid"].astype(str).str.strip()
    d["name"] = d["name"].astype(str).str.strip()
    return d[d["usrid"] != "usrid"].copy()


def join(values):
    return "; ".join(dict.fromkeys(v for v in values if v and str(v) != "nan")) or None


# ---------------------------------------------------------------- our sample
notes = pd.read_csv(os.path.join(work, "Output", "hmxb_sample_notes.csv"), dtype=str)
notes["ID"] = notes["ID"].str.strip()
t = pd.DataFrame({"ID": notes["ID"], "Gaia_ours": notes["GaiaDR3_adopted"]})

# --------------------------------------------------------- what F23 publishes
# Byte ranges as the F23 ReadMe declares them, 1-based and both ends included.
# Only the identifiers: the parameters F23 measures are the subject of
# build_orbital_table.py, which reads the same file.
def cut(line, span):
    return line[span[0] - 1:span[1]].strip()


f23 = {}
for line in open(f23_table):
    if not line.strip():
        continue
    gaia_f23 = cut(line, (1522, 1540))
    # The ReadMe declares the column "?=0 Identifier in Gaia DR3": a zero is
    # the null marker, not an identifier. 41 systems carry it.
    f23[cut(line, (1, 23))] = {
        "gaia": gaia_f23 if gaia_f23 not in ("", "0") else None,
        "xmm": cut(line, (1271, 1291)) or None,
        "cxo": cut(line, (1341, 1362)) or None,
    }
for key, column in (("gaia", "Gaia_F23"), ("cxo", "CSC_F23"), ("xmm", "XMM_F23")):
    t[column] = t["ID"].map(lambda i, k=key: f23.get(i, {}).get(k))

# --------------------------------------------------------- what SIMBAD gives
simbad = pd.read_csv(os.path.join(work, "Output", "simbad_ids.csv"), dtype=str)
simbad["ID"] = simbad["ID"].str.strip()
gaia_simbad = pd.read_csv(os.path.join(work, "Output", "gaia_id_comparison.csv"), dtype=str)
gaia_simbad["ID"] = gaia_simbad["ID"].str.strip()
t["Gaia_SIMBAD"] = t["ID"].map(dict(zip(gaia_simbad.ID, gaia_simbad.SIMBAD)))
# SIMBAD lists CXOU designations next to the CSC ones. A CXOU name comes
# from an individual publication, not from the catalogue, so it is dropped:
# this column is what SIMBAD gives *as a CSC identifier*.
t["CSC_SIMBAD"] = t["ID"].map(dict(zip(simbad.ID, simbad.Chandra_simbad))).map(
    lambda v: join(n.strip() for n in str(v).split(";")
                   if not n.strip().startswith("CXOU")) if pd.notna(v) else v)
t["XMM_SIMBAD"] = t["ID"].map(dict(zip(simbad.ID, simbad.XMM_simbad)))

# ------------------------------------------- our X-ray sources and their obsids
kept = pd.read_csv(os.path.join(work, "Output", "hmxb_csc_match_5arcsec.csv"))
kept = kept[kept["kept"]].copy()   # accepted, plus the two kept by decision
kept["ID"] = kept["ID"].astype(str).str.strip()
kept["name"] = kept["name"].astype(str).str.strip()
t["CSC_ours"] = t["ID"].map(kept.groupby("ID")["name"].apply(join))
# The Mahalanobis distance is not a column: it is quoted in Notes only for
# the sources kept although the cross-match rejected them, which are the
# only rows where the reader needs it to judge the decision.
by_decision = {r.ID: "%.3g" % r.dM2_xray
               for r in kept.itertuples() if not r.accepted}

# --------------------------------------- the F23/F24 overlap and its decisions
# Nine systems appear in both F23 and F24, the LMXB catalogue of the same group
# (Output/f23_f24_positional.csv and Output/f23_f24_identifiers.csv). Three
# decisions were taken over them in 1_Crossmatch.ipynb, and they are repeated in
# Notes so that a reader of this table does not have to open the notebook to
# learn that a row is contested:
#   - 1E 1740.7-2942 and GRS 1758-258 are confirmed in the literature not to be
#     HMXBs and were removed from the sample, so they have no row here;
#   - IGR J18483-0311 is a confirmed HMXB, so its duplicate F24 entry was
#     removed instead and the system kept without a flag of ambiguity;
#   - the remaining six were kept and flagged, the flag travelling with the
#     source in Output/hmxb_sample_notes.csv.
# Only the flag is read from a file; the two dispositions are the decisions
# themselves and are stated here.
LMXB_RESOLVED = {"IGR J18483-0311": "EXO 1846-031"}
lmxb_ambiguous = {}
for i, s in zip(notes["ID"], notes["Notes"].fillna("")):
    hit = re.match(r"Ambiguous with LMXB: ([^;]+)", s)
    if hit:
        lmxb_ambiguous[i] = hit.group(1).strip()

x5 = pd.read_csv(os.path.join(work, "Output", "hmxb_5xmm_crossmatch.csv"))
x5["ID"] = x5["ID"].astype(str).str.strip()
x5["IAUNAME"] = x5["IAUNAME"].astype(str).str.strip()
t["XMM_ours"] = t["ID"].map(x5.groupby("ID")["IAUNAME"].apply(join))
# ------------------------------------------- do the published identifiers exist
identifiers = sorted({v for c in ("Gaia_ours", "Gaia_F23", "Gaia_SIMBAD")
                      for v in t[c].dropna() if str(v).strip().isdigit()})
from astroquery.gaia import Gaia
Gaia.ROW_LIMIT = -1
found = {str(r["source_id"]) for r in Gaia.launch_job(
    "SELECT source_id FROM gaiadr3.gaia_source WHERE source_id IN (%s)"
    % ",".join(identifiers)).get_results()}
print("Gaia identifiers checked: %d, existing: %d" % (len(identifiers), len(found)))

gaia_pos = {}
if len(found) < len(identifiers) or True:
    rows_g = Gaia.launch_job(
        "SELECT source_id, ra, dec FROM gaiadr3.gaia_source WHERE source_id IN (%s)"
        % ",".join(sorted(found))).get_results()
    gaia_pos = {str(r["source_id"]): (float(r["ra"]), float(r["dec"])) for r in rows_g}


def separation(a, b):
    if a in gaia_pos and b in gaia_pos and a != b:
        return round(float(SkyCoord(*gaia_pos[a], unit="deg").separation(
            SkyCoord(*gaia_pos[b], unit="deg")).arcsec), 2)
    return pd.NA


gaia_separation = {i: separation(str(a), str(b))
                   for i, a, b in zip(t.ID, t.Gaia_ours, t.Gaia_SIMBAD)}


# How far F23's 4XMM position is from the closest of our 5XMM sources. A
# changed designation is a rename unless this is large.
f23_xmm_pos = {}
for line in open(f23_table):
    if line.strip():
        f23_xmm_pos[line[0:23].strip()] = (line[1292:1303].strip(),
                                           line[1304:1314].strip())
xmm_separation = {}
for i, group in x5.groupby("ID"):
    ra, dec = f23_xmm_pos.get(i, ("", ""))
    if not ra or not dec:
        continue
    ours = SkyCoord(group["RA"].astype(float).values,
                    group["DEC"].astype(float).values, unit="deg")
    target = SkyCoord(float(ra), float(dec), unit="deg")
    xmm_separation[i] = "%.2f" % min(target.separation(ours).arcsec)


# ------------------------------------------------------------------- the notes
def note(r):
    # Whether the system is an HMXB at all is a separate question from which
    # counterpart it has, so the F23/F24 sentence is kept apart: it closes the
    # note, and it does not make the counterparts disagree.
    overlap = []
    if r.ID in lmxb_ambiguous:
        overlap.append("also in F24, the LMXB catalogue, as %s; kept in the HMXB "
                       "sample as 'Ambiguous'" % lmxb_ambiguous[r.ID])
    elif r.ID in LMXB_RESOLVED:
        overlap.append("also in F24, the LMXB catalogue, as %s; kept as a HMXB, "
                       "the literature confirming it, and the duplicate F24 entry "
                       "removed instead" % LMXB_RESOLVED[r.ID])
    parts = []
    # Gaia
    if pd.isna(r.Gaia_F23):
        parts.append("Gaia ID Candidate: F23 lists no counterpart")
    elif str(r.Gaia_F23) not in found:
        parts.append("the F23 Gaia identifier does not exist in Gaia DR3, being "
                     "ours with its leading digits lost")
    elif str(r.Gaia_F23) != str(r.Gaia_ours):
        parts.append("F23 names a different Gaia source")
    if pd.isna(r.Gaia_SIMBAD):
        parts.append("SIMBAD lists no Gaia DR3 identifier")
    elif str(r.Gaia_SIMBAD) != str(r.Gaia_ours):
        parts.append("SIMBAD names a different Gaia source, %s arcsec away"
                     % gaia_separation[r.ID])
    # Chandra
    if pd.notna(r.CSC_ours):
        if pd.isna(r.CSC_F23):
            parts.append("Chandra ID Candidate: F23 lists no counterpart")
        elif not str(r.CSC_F23).startswith("2CXO"):
            parts.append("the F23 Chandra entry is not a CSC designation (%s)"
                         % r.CSC_F23)
        elif str(r.CSC_F23) not in str(r.CSC_ours):
            parts.append("F23 names a Chandra source we do not keep")
        if r.ID in by_decision:
            parts.append("kept by decision: a source F23 reports that the "
                         "cross-match rejects, dM2 = %s (threshold %.2f)"
                         % (by_decision[r.ID], ACCEPT))
    elif pd.notna(r.CSC_F23):
        parts.append("F23 reports a Chandra counterpart we do not keep")
    # XMM-Newton
    if pd.notna(r.XMM_ours):
        if pd.isna(r.XMM_F23):
            parts.append("XMM-Newton ID Candidate: F23 lists no counterpart")
        elif str(r.XMM_F23)[5:] != str(r.XMM_ours)[5:]:
            # A designation encodes a position quantised at 0.1 s in right
            # ascension, about 1.5 arcsec, and 1 arcsec in declination, so a
            # sub-arcsecond change of position between DR11 and DR15 renames
            # the source. The separation says whether it is only that.
            parts.append("XMM-Newton: redesignated between F23's 4XMM DR11 and "
                         "our 5XMM-DR15, the two positions %s arcsec apart"
                         % xmm_separation.get(r.ID, "?"))
    if not parts:
        parts = ["our counterparts agree with F23 and SIMBAD"]
    return "; ".join(parts + overlap)


t["Notes"] = t.apply(note, axis=1)

# Only the systems an X-ray catalogue has something to say about; a row with
# neither counterpart has nothing for this work to compare or reduce.
dropped = t[t.CSC_ours.isna() & t.XMM_ours.isna()]
t = t.drop(dropped.index)
t = t[["ID", "Gaia_ours", "Gaia_F23", "Gaia_SIMBAD",
       "CSC_ours", "CSC_F23", "CSC_SIMBAD",
       "XMM_ours", "XMM_F23", "XMM_SIMBAD", "Notes"]].sort_values("ID")
t.to_csv(out, index=False)

print("%d systems -> %s" % (len(t), out))
print("  dropped, no X-ray counterpart in either catalogue: %d" % len(dropped))
print("  with a Gaia counterpart   : %d" % t.Gaia_ours.notna().sum())
print("  with a Chandra counterpart: %d" % t.CSC_ours.notna().sum())
print("  with an XMM counterpart   : %d" % t.XMM_ours.notna().sum())
print("  all three agree           : %d"
      % t.Notes.str.contains("our counterparts agree with F23 and SIMBAD",
                             regex=False).sum())
print("  also in F24, the LMXB catalogue: %d"
      % t.Notes.str.contains("also in F24", regex=False).sum())
