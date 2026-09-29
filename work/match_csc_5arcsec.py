#!/usr/bin/env python
"""Match the 5 arcsec CSCview pull, and write the observation list.

    env/bin/python match_csc_5arcsec.py

The criterion is the one of notebook 2: the Mahalanobis distance of Pineau
et al. (2011), Appendix A, between the Gaia position with its covariance and
the CSC position with its error ellipse rescaled from the 95% region by
sqrt(5.99), accepted at dM2 < 9.21.

Two sources F23 reports are kept although the statistic rejects them; each
carries a Notes value saying so and giving its dM2, because a sample that
follows two rules has to say which applies to each row.

Writes Output/hmxb_csc_match_5arcsec.csv, every candidate pair with its
distance, its acceptance and its note, and rewrites
Observations/Chandra/obsids_chandra.txt with the observations of the sources
that are kept.
"""


import os
import re

import numpy as np
import pandas as pd
from astropy.table import Table

import crossmatch

base = "/home/marina/Doctorado/2026/HMXB_project"
pull = os.path.join(base, "hmxb_5arcsec_csc.tsv")
previous = os.path.join(base, "work", "Output", "hmxb_csc_crossmatch.csv")
out = os.path.join(base, "work", "Output", "hmxb_csc_match_5arcsec.csv")
obsid_file = os.path.join(base, "Observations", "Chandra", "obsids_chandra.txt")

ACCEPT = 9.21   # chi-square at 99%, two degrees of freedom
CV_CSC = 5.99   # the CSC error ellipse is a 95% region, two d.o.f.


def read_cscview(path):
    lines = open(path).read().split("\n")
    start = next(i for i, l in enumerate(lines) if l and not l.startswith("#") and "\t" in l)
    names = [n.strip() for n in lines[start].split("\t")]
    d = pd.read_csv(path, sep="\t", skiprows=start, header=None, names=names)
    d["usrid"] = d["usrid"].astype(str).str.strip()
    d["name"] = d["name"].astype(str).str.strip()
    return d[d["usrid"] != "usrid"].copy()


def sexagesimal(ra, dec):
    """CSCview writes positions as 'hh mm ss.ss' and '+dd mm ss.s'."""
    h, m, s = [float(v) for v in str(ra).split()]
    ra_deg = 15.0 * (h + m / 60.0 + s / 3600.0)
    sign = -1.0 if str(dec).strip().startswith("-") else 1.0
    dd, dm, ds = [abs(float(v)) for v in str(dec).split()]
    return ra_deg, sign * (dd + dm / 60.0 + ds / 3600.0)


d = read_cscview(pull)
d[["ra_decimal", "dec_decimal"]] = [sexagesimal(a, b) for a, b in zip(d["ra"], d["dec"])]

# The 5 arcsec pull does not carry the error ellipse, which the statistic
# needs. It is taken from the accepted table for the sources already there
# and fetched from the CSC 2.1 TAP service for the rest, so that no ellipse
# is invented and every one comes from the catalogue itself.
ell = pd.read_csv(previous)[["name", "err_ellipse_r0", "err_ellipse_r1",
                             "err_ellipse_ang"]]
ell["name"] = ell["name"].astype(str).str.strip()
ell = ell.drop_duplicates("name")
absent = sorted(set(d["name"]) - set(ell["name"]))
if absent:
    import pyvo
    q = ("SELECT m.name, m.err_ellipse_r0, m.err_ellipse_r1, m.err_ellipse_ang "
         "FROM csc21.master_source m WHERE m.name IN (%s)"
         % ",".join("'%s'" % n for n in absent))
    fetched = pyvo.dal.TAPService(
        "http://cda.cfa.harvard.edu/csctap").search(q).to_table().to_pandas()
    fetched["name"] = fetched["name"].astype(str).str.strip()
    print("error ellipses fetched from CSC 2.1: %d of %d missing"
          % (len(fetched), len(absent)))
    ell = pd.concat([ell, fetched], ignore_index=True)
d = d.merge(ell, on="name", how="left")
for c in ("err_ellipse_r0", "err_ellipse_r1", "err_ellipse_ang"):
    d[c] = pd.to_numeric(d[c], errors="coerce")
assert d["err_ellipse_r0"].notna().all(), "an ellipse is still missing"

# The Gaia position and its covariance, as in notebook 2
sample = Table.read(os.path.join(base, "work", "Output",
                                 "hmxb_sample_positions.vot")).to_pandas()
sample["ID"] = sample["ID"].astype(str).str.strip()
by_id = sample.set_index("ID")


def gaia_cov(row):
    sa = row["e_RA_ICRS"] / 1000.0 / 3600.0   # mas -> deg
    sd = row["e_DE_ICRS"] / 1000.0 / 3600.0
    return sa ** 2, 0.0, 0.0, sd ** 2


pairs = d[["usrid", "name", "ra_decimal", "dec_decimal",
           "err_ellipse_r0", "err_ellipse_r1", "err_ellipse_ang"]].drop_duplicates(
    subset=["usrid", "name"])
rows = []
for _, p in pairs.iterrows():
    if p["usrid"] not in by_id.index or pd.isna(p["err_ellipse_r0"]):
        continue
    s = by_id.loc[p["usrid"]]
    X = crossmatch.matrix_rotation(p["err_ellipse_r0"], p["err_ellipse_r1"],
                                   CV_CSC, p["err_ellipse_ang"])
    dM2 = crossmatch.Mahalanobis_distance_squared(
        *gaia_cov(s), *X, s["RA_ICRS"], s["DE_ICRS"],
        p["ra_decimal"], p["dec_decimal"])
    rows.append({"ID": p["usrid"], "name": p["name"],
                 # haversine_distance already returns arcseconds
                 "sep_arcsec": crossmatch.haversine_distance(
                     s["RA_ICRS"], s["DE_ICRS"],
                     p["ra_decimal"], p["dec_decimal"]),
                 "err_ellipse_r0": p["err_ellipse_r0"],
                 "err_ellipse_r1": p["err_ellipse_r1"],
                 "dM2_xray": dM2})

m = pd.DataFrame(rows)
m["accepted"] = m["dM2_xray"] < ACCEPT

# Reported by F23, returned by the 5 arcsec cone, rejected by the statistic,
# kept at the user's decision of 27 September.
KEEP_ANYWAY = {("IGR J13186-6257", "2CXO J131825.0-625815"),
               ("2MASS J22535512+6243368", "2CXO J225355.0+624337")}
m["kept"] = m["accepted"] | pd.Series(
    [(i, n) in KEEP_ANYWAY for i, n in zip(m.ID, m.name)], index=m.index)
m["Notes"] = ["accepted by the cross-match" if r.accepted else
              ("kept by decision: reported by F23, rejected by the cross-match "
               "at dM2 = %.1f (threshold 9.21)" % r.dM2_xray) if r.kept else
              "rejected by the cross-match at dM2 = %.1f" % r.dM2_xray
              for r in m.itertuples()]
m = m.sort_values(["ID", "dM2_xray"])
m.to_csv(out, index=False)

print("5 arcsec pull: %d rows, %d distinct (HMXB, 2CXO) pairs" % (len(d), len(m)))
print("accepted at dM2 < %.2f : %d pairs over %d systems"
      % (ACCEPT, m.accepted.sum(), m.loc[m.accepted, "ID"].nunique()))
print("rejected                : %d" % (~m.accepted).sum())
print("-> %s" % out)


print("%d candidate pairs: %d accepted, %d kept by decision, %d rejected"
      % (len(m), m.accepted.sum(), (m.kept & ~m.accepted).sum(), (~m.kept).sum()))
for r in m[m.kept & ~m.accepted].itertuples():
    print("    kept: %-24s %-24s sep %.2f\"  dM2 %.1f"
          % (r.ID, r.name, r.sep_arcsec, r.dM2_xray))
print("-> %s" % out)

pairs = set(zip(m.loc[m.kept, "ID"], m.loc[m.kept, "name"]))
rows = d[[(u, n) in pairs for u, n in zip(d["usrid"], d["name"])]]
obsids = sorted(set(int(o) for o in rows["obsid"].dropna()))
previous = sorted(int(x) for x in open(obsid_file)) if os.path.exists(obsid_file) else []
with open(obsid_file, "w") as fh:
    fh.write("\n".join(str(o) for o in obsids) + "\n")
print("\n%d kept sources over %d systems, %d observations -> %s"
      % (len(pairs), m.loc[m.kept, "ID"].nunique(), len(obsids), obsid_file))
print("   was %d; added %s" % (len(previous), sorted(set(obsids) - set(previous))))
