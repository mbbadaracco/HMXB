#!/usr/bin/env python
"""Build the two lists of catalogue detections whose published fits we repeat.

Run with the project environment, from work/:

    env/bin/python select_stage1_targets.py

Chandra: every detection of one of our sources that carries a power-law fit
in CSC 2.1.1, identified by observation and source name. The catalogue fits
one spectrum per observation detection, so each row is a separate target.

XMM-Newton: every detection that holds a pipeline spectrum and belongs to one
of our stacked sources that carries a SPEC_ fit in 5XMM-DR15. The fit results
live on the stacked-source rows and the spectra on the detection rows, so the
two are joined here by SRCID.

Nothing is cross-matched: the Chandra list is keyed by observation and source
name, the XMM-Newton list by SRCID, all of them already in our tables.
"""

import os

import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.table import Table

base_path = "/home/marina/Doctorado/2026/HMXB_project"
csc_pull = os.path.join(base_path, "hmxb_5arcsec_csc.tsv")
xmm_catalogue = os.path.join(base_path, "5XMM-DR15", "5XMM_DR15.fits.gz")
xmm_ours = os.path.join(base_path, "work", "Output", "hmxb_5xmm_crossmatch.fits")
out_csc = os.path.join(base_path, "work", "Output", "csc_fitted_detections.tsv")
out_xmm = os.path.join(base_path, "work", "Output", "xmm_spectra_detections.tsv")


def read_cscview(path):
    """The CSCview TSV opens with one commented description per column, then
    the header line, then the data. The header line is repeated among the
    rows, so it is dropped."""
    lines = open(path).read().split("\n")
    start = next(i for i, l in enumerate(lines) if l and not l.startswith("#") and "\t" in l)
    names = [n.strip() for n in lines[start].split("\t")]
    d = pd.read_csv(path, sep="\t", skiprows=start, header=None, names=names)
    d["name"] = d["name"].astype(str).str.strip()
    return d[d["name"] != "name"].copy()


# Chandra
d = read_cscview(csc_pull)
d["obsid"] = pd.to_numeric(d["obsid"], errors="coerce").astype("Int64")
d["powlaw_nh"] = pd.to_numeric(d["powlaw_nh"], errors="coerce")
fitted = d[d["powlaw_nh"].notna()][["obsid", "name"]].drop_duplicates()
fitted = fitted.sort_values(["obsid", "name"])
fitted.to_csv(out_csc, sep="\t", index=False)
print("Chandra: %d detections, %d observations, %d sources -> %s"
      % (len(fitted), fitted["obsid"].nunique(), fitted["name"].nunique(), out_csc))

# XMM-Newton
ours = Table.read(xmm_ours)
nh = np.asarray(ours["SPEC_NH_PL"], dtype=float)
srcids = set(np.asarray(ours["SRCID"])[np.isfinite(nh)].tolist())
print("XMM-Newton: %d of our %d stacked sources carry a SPEC_ fit"
      % (len(srcids), len(ours)))

with fits.open(xmm_catalogue, memmap=True) as h:
    cat = h[1].data
    keep = np.isin(cat["SRCID"], list(srcids)) & cat["SPECTRA"]
    # FITS columns are big-endian; pandas sorts only native byte order
    rows = pd.DataFrame({
        "srcid": cat["SRCID"][keep].astype("int64"),
        "obsid": cat["OBS_ID"][keep].astype(str),
        "srcnum": cat["PPS_SRCNUM"][keep].astype("int64"),
        "iauname": cat["IAUNAME"][keep].astype(str),
    })
rows = rows.sort_values(["srcid", "obsid"])
rows.to_csv(out_xmm, sep="\t", index=False)
print("XMM-Newton: %d detections with a spectrum, %d observations -> %s"
      % (len(rows), rows["obsid"].nunique(), out_xmm))
