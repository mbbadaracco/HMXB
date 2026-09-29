#!/usr/bin/env python
"""Download the CSC 2.1 per-detection spectral products of our sources.

Run with CIAO's python, from work/:

    ciaoinit
    python download_csc_products.py

For every (obsid, source name) pair that carries a power-law fit in the
catalogue, it retrieves the source spectrum, its background, the ARF, the
RMF and the region definition, so that the published fit can be repeated
on the catalogue's own data. It downloads nothing that is already there.
"""

import os
import sys

import ciao_contrib.cda.csccli as csc

base_path = "/home/marina/Doctorado/2026/HMXB_project"
target_file = os.path.join(base_path, "work", "Output", "csc_fitted_detections.tsv")
root = os.path.join(base_path, "work", "csc_products")
filetypes = "pha,arf,rmf,reg"
bands = "broad"
catalog = "csc2.1"

# (obsid, source name) of every detection with a published power-law fit
targets = {}
with open(target_file) as fh:
    next(fh)
    for line in fh:
        obsid, name = line.rstrip("\n").split("\t")
        targets.setdefault(obsid, set()).add(name)

print("%d observations, %d detections" % (len(targets), sum(len(v) for v in targets.values())))
os.makedirs(root, exist_ok=True)

for i, obsid in enumerate(sorted(targets, key=int), 1):
    wanted = targets[obsid]
    print("\n[%d/%d] obsid %s: %s" % (i, len(targets), obsid, ", ".join(sorted(wanted))))
    sys.stdout.flush()
    try:
        page = csc.search_src_by_obsid(obsid, csc.required_cols, catalog)
        found = csc.parse_csc_result(page)
    except Exception as exc:
        print("   query failed: %s" % exc)
        continue
    mine = [s for s in found if s["name"].strip() in wanted]
    if not mine:
        print("   none of the sources came back from the query")
        continue
    for src in mine:
        print("   %s  region_id %s  %s" % (src["name"], src["region_id"], src["instrument"]))
        sys.stdout.flush()
        try:
            csc.retrieve_files_per_src(src, root, filetypes, bands, catalog)
        except Exception as exc:
            print("   retrieval failed: %s" % exc)
