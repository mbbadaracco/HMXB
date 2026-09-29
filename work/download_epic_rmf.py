#!/usr/bin/env python
"""Download the canned EPIC response matrices our XMM-Newton spectra name.

Run with the project environment, from work/:

    env/bin/python download_epic_rmf.py

The XMM-Newton pipeline distributes the source spectrum, its background and
the ARF, but no RMF: each spectrum names a canned matrix in its RESPFILE
keyword. This reads that keyword from the spectra of the detections listed in
Output/xmm_spectra_detections.tsv and fetches each matrix once.

MOS matrices are taken from the 5eV directory, which is the binning of the
pipeline spectra (SPECDELT = 5, DETCHANS = 2400). The pn matrices in the
repository carry a version suffix that the spectra do not name, so the
current one is taken and saved under the name the spectrum expects.
"""

import csv
import glob
import os
import urllib.request

from astropy.io import fits

base_path = "/home/marina/Doctorado/2026/HMXB_project/work"
detections = os.path.join(base_path, "Output", "xmm_spectra_detections.tsv")
pps = os.path.join(base_path, "xmm_pps")
root = os.path.join(base_path, "epic_rmf")
repo = "http://xmm-tools.cosmos.esa.int/external/xmm_ccf/ccf/extras/responses"
pn_version = "v22.0"  # the version the repository offers, not the one the pipeline used

# Which spectra belong to our detections
wanted = {}
with open(detections) as fh:
    for r in csv.DictReader(fh, delimiter="\t"):
        wanted.setdefault(r["obsid"], set()).add(int(r["srcnum"]))

responses = set()
spectra = 0
for obsid in sorted(wanted):
    for srcnum in sorted(wanted[obsid]):
        tag = "%03X" % srcnum  # the source number is hexadecimal in the file name
        for f in glob.glob(os.path.join(pps, obsid, "pps", "*SRSPEC?%s.FTZ" % tag)):
            if os.path.basename(f)[11:13] not in ("PN", "M1", "M2"):
                continue  # EPIC only, the RGS spectra of the same source are not used
            spectra += 1
            header = fits.getheader(f, 1)
            if "RESPFILE" in header:
                responses.add(header["RESPFILE"].strip())

print("%d EPIC spectra, %d distinct response matrices" % (spectra, len(responses)))
os.makedirs(root, exist_ok=True)

for name in sorted(responses):
    out = os.path.join(root, name)
    if os.path.exists(out):
        print("Skipping %s, already downloaded" % name)
        continue
    if name.startswith("epn"):
        url = "%s/PN/%s_%s.rmf" % (repo, name[:-4], pn_version)
    else:
        url = "%s/MOS/5eV/%s" % (repo, name)
    print("Downloading %s from %s" % (name, url))
    try:
        urllib.request.urlretrieve(url, out)
    except Exception as exc:
        print("   FAILED: %s" % exc)
        if os.path.exists(out):
            os.remove(out)
