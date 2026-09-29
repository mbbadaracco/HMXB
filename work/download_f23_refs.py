#!/usr/bin/env python
"""Fetch the papers F23 cites for the orbital periods of our systems.

    env/bin/python download_f23_refs.py

They are read for two questions that F23's own table cannot answer: where
the zero point of each orbital period sits, and which systems have shown a
Type II outburst. Each bibcode is resolved to its arXiv eprint through the
ADS link gateway and the PDF is saved under f23_refs/ as <bibcode>.pdf.

They are not in Papers/ and so cannot be cited; they are read only.
"""

import os
import re
import subprocess
import time

import pandas as pd

base = "/home/marina/Doctorado/2026/HMXB_project"
table = os.path.join(base, "work", "Output", "hmxb_orbital_parameters.csv")
root = os.path.join(base, "work", "f23_refs")

os.makedirs(root, exist_ok=True)
d = pd.read_csv(table, dtype=str)
bibcodes = sorted({b.strip() for b in d["r_Porb"].dropna() if b.strip()})
print("%d distinct references for the orbital periods" % len(bibcodes))

# References already on disk under a name of their own, as the published
# article rather than the arXiv version. Fetching them again would leave two
# copies of the same paper.
ALREADY = {"2013MNRAS.434.2182G": "Goossens2013.pdf",
           "2024MNRAS.528..863B": "Bozzo2024.pdf"}

resolved = failed = 0
for bib in bibcodes:
    out = os.path.join(root, ALREADY.get(bib, "%s.pdf" % bib))
    if os.path.exists(out) and os.path.getsize(out) > 0:
        continue
    # The gateway answers with a redirect to arxiv.org/abs/<id> when an
    # eprint exists, and with nothing when it does not.
    location = subprocess.run(
        ["curl", "-s", "-o", "/dev/null", "-w", "%{redirect_url}", "--max-time", "25",
         "https://ui.adsabs.harvard.edu/link_gateway/%s/EPRINT_HTML" % bib],
        capture_output=True, text=True).stdout
    match = re.search(r"arxiv\.org/abs/(\S+)", location)
    if not match:
        print("   no eprint: %s" % bib)
        failed += 1
        continue
    subprocess.run(["curl", "-sL", "--max-time", "60", "-A",
                    "HMXB-project/1.0 (academic use)", "-o", out,
                    "https://arxiv.org/pdf/%s" % match.group(1)], check=False)
    if os.path.exists(out) and open(out, "rb").read(4) == b"%PDF":
        resolved += 1
    else:
        if os.path.exists(out):
            os.remove(out)
        print("   download failed: %s" % bib)
        failed += 1
    time.sleep(1)

have = len([f for f in os.listdir(root) if f.endswith(".pdf")])
print("fetched now: %d | without an eprint or failed: %d | on disk: %d"
      % (resolved, failed, have))
