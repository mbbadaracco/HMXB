#!/usr/bin/env python
"""Stage 1 for Chandra: repeat the CSC 2.1.1 power-law fits on their own data.

Run with CIAO's python, from work/:

    ciaoinit
    python refit_csc_stage1.py

For every detection of Output/csc_fitted_detections.tsv the catalogue's own
spectrum, background, ARF and RMF are loaded and fitted with the recipe of
Evans et al. (2024), Sect. 3.13: the PI spectrum grouped to a minimum of 16
counts per bin, the background subtracted, an absorbed power law with the
Balucinska-Church & McCammon cross-sections and the Anders & Grevesse
abundances, and chi-square with data variance as the statistic. Nothing is
extracted and nothing is re-reduced here: the data and the responses are the
catalogue's, so what this tests is our fitting configuration.

The detections are independent, so they are fitted njobs at a time. Each one
gets its own Sherpa session inside its own process, because the session and
the XSPEC settings behind xsphabs are global to a process.

The output table carries our value and theirs side by side, so that the
comparison is against their published 68% limits and not against an
arbitrary tolerance.
"""

import csv
import glob
import multiprocessing
import os
import sys
import traceback

base_path = "/home/marina/Doctorado/2026/HMXB_project"
detections = os.path.join(base_path, "work", "Output", "csc_fitted_detections.tsv")
products = os.path.join(base_path, "work", "csc_products")
csc_pull = os.path.join(base_path, "hmxb_5arcsec_csc.tsv")
out = os.path.join(base_path, "work", "Output", "stage1_csc_refit.tsv")

group_min = 16       # counts per channel bin, Evans et al. (2024), Sect. 3.13
band = (0.5, 7.0)    # the ACIS broad band, in which the catalogue quotes its flux
njobs = 12           # detections fitted at the same time


def number(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def read_cscview(path):
    """The CSCview TSV opens with one commented description per column, then
    the header line, then the data, with the header repeated among the rows."""
    lines = open(path).read().split("\n")
    start = next(i for i, l in enumerate(lines) if l and not l.startswith("#") and "\t" in l)
    names = [n.strip() for n in lines[start].split("\t")]
    table = {}
    for line in lines[start + 1:]:
        if not line.strip():
            continue
        values = line.split("\t")
        if len(values) != len(names):
            continue
        row = dict(zip(names, [v.strip() for v in values]))
        if row["name"] == "name":
            continue
        table.setdefault((number(row["obsid"]), row["name"]), row)
    return table


published = read_cscview(csc_pull)


def setup():
    """Once per worker process: the abundance table and the photoelectric
    cross-sections are XSPEC settings, and XSPEC is global to the process."""
    from sherpa.astro import ui
    ui.set_xsabund("angr")
    ui.set_xsxsect("bcmc")


def fit_one(target):
    """Fit one detection in its own Sherpa session and return a plain dict."""
    from sherpa.astro.ui.utils import Session
    from sherpa.astro.xspec import XSphabs
    from sherpa.models.basic import PowLaw1D

    obsid, name = int(target["obsid"]), target["name"]
    row = {"obsid": obsid, "name": name}
    pha = glob.glob(os.path.join(products, name.replace(" ", ""),
                                 "%05d_*" % obsid, "*_pha3.fits.gz"))
    if not pha:
        row["note"] = "no spectrum on disk"
        return row
    row["pha"] = os.path.basename(pha[0])

    try:
        s = Session()
        s.set_stat("chi2datavar")
        s.load_pha(1, pha[0])
        s.group_counts(1, group_min)
        s.notice(band[0], band[1])
        s.subtract(1)
        # The components are built here rather than named in a string,
        # because a bare Session does not register the model types
        absorption, powerlaw = XSphabs("abs1"), PowLaw1D("p1")
        s.set_source(1, absorption * powerlaw)
        nh, gamma, ampl = absorption.nH, powerlaw.gamma, powerlaw.ampl
        nh.val, gamma.val, ampl.val = 0.3, 2.0, 1e-4   # nH in 10^22 cm^-2
        s.set_method("levmar")
        s.fit(1)
        fit = s.get_fit_results()

        # Levenberg-Marquardt from a single starting point falls into a local
        # minimum for the most absorbed sources, where it ends with a negative
        # photon index and a reduced chi-square far above the published one.
        # The catalogue's own statistic says when that happened, so those are
        # fitted again with a global optimiser and the better fit is kept.
        target_stat = number(published.get((float(obsid), name), {}).get("powlaw_stat"))
        if target_stat is not None and fit.rstat > 1.2 * target_stat:
            row["retried"] = "moncar"
            best = (fit.statval, nh.val, gamma.val, ampl.val)
            s.set_method("moncar")
            nh.val, gamma.val, ampl.val = 0.3, 2.0, 1e-4
            s.fit(1)
            fit = s.get_fit_results()
            if fit.statval > best[0]:
                nh.val, gamma.val, ampl.val = best[1], best[2], best[3]
                s.set_method("levmar")
                s.fit(1)
                fit = s.get_fit_results()

        row.update({
            "our_nh": nh.val * 100.0,   # 10^22 -> the catalogue's 10^20 cm^-2
            "our_gamma": gamma.val,
            "our_flux": s.calc_energy_flux(band[0], band[1], id=1),
            "our_rstat": fit.rstat,
            "our_dof": fit.dof,
        })
    except Exception:
        row["note"] = "fit failed"
        row["traceback"] = traceback.format_exc().splitlines()[-1]
        return row

    cat = published.get((float(obsid), name))
    if cat is None:
        row["note"] = "no published row"
    else:
        for ours, theirs in (("csc_nh", "powlaw_nh"), ("csc_nh_lo", "powlaw_nh_lolim"),
                             ("csc_nh_hi", "powlaw_nh_hilim"), ("csc_gamma", "powlaw_gamma"),
                             ("csc_gamma_lo", "powlaw_gamma_lolim"),
                             ("csc_gamma_hi", "powlaw_gamma_hilim"),
                             ("csc_flux", "flux_powlaw"), ("csc_rstat", "powlaw_stat")):
            row[ours] = number(cat[theirs])
    return row


def inside(r, value, lo, hi):
    return all(r.get(k) is not None for k in (value, lo, hi)) and r[lo] <= r[value] <= r[hi]


if __name__ == "__main__":
    with open(detections) as fh:
        targets = list(csv.DictReader(fh, delimiter="\t"))
    print("%d detections to fit, %d at a time" % (len(targets), njobs))
    sys.stdout.flush()

    rows = []
    with multiprocessing.Pool(njobs, initializer=setup) as pool:
        for i, row in enumerate(pool.imap_unordered(fit_one, targets), 1):
            print("[%d/%d] %s %s %s"
                  % (i, len(targets), row["obsid"], row["name"],
                     row.get("note", row.get("retried", ""))))
            sys.stdout.flush()
            rows.append(row)

    rows.sort(key=lambda r: (r["obsid"], r["name"]))
    fields = ["obsid", "name", "pha", "our_nh", "our_gamma", "our_flux", "our_rstat", "our_dof",
              "csc_nh", "csc_nh_lo", "csc_nh_hi", "csc_gamma", "csc_gamma_lo", "csc_gamma_hi",
              "csc_flux", "csc_rstat", "retried", "note"]
    with open(out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("%.6g" % v if isinstance(v, float) else v) for k, v in r.items()})
    print("\n%d rows written to %s" % (len(rows), out))

    done = [r for r in rows if r.get("our_nh") is not None and r.get("csc_nh_lo") is not None]
    n_nh = sum(1 for r in done if inside(r, "our_nh", "csc_nh_lo", "csc_nh_hi"))
    n_g = sum(1 for r in done if inside(r, "our_gamma", "csc_gamma_lo", "csc_gamma_hi"))
    print("inside the published 68%% limits: N_H %d/%d, Gamma %d/%d"
          % (n_nh, len(done), n_g, len(done)))
