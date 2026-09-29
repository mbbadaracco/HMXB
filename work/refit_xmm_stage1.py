#!/usr/bin/env python
"""Stage 1 for XMM-Newton: repeat the 5XMM-DR15 fits on the catalogue's data.

Run with the BXA environment, from work/, after initialising HEASoft:

    heainit
    ./env_bxa/bin/python refit_xmm_stage1.py [srcid ...]

For each of our stacked sources that carries a SPEC_ fit, the pipeline
spectra of its detections are read, the highest signal-to-noise detection of
each instrument is kept, and the pair is fitted exactly as Webb et al. (2026)
Sect. 6 describe: cflux*phabs*powerlw explored with BXA/UltraNest under a
log-uniform prior on N_H over 1e19-1e25 cm^-2, a uniform prior on Gamma over
[1,3] and a uniform prior on log10(flux) over [-15,-9], with the Cash
statistic and, when pn and MOS are both present, a free MOS normalisation.

The spectra, the backgrounds and the ARFs are the pipeline's; the RMFs are
the canned ones the spectra name, fetched by download_epic_rmf.py. Nothing is
extracted here, so what this tests is our fitting configuration.

The sources are independent, so they run njobs at a time, each in its own
process: XSPEC keeps one global state per process.
"""

import csv
import glob
import multiprocessing
import os
import shutil
import subprocess
import sys
import traceback

import numpy as np
from astropy.io import fits

base_path = "/home/marina/Doctorado/2026/HMXB_project"
detections = os.path.join(base_path, "work", "Output", "xmm_spectra_detections.tsv")
pps = os.path.join(base_path, "work", "xmm_pps")
rmf_path = os.path.join(base_path, "work", "epic_rmf")
ours = os.path.join(base_path, "work", "Output", "hmxb_5xmm_crossmatch.fits")
workdir = os.path.join(base_path, "work", "stage1_xmm")
out = os.path.join(base_path, "work", "Output", "stage1_xmm_refit.tsv")

band = tuple(float(x) for x in os.environ.get("FIT_BAND", "0.3,10.0").split(","))
# The catalogue states 0.3-10 keV only for its Levenberg-Marquardt pre-screen and
# never states the band of the Bayesian fit; FIT_BAND exists to test both.
# SPEC_FLUX_PL is quoted over the fitted band, not over the catalogue's
# 0.2-12 keV detection band: rescaling our 0.2-12 keV flux to 0.3-10 keV with
# the fitted model reproduces the published value to better than 1% for every
# source tested, which is what removed a 16% offset.
flux_band = band
njobs = 6              # sources fitted at the same time


def signal_to_noise(spectrum):
    """N / sqrt(2T - N), the definition of Webb et al. (2026), Sect. 6, with
    N the net counts and T the total counts of the extracted spectrum."""
    header = fits.getheader(spectrum, 1)
    counts = fits.getdata(spectrum, 1)["COUNTS"]
    background = os.path.join(os.path.dirname(spectrum), header["BACKFILE"])
    hb = fits.getheader(background, 1)
    cb = fits.getdata(background, 1)["COUNTS"]
    total = float(counts.sum())
    scaled = float(cb.sum()) * header["BACKSCAL"] / hb["BACKSCAL"]
    net = total - scaled
    if total <= 0 or 2 * total - net <= 0:
        return None, total, net, background
    return net / np.sqrt(2 * total - net), total, net, background


def spectra_of(srcid, rows):
    """The EPIC spectra of every detection of one stacked source."""
    found = []
    for r in rows:
        tag = "%03X" % int(r["srcnum"])  # the source number is hexadecimal in the name
        for f in sorted(glob.glob(os.path.join(pps, r["obsid"], "pps",
                                               "*SRSPEC?%s.FTZ" % tag))):
            camera = os.path.basename(f)[11:13]
            if camera not in ("PN", "M1", "M2"):
                continue  # EPIC only, the RGS spectra of the same source are not used
            sn, total, net, background = signal_to_noise(f)
            if sn is None or total <= 0 or net <= 0:
                continue  # the catalogue requires positive total, background and net counts
            found.append({"file": f, "obsid": r["obsid"], "camera": camera,
                          "instrument": "pn" if camera == "PN" else "MOS",
                          "sn": sn, "total": total, "net": net, "background": background})
    return found


def group_spectrum(spectrum, background, destination):
    """One or more counts per bin, as the catalogue re-bins its spectra.

    The GROUPING column is rewritten here rather than with ftgrouppha, which
    in this installation insists on a terminal and exits before doing any
    work. A group opens at the first channel that follows a closed one and
    closes as soon as it holds one count, which is the OGIP convention: 1
    marks the first channel of a group and -1 every channel after it."""
    os.makedirs(destination, exist_ok=True)
    grouped = os.path.join(destination,
                           os.path.basename(spectrum).replace(".FTZ", "_grp.fits"))
    shutil.copy(background, os.path.join(destination, os.path.basename(background)))
    if os.path.exists(grouped):
        return grouped
    with fits.open(spectrum) as h:
        counts = np.asarray(h[1].data["COUNTS"], dtype=float)
        grouping = np.full(counts.size, -1, dtype=np.int16)
        accumulated = 0.0
        opening = True
        for i, c in enumerate(counts):
            if opening:
                grouping[i] = 1
                accumulated = 0.0
                opening = False
            accumulated += c
            if accumulated >= 1.0:
                opening = True
        h[1].data["GROUPING"] = grouping
        # XSPEC resolves BACKFILE, RESPFILE and ANCRFILE against the working
        # directory and prompts for whatever it cannot find, which blocks for
        # ever in a worker process. They are blanked here and the three files
        # are attached explicitly, by absolute path, after loading.
        for keyword in ("BACKFILE", "RESPFILE", "ANCRFILE"):
            h[1].header[keyword] = "none"
        h.writeto(grouped, overwrite=True)
    return grouped


published = {}


def fit_source(job):
    srcid, rows = job
    result = {"srcid": srcid}
    try:
        import xspec
        from bxa.xspec.solver import BXASolver
        import bxa.xspec as bxa

        xspec.Xset.chatter = 0
        xspec.Xset.logChatter = 0
        xspec.Fit.statMethod = "cstat"

        found = spectra_of(srcid, rows)
        if not found:
            result["note"] = "no usable spectrum"
            return result
        # one pn and one MOS, the highest signal-to-noise of each
        chosen = []
        for instrument in ("pn", "MOS"):
            candidates = [s for s in found if s["instrument"] == instrument]
            if candidates:
                chosen.append(max(candidates, key=lambda s: s["sn"]))
        # Which camera the catalogue used is not published, but SPEC_DOF_PL
        # is: it equals the number of fitted bins minus the free parameters,
        # three for one instrument and four for two. Counting our own bins in
        # each combination identifies the catalogue's choice exactly, for all
        # eleven sources, so the selection follows it instead of the
        # highest-signal-to-noise rule when the two disagree.
        destination = os.path.join(workdir, str(srcid))
        published_dof = published.get(srcid, {}).get("xmm_dof")
        if published_dof is not None and chosen:
            counts = {}
            for s in chosen:
                xspec.AllData.clear()
                grouped = group_spectrum(s["file"], s["background"], destination)
                spec = xspec.Spectrum(grouped)
                header = fits.getheader(s["file"], 1)
                spec.response = os.path.join(rmf_path, header["RESPFILE"].strip())
                spec.ignore("**-%f %f-**" % band)
                counts[s["instrument"]] = len(spec.noticed)
            xspec.AllData.clear()
            options = []
            for s in chosen:
                options.append(([s], counts[s["instrument"]] - 3))
            if len(chosen) > 1:
                options.append((chosen, sum(counts.values()) - 4))
            best, dof = min(options, key=lambda o: abs(o[1] - published_dof))
            result["dof_match"] = dof - published_dof
            chosen = best
        result["n_spectra"] = len(chosen)
        result["chosen"] = ";".join("%s:%s:%s:%.2f" % (s["instrument"], s["obsid"],
                                                       s["camera"], s["sn"])
                                    for s in chosen)

        xspec.AllData.clear()
        xspec.AllModels.clear()
        # One data group per instrument, so that the two can carry different
        # normalisations; with a single group they would share every parameter.
        grouped = [group_spectrum(s["file"], s["background"], destination)
                   for s in chosen]
        xspec.AllData(" ".join("%d:%d %s" % (i, i, g)
                               for i, g in enumerate(grouped, 1)))
        for i, s in enumerate(chosen, 1):
            spec = xspec.AllData(i)
            header = fits.getheader(s["file"], 1)
            spec.response = os.path.join(rmf_path, header["RESPFILE"].strip())
            spec.response.arf = os.path.join(os.path.dirname(s["file"]),
                                             header["ANCRFILE"].strip())
            spec.background = os.path.join(destination, os.path.basename(s["background"]))
            spec.ignore("**-%f %f-**" % band)

        model = xspec.Model("constant*cflux*phabs*powerlaw")
        model.constant.factor.values = [1.0, -1]      # pn fixed at unity
        model.cflux.Emin.values = [flux_band[0], -1]
        model.cflux.Emax.values = [flux_band[1], -1]
        model.cflux.lg10Flux.values = [-12.0, 0.05, -15, -15, -9, -9]
        model.phabs.nH.values = [1.0, 0.01, 1e-3, 1e-3, 1e3, 1e3]   # 10^22 cm^-2
        model.powerlaw.PhoIndex.values = [2.0, 0.01, 1.0, 1.0, 3.0, 3.0]
        model.powerlaw.norm.values = [1.0, -1]        # cflux carries the normalisation

        parameters = [model.phabs.nH, model.powerlaw.PhoIndex, model.cflux.lg10Flux]
        priors = [bxa.create_loguniform_prior_for(model, model.phabs.nH),
                  bxa.create_uniform_prior_for(model, model.powerlaw.PhoIndex),
                  bxa.create_uniform_prior_for(model, model.cflux.lg10Flux)]
        if len(chosen) > 1:
            # XSPEC ties the parameters of the second data group to the first;
            # only the MOS normalisation is released, the pn one stays at unity.
            mos = xspec.AllModels(2)
            mos.constant.factor.link = ""
            mos.constant.factor.values = [1.0, 0.01, 0.1, 0.1, 10.0, 10.0]
            mos.constant.factor.frozen = False
            parameters.append(mos.constant.factor)
            priors.append(bxa.create_uniform_prior_for(mos, mos.constant.factor))

        # The prior helpers already return the complete transformation
        solver = BXASolver(transformations=priors,
                           outputfiles_basename=os.path.join(destination, "bxa"))
        analysis = solver.run(resume=True)

        posterior = np.array(analysis["samples"])
        names = [t["name"] for t in solver.transformations]
        for i, n in enumerate(names):
            column = posterior[:, i]
            if n == "log(nH)":
                # BXA stores the log-uniform parameter as log10 of its value:
                # its transform returns x*spread + log10(bottom), and the
                # 10**x that XSPEC receives is the aftertransform. phabs
                # carries nH in 10^22 cm^-2, the catalogue quotes cm^-2.
                column = 10 ** column * 1e22
            lo, med, hi = np.percentile(column, [15.87, 50, 84.13])
            key = {"log(nH)": "our_nh", "nH": "our_nh", "PhoIndex": "our_gamma",
                   "lg10Flux": "our_lg10flux", "factor": "our_iin"}.get(n, "our_" + n)
            result[key] = med
            result[key + "_lo"] = lo
            result[key + "_hi"] = hi
        result["our_flux"] = 10 ** result["our_lg10flux"]
        result["logz"] = analysis["logz"]
    except Exception:
        result["note"] = "fit failed"
        result["traceback"] = traceback.format_exc().splitlines()[-1]
    return result


def published_values():
    """The catalogue's own values, read once and inherited by the workers."""
    from astropy.table import Table
    t = Table.read(ours)
    out = {}
    for r in t:
        nh = float(r["SPEC_NH_PL"])
        if not np.isfinite(nh):
            continue
        out[int(r["SRCID"])] = {
            "xmm_nh": nh, "xmm_gamma": float(r["SPEC_GAMMA_PL"]),
            "xmm_flux": float(r["SPEC_FLUX_PL"]), "xmm_iin": float(r["SPEC_IIN_PL"]),
            # The ERR_LO and ERR_UP columns are the lower and upper bounds of
            # the credible interval themselves, not offsets from the median:
            # for 3060385010100002 they read 1.02e22 and 1.56e22 around a
            # median of 1.27e22, and the Gamma pair brackets its median too.
            "xmm_nh_lo": float(r["SPEC_NH_ERR_LO_PL"]),
            "xmm_nh_hi": float(r["SPEC_NH_ERR_UP_PL"]),
            "xmm_gamma_lo": float(r["SPEC_GAMMA_ERR_LO_PL"]),
            "xmm_gamma_hi": float(r["SPEC_GAMMA_ERR_UP_PL"]),
            "xmm_pvalue": float(r["SPEC_PVALUE_PL"]), "xmm_dof": float(r["SPEC_DOF_PL"]),
        }
    return out


if __name__ == "__main__":
    with open(detections) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    by_source = {}
    for r in rows:
        by_source.setdefault(int(r["srcid"]), []).append(r)
    if len(sys.argv) > 1:
        wanted = {int(a) for a in sys.argv[1:]}
        by_source = {k: v for k, v in by_source.items() if k in wanted}

    os.makedirs(workdir, exist_ok=True)
    published.update(published_values())
    jobs = sorted(by_source.items())
    print("%d sources to fit, %d at a time" % (len(jobs), njobs))
    sys.stdout.flush()

    results = []
    with multiprocessing.Pool(min(njobs, len(jobs))) as pool:
        for i, r in enumerate(pool.imap_unordered(fit_source, jobs), 1):
            print("[%d/%d] %s %s" % (i, len(jobs), r["srcid"],
                                     r.get("note", r.get("chosen", ""))))
            sys.stdout.flush()
            results.append(r)

    for r in results:
        r.update(published.get(r["srcid"], {}))
    results.sort(key=lambda r: r["srcid"])

    fields = ["srcid", "n_spectra", "chosen",
              "our_nh", "our_nh_lo", "our_nh_hi", "our_gamma", "our_gamma_lo",
              "our_gamma_hi", "our_lg10flux", "our_flux", "our_iin", "logz",
              "xmm_nh", "xmm_nh_lo", "xmm_nh_hi", "xmm_gamma", "xmm_gamma_lo",
              "xmm_gamma_hi", "xmm_flux", "xmm_iin", "xmm_pvalue", "xmm_dof",
              "dof_match", "note", "traceback"]
    with open(out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for r in results:
            w.writerow({k: ("%.6g" % v if isinstance(v, float) else v)
                        for k, v in r.items()})
    print("\n%d rows written to %s" % (len(results), out))
