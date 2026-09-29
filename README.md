# X-ray properties of Galactic HMXBs corrected by interstellar absorption

Marina Belén Badaracco · Professors: Dr. Matías Zaldarriaga & Dr. Nahuel Miron Granese
Work carried out with Claude Code 2.1.278 (Opus 5, effort high).

**Presented from two things:** [`index.html`](index.html) and [`summary.pdf`](summary.pdf).
Everything else in this repository is behind them, and is meant to be readable by an agent that
has never spoken to the author.

> [!WARNING]
> **Change the paths before running anything.** The scripts and notebooks were run on the
> author's machine and contain absolute paths to it. Two prefixes appear, in `work/` and
> `Observations/`:
>
> | Path in the files | What it stands for |
> |---|---|
> | `/home/marina/Doctorado/2026/HMXB_project` | the root of this repository |
> | `/home/marina/Software/SAS/CCF` | the SAS calibration files (`SAS_CCFPATH`) |
>
> From the root of your copy, this replaces both (set your own CCF location first):
>
> ```
> CCF=/path/to/your/SAS/CCF
> grep -rlE '/home/marina' work Observations | xargs sed -i \
>   -e "s|/home/marina/Doctorado/2026/HMXB_project|$PWD|g" \
>   -e "s|/home/marina/Software/SAS/CCF|$CCF|g"
> grep -rn '/home/marina' work Observations   # should print nothing
> ```

> [!IMPORTANT]
> **Software versions the results were produced with.** Other versions may give different
> reductions and fits.
>
> | Software | Version | Calibration |
> |---|---|---|
> | CIAO | 4.17.0 | CALDB 4.12.2 |
> | SAS | 22.1.0-a8f2c2afa-20250304 | CCF set dated 2 February 2026 |
> | HEASoft | 6.36, with XSPEC 12.15.1 | — |
> | BXA / UltraNest | 5.1.1 / 4.5.2 (in `work/env_bxa`) | — |
> | SAOImage ds9 | 8.6 | — |
> | Python | 3.12 (`work/environment.yml`) | — |

The goal is the local (circumstellar) column density of Galactic high-mass X-ray binaries as a
function of **orbital phase**, obtained by fitting our own spectra and subtracting the
interstellar column from the 3D reddening cube of Doroshenko et al. (2024) through a Monte Carlo.
Reaching it needs two homogenizations that are independent of each other: of the X-ray
measurement, because the two archival catalogues fit different models over different bands with
different statistics; and of the time coordinate, because the published ephemerides use ten
different definitions of phase zero.

## Where to start reading

| File | What it is |
|---|---|
| `index.html` | the page this is presented from |
| `summary.pdf` | the four-page summary (source: `work/summary/summary.tex`) |
| `work/PROVENANCE.md` | **the record**: every result, the file and function that produced it, the check behind it, and what was *not* checked |
| `work/README.md` | what each script in `work/` does |
| `work/1_Crossmatch.html`, `work/2_XrayCrossmatch.html` | the two notebooks, exported |

`PROVENANCE.md` is the entry point for anyone checking the work. It is organised by result
identifier (`S-…` sample, `X-…` cross-match and images, `C-…` what the catalogues report,
`R-…` reduction and fits), and every entry carries a **Checked** and, where it
applies, a **Not checked** paragraph.

## Environment

```
conda env create -f work/environment.yml -p work/env
work/env/bin/python -m ipykernel install --user --name hmxb
```

Three external packages are used but not installed by that file, because they are large and
mission-specific: **CIAO 4.17** (`chandra_repro`, `dmcopy`, `wavdetect`, `mkpsfmap`), **SAS 22.1.0**
(`emproc`, `epproc`, `evselect`, `edetect_chain`) and **HEASoft 6.36/XSPEC 12.15.1** with **BXA 5.1.1**
and **UltraNest 4.5.2** in a second venv, `work/env_bxa`, created on HEASoft's python with
`--system-site-packages`. `work/setup_texlive.sh` installs the LaTeX toolchain into
`work/texlive` for the summary. The figure tool needs **SAOImage ds9 8.6**
and its XPA utilities.

## Reproducing it

From `work/`, in this order. Each script skips what is already on disk, so it is safe to rerun.

```
env/bin/python build_notebook.py          # 1_Crossmatch.ipynb, the sample
env/bin/python build_notebook_2.py        # 2_XrayCrossmatch.ipynb, the counterparts
env/bin/python match_csc_5arcsec.py       # the Chandra match and the obsid list
env/bin/python build_crossmatch_table.py  # Output/hmxb_crossmatch.csv
env/bin/python build_orbital_table.py     # Output/hmxb_orbital_parameters.csv
env/bin/python download_f23_refs.py       # the papers F23 cites
env/bin/python select_stage1_targets.py   # which catalogue fits to repeat
env/bin/python download_csc_products.py   # (CIAO)      the CSC spectra
./download_xmm_pps.sh                     #             the XMM pipeline spectra
env/bin/python download_epic_rmf.py       #             the canned EPIC responses
env/bin/python refit_csc_stage1.py        # (CIAO)      Stage 1, Chandra
./env_bxa/bin/python refit_xmm_stage1.py  # (HEASoft)   Stage 1, XMM-Newton
env/bin/python build_image_tables.py      # Observations/Images/, the per-system tables
env/bin/python build_mosaic_list.py       # Output/mosaic_systems*.txt
env/bin/python dedup_image_products.py --apply   # share the products of shared observations
```

Then, under `Observations/`:

```
Chandra/download_chandra.sh               XMM-Newton/download_xmm-newton.sh
Chandra_repro/chandra_repro.sh            XMM-Newton_repro/xmm-newton_repro.sh
                                          XMM-Newton_repro/xmm-newton_flag-pattern.sh
Images/chandra_images.sh                  Images/xmm-newton_images.sh
Images/ds9_phase_mosaic.sh <system>       # one system's observations, in phase order
```

`build_crossmatch_table.py` needs the network once, to ask Gaia DR3 whether each published
identifier exists. The summary is rebuilt with
`work/texlive/bin/x86_64-linux/pdflatex` from `work/summary/summary.tex`.

## What is in here, and what is not

Included: all code, the curated inputs, the derived tables the pipeline consumes, the two
notebooks and their exports, the summary source and its PDF, and the full provenance record.
Every file under `work/Output/` is read by at least one script or notebook; run logs and
diagnostic tables are left out and are regenerated by rerunning the step that made them.

Not included, with how to obtain each:

| Not here | Size | How to get it |
|---|---|---|
| `Observations/*_repro/`, `Observations/Images/` products | 367 GB | the download and reprocessing scripts above |
| `5XMM-DR15/5XMM_DR15.fits.gz` | 2.4 GB | <https://xmmssc.irap.omp.eu/Catalogue/5XMM-DR15/> |
| `Papers/`, `work/f23_refs/` | 186 MB | published papers: those in `Papers/` are listed in [`Papers.txt`](Papers.txt); `download_f23_refs.py` fetches the F23 references |
| `work/env`, `work/env_bxa`, `work/texlive` | 3.7 GB | `environment.yml` and `setup_texlive.sh` |
| `work/csc_products`, `work/epic_rmf`, `work/xmm_pps`, `work/stage1_xmm` | 847 MB | the three download scripts above |

`catalogues/` holds the two inputs that are small enough to ship: the CSCview 5-arcsec pull
(`hmxb_5arcsec_csc.tsv`, made with `cscquery.prop`) and the F23 machine-readable table
(`f23_tablea_catalog.csv`). The byte-range version F23 publishes, with its ReadMe, is in
`work/Input/HMXB/`.

**`Papers.txt` lists the papers that should be stored in `Papers/`.** They are the papers this
work cites, and the only ones it may cite (see `Agreement.md`); they are not shipped because they
are not ours to redistribute. Download each one (by the DOI or arXiv identifier given in the list, where there is one) and
save it as a PDF in a `Papers/` folder at the root of the repository. The references for the
individual orbital parameters are not in that list: they are F23's, and are carried per system in
`work/Output/hmxb_orbital_parameters.csv`; `download_f23_refs.py` fetches them into
`work/f23_refs/`.

## A note on how it was made

`Agreement.md` is the working agreement the agent operated under: papers are cited only if the PDF
is in `Papers/`; everything the agent produced lives in `work/`; nothing is installed into the base
environment; a check is preferred to a claim; and it never claims to have run something it did not
run. `PROVENANCE.md` is written to the same standard, which is why it records failures — a
withdrawn preprint whose six values were discarded, a catalogue value attributed to the wrong
system, a detection run that filled `/tmp`, a display that hung on a gigapixel image — alongside
the results.
