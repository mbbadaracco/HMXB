# work/

Everything produced for this project lives here.

## Layout

```
1_Crossmatch.ipynb       sample construction: F23 x F24 duplicates,
                         Gaia DR3 counterparts, SIMBAD identifiers.
1_Crossmatch.html        the same, exported.
build_notebook.py        emits 1_Crossmatch.ipynb.
2_XrayCrossmatch.ipynb   Chandra (CSC 2.1.1) and XMM-Newton (5XMM-DR15)
                         counterparts, and their identifiers against F23
                         and SIMBAD.
2_XrayCrossmatch.html    the same, exported.
build_notebook_2.py      emits 2_XrayCrossmatch.ipynb.
_helpers_cell.py         two helper functions of the notebook, kept in a
                         file only because they contain nested quotes.
crossmatch.py            positional cross-match: haversine separation,
                         rotation of the error ellipses, Mahalanobis
                         distance.  Pineau et al. (2011), Appendix A.
check_rotation_formula.py  the sympy check of that rotation.
match_csc_5arcsec.py     the Chandra match on the 5 arcsec pull, and the
                         observation list it implies.
build_crossmatch_table.py  Output/hmxb_crossmatch.csv, the final table:
                         our counterparts, F23's and SIMBAD's, the
                         obsids, the orbital period and the notes.
download_f23_refs.py     the papers F23 cites for the orbital periods.
select_stage1_targets.py which catalogue fits Stage 1 repeats.
download_csc_products.py the CSC spectra and responses.  Needs CIAO.
download_xmm_pps.sh      the XMM-Newton pipeline spectra.
download_epic_rmf.py     the canned EPIC responses the pipeline omits.
refit_csc_stage1.py      Stage 1 for Chandra, 12 at a time.  Needs CIAO.
refit_xmm_stage1.py      Stage 1 for XMM-Newton, 6 at a time, with BXA
                         and UltraNest in env_bxa.  Needs HEASoft.
environment.yml          the conda environment.
setup_texlive.sh         installs the LaTeX toolchain into work/texlive.
PROVENANCE.md            where every result comes from, and what has
                         NOT been checked.
Input/                   catalogues and VizieR query results (read-only).
Output/                  the tables the scripts and notebooks write.
f23_refs/                the papers F23 cites, fetched for reading only.
csc_products/            CSC 2.1 per-detection spectra and responses.
xmm_pps/                 XMM-Newton pipeline spectra and ARFs.
epic_rmf/                canned EPIC response matrices.
stage1_xmm/              grouped spectra and BXA chains, one per source.
env/                     the conda environment.
env_bxa/                 XSPEC, BXA and UltraNest, on HEASoft's python.
texlive/                 TeX Live, installed by setup_texlive.sh.
```

## Reproducing

```bash
cd work
conda env create -f environment.yml -p ./env
./env/bin/python -m nbconvert --to notebook --execute --inplace 1_Crossmatch.ipynb
./env/bin/python -m nbconvert --to html 1_Crossmatch.ipynb
./env/bin/python -m nbconvert --to notebook --execute --inplace 2_XrayCrossmatch.ipynb
./env/bin/python -m nbconvert --to html 2_XrayCrossmatch.ipynb
```

`2_XrayCrossmatch.ipynb` decompresses `5XMM_DR15.fits.gz` (8.5 GB) into
the system temporary directory and deletes it again at the end.  Set
`XMM_FITS_CACHE` to put it elsewhere, or `KEEP_XMM_FITS = True` in the
notebook to keep it.

and for the summary

```bash
./setup_texlive.sh          # only once
cd summary && ../texlive/bin/x86_64-linux/pdflatex summary.tex
```

## Notes

- The notebook queries the Gaia TAP service, VizieR and SIMBAD, so it
  needs a network. `gaia_data()` sends a `LEFT JOIN` against
  `gaiadr3.astrophysical_parameters` for a few hundred `source_id`s, and
  Gaia's *synchronous* endpoint answers that with
  `Error 408: Job timeout/aborted` — reproducibly, not occasionally. The
  notebook therefore redirects `Gaia.launch_job` to the asynchronous
  endpoint, which is the one meant for queries of that size, with a
  retry on top. Same query, same service, longer server-side budget.
- Results that cost a query are cached under `Output/`. Delete the
  relevant file to force a re-query.
- `crossmatch.py` is a byte-identical copy of `AppendixB.py` from
  `~/Doctorado/2026/FotometriaSintetica`; only the file name differs.
  See `PROVENANCE.md`, entry S-01.
- `Input/` is treated as read-only. The VizieR `.tsv` files carry the
  request URL and the query date in their own headers.
- Any Gaia `source_id` read from a file needs an explicit `dtype=str`.
  A 19-digit identifier does not survive `float64`, and nothing raises
  when it is lost — see `PROVENANCE.md`, entry S-06.
