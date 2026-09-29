# Provenance

Where every statement in the manuscript comes from, and what has **not** been
checked. It is organised by what the statement is about, not by when the work
was done: an entry is the current account of one thing, rewritten when that
thing changes rather than appended to.

**Evidence labels.** *read* means a cited prescription was read in the paper
named. *checked* means the evidence was produced here, by a script that is in
`work/` and can be rerun. *not checked* names something that was not
established, and every entry that has one ends with it.

**Standing rules.** No result rests on the agreement of a single fit
statistic. An identifier is read as a string, never through a float. Every
process that can run in parallel does. No file is overwritten, and a step
skips an output that already exists. Nothing is cited that is not a PDF in
`Papers/`.

## How to reproduce

From `work/`, with the environment of `environment.yml` in `env/`:

```
env/bin/python build_notebook.py        # 1_Crossmatch.ipynb, the sample
env/bin/python build_notebook_2.py      # 2_XrayCrossmatch.ipynb, the counterparts
env/bin/python match_csc_5arcsec.py     # the Chandra match, and the obsid list
env/bin/python build_crossmatch_table.py  # Output/hmxb_crossmatch.csv, the counterparts
env/bin/python build_orbital_table.py   # Output/hmxb_orbital_parameters.csv, the orbits
env/bin/python download_f23_refs.py     # the papers F23 cites, into f23_refs/
env/bin/python select_stage1_targets.py # which catalogue fits to repeat
env/bin/python download_csc_products.py # (needs CIAO)  the CSC spectra
./download_xmm_pps.sh                   # the XMM-Newton pipeline spectra
env/bin/python download_epic_rmf.py     # the canned EPIC responses
env/bin/python refit_csc_stage1.py      # (needs CIAO)  Stage 1, Chandra
./env_bxa/bin/python refit_xmm_stage1.py  # (needs HEASoft) Stage 1, XMM-Newton
env/bin/python build_image_tables.py    # Observations/Images/, the per-system tables
env/bin/python build_mosaic_list.py     # Output/mosaic_systems*.txt, what can be drawn
env/bin/python dedup_image_products.py --apply   # share the products of observations covering two systems
```

and, under `Observations/`, `download_chandra.sh`, `chandra_repro.sh`,
`download_xmm-newton.sh`, `xmm-newton_repro.sh` and
`xmm-newton_flag-pattern.sh`; then, under `Observations/Images/`,
`chandra_images.sh` and `xmm-newton_images.sh` (X-11) to make the images and
detect the sources, and `ds9_phase_mosaic.sh <system>` (X-12) to look at one
system's observations in phase order.


`build_crossmatch_table.py` needs the network once, to ask Gaia DR3 whether
each published identifier exists.


## 1. What the project is for

### I-01 — what F23 and N23 contain

- **Claim.** F23: 152 HMXBs, 111 Gaia DR3 counterparts, no X-ray flux /
  column / spectral slope. N23: 169 systems, fluxes from five
  instruments in five bands plus one column taken from the *Swift*
  catalogue.
- **Source.** `Papers/Fortin2023.pdf` (abstract) and
  `Papers/Neumann2023.pdf` (abstract; Sect. 2.2 *X-ray properties*; the
  column-by-column description of the catalogue, Col. 30 for the column
  density and Cols. 31–41 for the fluxes).
- **How checked.** read. N23 Sect. 2.2 names the bands: 0.2–12 keV for
  *XMM-Newton*, ACIS broad 0.5–7.0 keV or HRC wide 0.1–10.0 keV for
  *Chandra*, 0.3–10.0 keV for *Swift*/XRT, plus *Swift*/BAT and
  *INTEGRAL*. The column density in Col. 30 is described as taken from
  the *Swift* point-source catalogue, where it was obtained by
  interpolating a power law through hardness ratios.
- **Not checked.** The claim about F23 was read from the paper. The
  machine-readable catalogue (`Extra/catalog.csv`) was *not* gone
  through column by column for this draft.

### I-02 — F23 distances come from Bailer-Jones et al. (2021)

- **Claim.** The distances F23 tabulates are those of Bailer-Jones et al.
  (2021), queried through the Gaia identifier of the counterpart.
- **Source.** `Papers/Fortin2023.pdf`, Sect. 3 and the notes to
  Table A.1: "Distance is queried from Bailer-Jones et al. (2021)".
- **How checked.** read.

### I-03 — what the 3D N_H tool is

- **Claim.** `\citet{doroshenko2024}` combines a 3D optical reddening
  map with far-infrared dust emission and radio dispersion-measure
  modelling into a pan-Galactic cumulative-reddening cube, calibrated
  against independent X-ray absorption measurements, and distributed
  with the E(B−V) → N_H conversion factors for two abundance tables.
- **Source.** `Papers/Doroshenko2024.pdf`, abstract, Sect. 2 and Sect. 4.
- **How checked.** read.
- **Not checked.** Nothing has been queried from the cube for this
  draft. The numerical calibration constants are *not* quoted in the
  manuscript yet, and will carry their own entry when they are.

---

### T-03 — project and journal requirements

`../final-project.html` requires reproducible code, data, environment and
checks, with provenance, an HTML presentation and a PDF, ultimately in a
GitHub repository. No commit or push is authorized by this task.
The earlier manuscript used the supplied `aa.cls` v9.4, `aa.bst` and official
example. The earlier request to the A&A author page returned HTTP 403;
compliance with all prose author guidelines has not been established.

## 2. Inputs, software and bibliography

### T-01 / T-02 — environments and manuscript build

The existing project environment is `env/`; the XMM fitter uses `env_bxa/`
with HEASoft/PyXspec; the Chandra fitter uses CIAO's Python and Sherpa.
No package is installed into a base environment: `env_bxa/` is a venv on
HEASoft's python with `--system-site-packages`, holding XSPEC 12.15.1,
UltraNest 4.5.2 and BXA 5.1.1. TeX Live is installed in `texlive/`.
The earlier session recorded a successful `../manuscript_aanda/Makefile`
`make check` (no undefined references/citations or LaTeX errors).
That is a dated build check, not a claim that the manuscript reflects the
fits. The manuscript has not been rebuilt or updated since.

### T-04 — reduction calibration

User-supplied versions: CIAO 4.17.0 / CALDB 4.12.2; SAS
22.1.0-a8f2c2afa-20250304; CCF set dated 2 February 2026. The previous session
recorded 548 CCF files (2.8 GB). No CCF update was completed. That calibration
record concerns our event reprocessing; Stage 1 uses the catalogue responses.
Whether subsequent CCF changes affect this sample remains unchecked.

### B-01 — bibliography audit recorded on 25 September 2026

- **How checked.** Each of the 21 entries of
  `manuscript_aanda/references.bib` carries a `pdf =` field, and every
  one of those files exists in `Papers/`. The converse does not hold and
  is not required: `Papers/` now holds 63 PDFs, most of them read but not
  cited. Two entries (`haberl2016`, `prokhorenko2026`) are present but
  not yet cited.

### B-02 — publication metadata resolved in the earlier session

- **How checked.** resolved. The PDFs of Evans et al., Doroshenko,
  Sazonov & Khabibullin, Haberl & Sturm and Prokhorenko et al. are
  preprints or carry no journal reference on their title page. Their
  publication data were resolved against the CrossRef REST API and the
  arXiv API:
  - Evans et al. 2024 → arXiv:2407.10799, **no journal reference**;
    cited as a preprint.
  - Doroshenko 2024 → arXiv:2403.03127, **no journal reference**;
    cited as a preprint.
  - Sazonov & Khabibullin 2017 → MNRAS, 468, 2249.
  - Haberl & Sturm → A&A, 586, A81, published **2016**, although the PDF
    is the 2015 preprint and the file is named `Haberl2015.pdf`.
  - Prokhorenko et al. → J. High Energy Astrophys., 53, 100611 (2026).
  - Webb et al. (5XMM-DR15) is a draft dated 5 June 2026 with no journal
    reference; cited as in press.

---

## 3. The sample

### S-01 / S-02 — recovered code and inputs

`1_Crossmatch.ipynb` and `build_notebook.py` recover the user's
`FotometriaSintetica/0_GaiaPhotometry+.ipynb`. `crossmatch.py` is the renamed
`AppendixB.py`; byte identity was checked again on 26 September. Inputs are
F23 `Input/HMXB/{tablea.dat,ReadMe}`, F24 `Input/LMXB/{lmxbcat.dat,ReadMe}`
and the saved VizieR cone queries in `Input/OnlyCoords/`. Their query headers
are retained. Confirmed non-HMXBs 1E 1740.7-2942 and GRS 1758-258 are removed;
ambiguous pairs remain with `Notes: Ambiguous with LMXB: <ID>`, which X-08
carries into the final table in words. The original literature classification
has not been independently re-adjudicated.
Gaia TAP uses the asynchronous endpoint because the synchronous join timed out.

### S-03 / S-04 — positional statistic and F23/F24 overlap

`crossmatch.py` implements the haversine separation, covariance rotation and
Mahalanobis statistic associated with Pineau (2011), Appendix A. Error axes
enter in arcseconds; coordinates enter in degrees. `ePos` is converted from
degrees to arcseconds in the notebook before covariance construction.
`check_rotation_formula.py` compares the sine-law implementation against
vector geometry. The archived run recorded maximum angular errors of
1.975e-7 degrees for the implemented sine law and 88.46 degrees for the printed
small-angle expression on its test sample. It is a numerical test, not a
symbolic proof; it was not rerun for Stage 1.
The saved F23/F24 comparisons contain seven positional candidates and nine
identifier pairs, with the two extra pairs separated by 1003 and 2156 arcsec.
These are candidate/identifier associations, not seven statistically accepted
positional matches. Source tables: `Output/f23_f24_{positional,identifiers}.csv`.

### S-05 / S-06 / S-07 — Gaia sample and adoption

The current saved comparison and adoption tables have **113** rows; the
distance-qualified sample and `hmxb_sample_positions.vot` have **109**.
These are different selections. The Gaia crossmatch is not restricted to
F23's pre-existing Gaia identifications or to availability of a Bailer-Jones
distance. The X-ray analysis uses the latter 109-system sample.
The acceptance criterion is dM2 < 9.21; F23's 90% position error is rescaled
with 4.61, and Gaia's mas errors are converted to arcsec, retaining RA_DEC_CORR.
All adopted Gaia identifiers are the user's crossmatch identifiers, including
3A 0656-072 and sources with no F23/SIMBAD identification. `Notes` records
agreement, disagreement and absence separately. Gaia identifiers are strings:
float64 cannot preserve arbitrary 19-digit identifiers.

`Output/simbad_ids.csv` caches `Simbad.query_objectids` for the 113 sources;
`1_Crossmatch.ipynb` contains the query and three-way identifier comparison.
The earlier provenance's detailed float-rounding and Gaia-existence audit is
historical evidence in the archive, not a claim that all those diagnostic
queries are present in the current notebook. No SIMBAD/Gaia query was rerun.
The VOTable contains ID, coordinates and coordinate errors, but no correlation;
notebook 2 obtains that correlation from
`Output/gaia_data_Tables/hmxbs_gaia_data.csv`. No external schema validation
or false-association-rate estimate has been performed.

## 4. Counterparts: ours, F23's and SIMBAD's

### X-01 / X-02 — X-ray crossmatches

`2_XrayCrossmatch.ipynb` and `build_notebook_2.py` use adopted Gaia positions.
CSC master-source ellipses are 95% regions (factor 5.99), as described in the
CSC query headers. Sexagesimal coordinates are converted to degrees using
the user's SkyCoord helper. 5XMM acceptance uses per-axis
sqrt(RADEC_ERR^2/2 + 0.88^2) arcsec, following the local draft §4.1.2.
Its ellipse statistic is supplementary because the position-angle convention
has not been established. The ellipse contour uses ΔC=2.30.
The saved results contain 65 CSC sources / 243 detection rows for 58 HMXBs,
and 64 5XMM sources for 60 HMXBs, all within the 109-system sample.
All passing matches are retained. The secondary `dM2_xray_F23` column does
not define the adopted sample; equal acceptance is not proof that its errors
are negligible or that either match is independently validated.
CSC completeness is limited by the supplied 2-arcsec candidate query. XMM
uses the full catalogue with an RA/Dec prefilter. No new match was performed.

### X-07 — the Chandra cross-match repeated on a 5 arcsec cone

**The pull.** `../hmxb_5arcsec_csc.tsv`, made in CSCview 1.3.11 on 27
September and recorded by `../cscquery.prop`: a crossmatch query against
`Output/hmxb_sample_positions.vot` — our Gaia positions, `ID`, `RA_ICRS`,
`DE_ICRS` — with `pos.cm.radius = 5.0` arcsec, catalogue release 2.1, all
rows. It returns 289 rows, 73 distinct sources over 62 systems, with the
same 239 columns as the earlier 2 arcsec pull. It **replaces**
`hmxb_2arcsec_csc.tsv`, which was deleted as superseded once it was
verified that all 69 of its pairs are among the 73 of the wider pull.
`select_stage1_targets.py`, `refit_csc_stage1.py` and `build_notebook_2.py`
all read the new file.

**One column is missing.** The 5 arcsec pull does not carry
`err_ellipse_r0/r1/ang`, which the Mahalanobis statistic needs. For the 65
sources already in `Output/hmxb_csc_crossmatch.csv` the ellipse was taken
from there; the remaining 8 were fetched from the CSC 2.1 TAP service
(`csc21.master_source` at `http://cda.cfa.harvard.edu/csctap`), all 8
returned. No ellipse is invented and every one comes from the catalogue.

**The match, by `match_csc_5arcsec.py`.** Same statistic and same
threshold as notebook 2: the covariance of the Gaia position against the CSC
error ellipse rescaled from its 95% region by sqrt(5.99), accepted at
dM2 < 9.21. Of the 73 pairs, **67 are accepted over 59 systems and 6 are
rejected**. Against the 2 arcsec result — 65 pairs over 58 systems —
**nothing is lost and two are gained**:

| system | source | separation | dM2 |
|---|---|---|---|
| `AX J1841.0-0536` | `2CXO J184107.4-053450` | 4.49 arcsec | 1.42 |
| `GX 301-2` | `2CXO J122638.2-624610X` | 4.82 arcsec | 5.28 |

The first passes on a large ellipse, `err_ellipse_r0` 12.15 arcsec; the
second is the extended-hull designation of a source already in the sample.

**What this settles about the three F23 sources of X-06.** The wider cone
returns all three, so they were absent because of the 2 arcsec radius, not
because the catalogue lacks them. Two of the three then fail the statistic:

- `2CXO J184107.4-053450`, `AX J1841.0-0536`: 4.49 arcsec, ellipse 12.15
  arcsec, dM2 = 1.42 — **accepted**, now in the sample.
- `2CXO J131825.0-625815`, `IGR J13186-6257`: 3.91 arcsec against an
  ellipse of 0.37 arcsec, dM2 = 753 — **rejected**, decisively.
- `2CXO J225355.0+624337`, `2MASS J22535512+6243368`: 0.46 arcsec, but the
  ellipse is 0.31 by 0.30 arcsec, dM2 = 13.8 — **rejected**, and only
  marginally above the threshold. A half-arcsecond offset is significant
  when the position is that good.

So F23 and we differ on two of them because the statistic rejects what a
fixed radius would have accepted, not because either catalogue is wrong.

**Two rejected sources are kept by decision.** On 27 September the user
chose to keep the two sources F23 reports that the statistic rejects. They
are kept, not re-accepted: `match_csc_5arcsec.py` writes
`Output/hmxb_csc_match_5arcsec.csv`, all 73 candidate pairs, where each row
carries a `Notes` value that is `accepted by the cross-match` (67 rows),
`rejected by the cross-match at dM2 = ...` (4 rows) or, for these two,
`kept by decision: reported by F23, rejected by the cross-match at
dM2 = ... (threshold 9.21)` with its own value — 13.8 for
`2CXO J225355.0+624337` and 752.7 for `2CXO J131825.0-625815`. The sample
therefore follows two rules, and the table says which applies to each row.

**The observation list was rewritten.**
`../Observations/Chandra/obsids_chandra.txt` now holds **193** observations,
up from 187; nothing was removed. The six added are 4649
(`AX J1841.0-0536`), 9049 (`IGR J13186-6257`) and 9919, 9920, 10811, 10812
(`2MASS J22535512+6243368`) — that is, four of the six come from a system
kept by decision rather than by the criterion. None of the six is on disk,
raw or reprocessed, so `download_chandra.sh` and `chandra_repro.sh` have to
be rerun; both skip what is already there.

**Not checked.** The new acceptances are not propagated further: the
accepted tables of X-01, the counterpart table and the identifier tables
still rest on the 2 arcsec match, except that the Stage 1 scripts now read
the 5 arcsec pull for the fit columns, which took the fitted-detection count
from 115 to 116 over 48 sources. Rerunning notebook 2 on the wider pull, and
everything downstream of it, has not been done. `build_notebook_2.py` reads the 5 arcsec pull
but has not been rerun, so `Output/hmxb_csc_crossmatch.csv` and the
comparison tables it writes still hold the 2 arcsec result.

### X-03 / X-04 — names and Notes

Comparisons use catalogue names. The implemented redesignation rule also
decodes the coordinates embedded in J-names and accepts roots within 3 arcsec;
this is a heuristic, not an independently verified catalogue identifier map.
It is preserved, not silently replaced. Prefix-only and accepted redesignation
changes do not count as disagreements; the adopted names are CSC 2.1.1/5XMM.
Neither `ePos12` nor `ePos13` determines identifier agreement.
All accepted names enter comparisons, not just the closest one. The saved
tables and `hmxb_sample_notes.csv` hold the result. AX J1700-419 remains the
unadjudicated XMM disagreement. The old comparison counts remain readable in
the notebook outputs; their science has not been rerun during Stage 1.

### S-08 — the Gaia identifier of each system, ours against F23 and SIMBAD

`build_crossmatch_table.py` writes the Gaia columns of `Output/hmxb_crossmatch.csv`, 109
systems: our adopted Gaia DR3 source, the one F23 publishes, the one SIMBAD
attaches, whether each of those exists in Gaia DR3, the three agreement
flags, the separation between our source and SIMBAD's where they differ,
both G magnitudes, and a Notes column that states the cause of every
disagreement. Identifiers are read and written as strings throughout: a
column holding blanks is promoted to float by pandas, which rounds an
18-19 digit `source_id` silently, and that is how `gaia_id_comparison.csv`
came to display 4.272350e+17 for 427234969757165952.

**Checked now against `gaiadr3.gaia_source`.** 123 distinct identifiers were
queried; 117 exist. All 113 of ours exist and all 103 SIMBAD ones exist.
**Six of the 110 F23 identifiers do not exist.** In each the published value
is our identifier with its leading digits removed: one digit for
`1A 0535+262`, `AAO+28 342`, `IGR J06074+2205`, `HD 259440` and
`SAX J0635.2+0533`, two digits for `3A 0656-072`. The truncation is in the
CDS table itself, not in our reading of it: the ReadMe declares the field
`1522-1540 I19` and the raw bytes for `1A 0535+262` are
` 441207615229815040`, eighteen digits with a leading blank in a
nineteen-character field. A value that names no Gaia DR3 source is therefore
reported as absent rather than compared as a different counterpart.

**Where ours and SIMBAD genuinely differ**, four systems, both identifiers
exist and the two sources are far apart. Whenever the SIMBAD identifier does
not coincide with ours, the separation between the two Gaia sources lies
**between 3.28 and 114.43 arcsec**: `IGR J18462-0223` 3.28, `3A 0656-072`
42.54, `IGR J19113+1533` 47.76 and `AX J1841.0-0536` 114.43 arcsec. None is
a near-duplicate; they are different stars, and SIMBAD's is the brighter in
three of the four.

**Three systems carry a Gaia ID Candidate.** F23 publishes no Gaia
counterpart for `GS 0834-430`, `IGR J12341-6143` and `GS 1839-06`, while our
cross-match returns one; the Notes column marks these `Gaia ID Candidate`
rather than treating them as disagreements. `GS 0834-430` is corroborated by
SIMBAD, which gives the same source; for the other two SIMBAD lists no Gaia
DR3 identifier either, so the identification rests on our cross-match alone.
Counting these, our adopted identifier differs from what F23 publishes for
**nine** systems, not seven: six truncated values and three absent ones.

**How F23 assigned the Gaia counterpart**, read from
`../Papers/Fortin2023.pdf`, Sect. 2.2 "Finding an unambiguous chain of
counterparts", and Table 1. It is not a match of Gaia against an X-ray
position; it is the last link of a chain.

1. **The premise.** "We considered that a secure identification of an HMXB
   partly comes from having an unambiguous list of its detections from hard
   X-rays down to the near-infrared." The chain also serves to drop sources
   listed as HMXBs on the strength of a single detection decades ago.
2. **The ladder of catalogues**, in increasing astrometric precision, with
   the radius of Table 1: HEAO 1, Uhuru 4, Ariel V 3, INTEGRAL and Fermi at
   20 arcmin; BeppoSAX 6 arcmin; Einstein 2E 4 arcmin; ROSAT 35 arcsec;
   Swift 2SXPS 8 arcsec; 4XMM DR11 4 arcsec; Chandra CSC 2 3 arcsec; 2MASS
   120 mas; **Gaia DR3 20 mas**. Each radius is "about twice of the worst
   astrometric performance in the corresponding catalogue".
3. **Which radius is used.** The Table 1 value applies only when the
   starting position is more accurate than the catalogue being queried;
   when the catalogue is the more accurate of the two, "the cone size was
   set to the error available in the positional data".
4. **The chain is walked recursively.** "after reviewing the counterparts
   found at high energies, we performed a recursive search, from poorly
   accurate counterparts to the most accurate catalogues (2MASS and Gaia)",
   which recovers "the chain of detection from high energies down to the
   optical/nIR wavelengths". The Gaia source is therefore reached from the
   preceding link, typically 2MASS or a soft X-ray position, not from the
   hard X-ray position.
5. **A systematic allowance.** For Swift, XMM-Newton, Chandra, 2MASS and
   Gaia, "we added 0.5 arcsec to the positional uncertainty when validating
   the chain", because detections of one source can otherwise be formally
   incompatible; they argue "it is unlikely that two separate sources lie
   closer than 0.5 arcsec". The same value was used by Fortin et al.
   (2022b) to find unambiguous Gaia counterparts to 2MASS sources.
6. **Manual verification.** "We verified each individual result of this
   automatic counterpart search. We manually removed false detections of
   counterparts, and searched for actual counterparts in the literature
   when necessary." Coordinates taken by hand from a publication, often an
   Astronomer's Telegram, carry a reference in the online catalogue.
7. **Their own caveat.** The automatic query "can generate false
   counterparts because the typical astrometrical accuracies that we used
   are based on the worst performance" of each facility.

Their distances are a separate operation: Bailer-Jones et al. (2021) is
built on EDR3, so the DR3 identifiers "cannot be directly retrieved" against
it and F23 first obtained EDR3 identifiers "using a cone sky match".

**SIMBAD.** Its documentation index lists no page describing how Gaia
identifiers are attached to existing objects, and the database exposes no
per-identifier provenance: `ident` carries an identifier and an object
reference and no bibcode. The origin of a SIMBAD Gaia DR3 association
therefore **cannot** be recovered from the database, and none is claimed
here.

**Not checked.** Whether the four genuine disagreements are resolved in
favour of our counterpart or SIMBAD's; that needs the spectroscopic
identification of the donor, not astrometry. Whether F23's published
identifiers are truncated in the dynamic version of the catalogue or only in
the frozen VizieR one.

### X-06 — the X-ray counterpart of each system, ours against F23 and SIMBAD

`build_crossmatch_table.py` writes the X-ray columns of `Output/hmxb_crossmatch.csv`, the same
study made for Gaia in S-08, for the 109 systems of the sample. It carries,
per mission, the designation we adopt, the one F23 publishes, the one SIMBAD
attaches, whether the F23 entry is a designation of the catalogue its column
claims, the agreement flag and two separations: against the source the
one-row-per-system comparison picked, and against the **closest** of our
accepted sources.

**Chandra.** We have a counterpart for 58 systems; F23 gives one for 50 and
none for 8, which are flagged `Chandra ID Candidate`. Of the 50, **48 are
the same designation** and **2 are not CSC designations at all**: F23's
Chandra column, declared "Identifier in Chandra CSC2", gives
`2SXPS J075542.5-293353` for `SGR 0755-2933`, which is a Swift designation,
and `4U 1901+03` for `4U 1901+03`, which is the system's own historical
name. Over the whole F23 catalogue the column holds 78 entries with four
prefixes that are not CSC — 2SXPS, CXOU, GRO and 4U — of which these two
fall in our sample. Where both exist the positions agree closely: the
separation between the F23 position and the nearest of ours runs from 0.05
to 0.95 arcsec, median 0.20.

**XMM-Newton.** We have a counterpart for 60 systems; F23 gives one for 52
and none for 8, flagged `XMM-Newton ID Candidate`. All 87 entries of F23's
XMM column are well-formed 4XMM designations, so nothing is of the wrong
kind. Of the 52, 14 carry the identical name, **37 are the same source
redesignated** between 4XMM DR11 and 5XMM-DR15 — F23 works with DR11 and we
with DR15, so a changed name is expected and is not a disagreement — and
**one is a genuine disagreement**: `AX J1700-419`, ours
`5XMM J170004.2-415804` against F23's `4XMM J170004.5-415809`, 5.80 arcsec
apart. Separations run from 0.20 to 5.80 arcsec, median 0.69.

**Two apparent disagreements were artifacts, and are removed.** Comparing
F23's identifier against a single one of our sources invents a disagreement
when a system has several. `4U 1700-377` and `IGR J16328-4726` each have two
accepted 5XMM sources, and the first version of this table compared F23
against the wrong one, giving separations of 18.76 and 17.37 arcsec. Against
the closest of our sources they are 1.81 and 0.72 arcsec. The table now
carries both columns so the difference is visible rather than hidden.

**One table for the three missions.** `build_crossmatch_table.py` writes
`Output/hmxb_crossmatch.csv`, 109 systems and 19 columns: for each of
Gaia, Chandra and XMM-Newton, the designation we adopt, the one F23
publishes, the one SIMBAD attaches, the validity and agreement flags, and
the angular separation. The X-ray separations are always to the **closest**
of our accepted sources, for the reason above. Four of the 113 rows carry no
X-ray columns — `AX J1714.1-3912`, `AX J1838.0-0655`, `H 1553-542` and
`gam Cas` — because they are outside the 109-system sample that the X-ray
comparison was built on.

**The Chandra position settles `AX J1700-419`.** It is the one genuine XMM
disagreement, 5.80 arcsec between our 5XMM source and F23's 4XMM one. For
the same system the Chandra designations agree exactly, ours and F23's both
`2CXO J170004.3-415805`, separated by 0.07 arcsec. That Chandra position
lies **1.61 arcsec from our 5XMM source and 4.33 arcsec from F23's 4XMM
one**, so the more precise measurement supports our choice. This is
evidence, not proof: it assumes the Chandra source is the same object, which
is what the chain of counterparts asserts rather than demonstrates.

**How F23 assigned the X-ray counterparts.** By the same chain of cone
searches as the Gaia counterpart, documented in S-08: the ladder of Table 1
of `../Papers/Fortin2023.pdf` gives 4 arcsec for 4XMM DR11 and 3 arcsec for
Chandra CSC 2, each "about twice of the worst astrometric performance in the
corresponding catalogue", with the cone reduced to the error on the starting
position when that is the less accurate of the two, 0.5 arcsec added to the
positional uncertainty when validating the chain, and every result inspected
by hand. The soft X-ray link is not the end of the chain: F23 note that the
astrometry of Chandra and XMM-Newton "can rival optical telescopes", so
these positions are what the search then follows down to 2MASS and Gaia.

**SIMBAD.** As for Gaia, no per-identifier provenance is published, so the
origin of a SIMBAD X-ray association cannot be recovered from the database.

**Not checked.** Whether F23's 4XMM DR11 designations exist in DR11 itself;
that release is not held here, so a designation was tested for being
well-formed and for its encoded position, not for existing. Which of
`AX J1700-419`'s two candidates is the counterpart. Whether the two non-CSC
entries are errors or deliberate substitutions where F23 found no CSC
source.

### X-08 — the final cross-match table

`build_crossmatch_table.py` writes `Output/hmxb_crossmatch.csv`, one row per
HMXB with an X-ray counterpart, and the single place where the cross-match is
to be read. It replaces
the per-mission comparison tables, whose contents are all here.

**Columns.** `ID`; for Gaia, `Gaia_ours`, `Gaia_F23`, `Gaia_SIMBAD`; for
Chandra, `CSC_ours`, `CSC_F23` and `CSC_SIMBAD`; for XMM-Newton, `XMM_ours`,
`XMM_F23` and `XMM_SIMBAD`; and `Notes`. Eleven columns, all of them about
which source in each catalogue is this system. The byte ranges of the three F23
identifier columns are read through one helper, `cut(line, span)`, which takes
them 1-based and inclusive exactly as the ReadMe writes them, so the table and
the ReadMe can be compared line by line. Nothing about the orbit is here: that
is X-09.
Separations and Mahalanobis distances are given inside `Notes` rather than in
columns of their own, so a number appears where the sentence that needs it
is, and only there: the dM2 is quoted for the two sources kept although the
cross-match rejected them, which are the only rows where a reader needs it
to judge the decision.

**The positions behind the separations.** The 4XMM-to-5XMM distance uses the
positions each catalogue tabulates, not positions decoded from the
designations: F23's `RA12deg` and `DE12deg` at bytes 1293-1303 and 1305-1314
of `tablea.dat`, declared "4XMM DR11 Right Ascension/Declination", against
the `RA` and `DEC` of `Output/hmxb_5xmm_crossmatch.csv`, which are the
5XMM-DR15 stacked-source positions, taking the nearest when a system has
more than one accepted source. Decoding the names instead would give only
the quantised position, 0.1 s in right ascension and 1 arcsec in
declination, which is the very rounding that causes the renaming. The
Gaia separation uses `ra` and `dec` from `gaiadr3.gaia_source` for both
identifiers. The observations are deliberately not here: they
live in `../Observations/Chandra/obsids_chandra.txt` and
`../Observations/XMM-Newton/obsids_xmm.txt`, which the reduction scripts
read, and the 5XMM `SRCID` that joins a stacked source to its detections is
in `Output/hmxb_5xmm_crossmatch.csv`.

**What Notes says.** One sentence per disagreement, so that the reason is
readable without going back to the catalogues. One sentence per counterpart
disagreement: a Gaia identifier F23
publishes that does not exist in Gaia DR3; a counterpart F23 does not carry
at all, marked *ID Candidate* — its `GaiaDR3` column is declared `?=0` in
the ReadMe, so a zero there is the null marker and is read as absent, which
41 systems of the catalogue carry; a SIMBAD Gaia source different from ours, with
its separation; an F23 Chandra entry that is not a CSC designation; a source redesignated between
F23's 4XMM DR11 and our 5XMM-DR15, with the distance between the two
positions, which is what says whether it is only a rename; and a source *kept by
decision* with the dM2 at which the cross-match rejected it. A system where
nothing disagrees reads `our counterparts agree with F23 and SIMBAD`.

The note closes with a sentence that is not about a counterpart but about the
system itself: whether F23 and F24 both claim it, and what was decided. Of the
nine systems the two catalogues share (S-03), three reach this table.
`IGR J18483-0311` reads *kept as a HMXB, the literature confirming it, and the
duplicate F24 entry removed instead*; `SAX J1819.3-2525` and
`IGR J21347+4737` read *kept in the HMXB sample as 'Ambiguous'*, naming the
F24 designation in each case. The other six either were removed from the
sample as confirmed non-HMXBs (`1E 1740.7-2942`, `GRS 1758-258`) or have no
X-ray counterpart here and so no row. The flag itself is read from
`Output/hmxb_sample_notes.csv`, where the notebook wrote it; the two
dispositions are decisions and are stated in `build_crossmatch_table.py`
rather than derived.

**Checked.** The table holds the **89** systems that an X-ray catalogue has
something to say about; the 20 with neither a Chandra nor an XMM-Newton
counterpart are dropped, having nothing for this work to compare or reduce.
Of the 89: all 89 carry a Gaia counterpart, 61 a Chandra one, 60 an
XMM-Newton one, 32 both, 3 are also in F24, and 32 agree with F23 and SIMBAD
on every counterpart. The 120 distinct Gaia identifiers of the full sample
were queried against `gaiadr3.gaia_source`; 113 exist, and the seven that do
not are all published values, never ours.

**CXOU designations are excluded from the CSC columns.** SIMBAD lists them
beside the CSC ones for 14 of the 109 systems, but a `CXOU` name comes from
an individual publication and not from the catalogue, so `CSC_SIMBAD` holds
only what SIMBAD gives as a CSC identifier. For four systems —
`GRO J2058+42`, `IGR J13020-6359`, `IGR J13186-6257` and `IGR J18256-1035` —
the CXOU name was the only one SIMBAD had, so their `CSC_SIMBAD` is now
empty, and the count of SIMBAD CSC identifiers falls from 50 to 46. The same
rule is what makes the F23 entry for `AX J1714.1-3912` a non-CSC value in
X-06, although that system is outside the 109.

### X-09 — the orbital parameters F23 reports

`build_orbital_table.py` writes `Output/hmxb_orbital_parameters.csv`, one row
per row of X-08's table, 89 systems and 23 columns. For each of the seven
quantities F23 measures it carries the value, the error and the reference, in
the catalogue's own order: `Mx_Msun`, `Mopt_Msun`, `Porb_d`, `Psuporb_d`,
`ecc`, `Pspin_s`, `RV_kms`. How many systems carry each: `Porb_d` 72,
`Mopt_Msun` 66, `Pspin_s` 57, `ecc` 41, `RV_kms` 28, `Mx_Msun` 12,
`Psuporb_d` 5. The byte ranges come from the F23 ReadMe through the same
`cut(line, span)` helper, 1-based and inclusive as written there.

The catalogued distance is deliberately absent: it is not an orbital
parameter, and the distance this work uses is Bailer-Jones (2021).

These parameters are a separate table from the cross-match because they answer
a separate question. X-08 is about which source in each catalogue is this
system; this one is about what is known of its orbit. The two join on `ID`, and
each stays readable on its own. `download_f23_refs.py` reads its `r_Porb`
column.

**The phase zero point and the longitude of periastron.** Neither is in F23.
The catalogue gives periods without epochs, so no orbital phase can be computed
from it for any system, and it gives no argument of periastron, without which
epochs measured against different reference points cannot be brought onto a
common phase. Both are read out of the literature by hand into
`Input/Ephemerides/porb_ephemerides.tsv` and merged into the table on the
system name as `MJD_T0`, `e_MJD_T0`, `T0_kind`, `r_T0`, `omega_deg`,
`e_omega_deg` and `r_omega`. `r_T0` and `r_omega` are separate columns and
frequently differ, and both frequently differ from `r_Porb`: the paper that
measured the period is often not the one that published an epoch.

Every row carries the sentence or table row it was read from, in a `quote`
column, so that an entry can be checked against its source without re-reading
the paper. Epochs given as JD or HJD are converted to MJD in the file and the
original form is kept in the quote. Where an uncertainty is asymmetric the
larger side is recorded and both appear in the quote.

**Eccentricities read from the literature.** F23 reports an eccentricity for
only 41 of the 72 systems with a period. Where it reports none and a paper
gives one, the value goes into the same file as `ecc_lit`, `e_ecc_lit` and
`r_ecc_lit`, and the merge writes it into `ecc`, `e_ecc` and `r_ecc` **only
where F23 is empty** — F23's own value is never overwritten, and `r_ecc`
therefore always names the source of whatever number is in the column. Three
were added this way: `XTE J1855-026` 0.04 +/- 0.02 (Falanga et al. 2015),
`XTE J0421+560` 0.62 (Barsukova et al. 2006) and `1A 1118-615` 0, the orbit
being consistent with circular (Staubert et al. 2011). The column now stands at
44. Bounds and ranges are not written into the column, because they are not
values: the `< 0.37` of `IGR J18450-0435`, the `< 0.0016` of `Cen X-3` and the
`0.3-0.4` estimate Bozzo et al. (2024) quote for `SAX J1818.6-1703` are stated
in `Notes` instead.

**Searching under F23's other identifiers.** F23 carries, besides the system
name, up to thirteen catalogue identifiers per source — `HEAO`, `UHURU4`,
`ARIEL3`, `IGR`, `ROSAT`, `SAX`, `Swift`, `XMM`, `Chandra`, `2MASS` and others
— and the literature frequently publishes an orbital solution under one of
those rather than under the name F23 uses. Searching them is what found
`H 1417-624`: its own period reference yields 79 characters of readable text,
but the system is `2S 1417-624`, whose full solution (`Tomega`, `omega`, `e`)
is in Raichur & Paul (2010), a paper already on disk for two other systems.
The same route identified `IGR J18450-0435` = `AX J1845.0-0433`,
`3A 0726-260` = `4U 0728-25`, `2E 1145.5-6155` = `2S 1145-619` = `H 1145-619`,
`IGR J16195-4945` = `AX J161929-4945`, `Cep X-4` = `GS 2138+56`,
`1H 1249-637` = `HD 110432` and `RX J0146.9+6121` = `LS I +61 235`. Those
searches were run and the papers read; none of the last six publishes an epoch,
which is why they remain in the list below rather than being absent from it.

**`T0_kind` matters as much as `MJD_T0`.** The epochs are not of one kind:
times of periastron passage, mean longitude 90 degrees (`Tpi/2`), mid-eclipse,
superior and inferior conjunction, the maximum of an outburst, an ellipsoidal
minimum, one spectroscopic phase zero, and, in eight cases, a folding origin
with no physical meaning — the start of a monitoring campaign, used only to
fold a light curve. A phase computed against one of these is not the same phase
as one computed against another, so the column names which, in the paper's own
terms. Converting between them is what the eccentricity and `omega_deg` are for.

**Whose reference is it.** The question only arises for a quantity F23 has a
column for. F23 publishes neither an epoch nor a longitude of periastron for
any system, so `r_T0` and `r_omega` are by construction not F23 columns and
`Notes` names those references plainly, without qualification; a reference
named a second time in the same note reads *the same paper*. Where F23 does
have a column and the value in it is not F23's, the note says so: the
eccentricity sentence reads *the eccentricity is not F23's, which reports none
for this system, but that of <bibcode>*, and the orbital period of
`IGR J18450-0435`, the only period in the table that is not F23's, is explained
at length in its own sentence.

**What `Notes` says, and what it does not.** The epoch, its uncertainty, its
kind, its reference, the longitude of periastron and its reference are columns
of this table, so `Notes` does not repeat them in prose: reading `MJD_T0`,
`T0_kind` and `r_T0` across a row says the same thing more compactly and can be
sorted and filtered. `Notes` carries only what no column can — why a value is
absent, where one came from when the column cannot say it, and what is wrong
with one that is there:

- the withdrawn-reference sentence, for the two systems whose quantities were
  discarded, with what replaces them;
- the misattribution sentence for `RX J2030.5+4751`;
- `no phase zero: <reason>`, for each of the 15 systems that have a period and
  no epoch, naming what the paper does instead;
- `the eccentricity is not F23's, which reports none for this system`, where a
  literature value fills the column, since `r_ecc` alone cannot say whether F23
  had one;
- a per-system caveat where one applies: a bound rather than a value
  (`< 0.37`, `< 0.0016`), an orbit consistent with circular, an omega that
  could not be read from its table.

Twenty-five of the 89 rows carry one; the other 64 are fully described by their
columns. The rule to read alongside them is that an epoch cannot be converted
to another convention without the eccentricity **and** `omega_deg`, so a row
with an eccentric orbit, an epoch and no omega is one whose phase cannot yet be
homogenised.

**The epochs brought to one convention.** `build_orbital_table.py` converts
the published epochs to a common reference, the **periastron passage**, and
writes `MJD_Tper`, `e_MJD_Tper` and `Tper_from`. Periastron is the target
because it is defined for every eccentric orbit, is a single instant shared by
both components, and is already the convention of 19 of the 57 epochs.

Two conversions cover everything convertible. From `Tpi/2`, the mean longitude
is `l = M + omega` and `M` advances linearly, so
`T_per = T0 - (pi/2 - omega) * P / 2pi`, exact and needing neither the
eccentricity nor an iterative solve. From a conjunction, which is fixed in true
anomaly rather than mean anomaly, `nu = pi/2 - omega` for a superior
conjunction and `-pi/2 - omega` for an inferior one, then `E` from `nu` without
iteration and `M = E - e sin E`, and `T_per = T0 - M * P / 2pi`. The offset is
an angle, defined modulo one orbit, so the periastron nearest the published
epoch is taken, `|T_per - T0| <= P/2`.

**A third route, added on 28 September: an outburst maximum of a Be X-ray
binary is a periastron passage.** Fornasini et al. (2023), Sect. 5.1.2 and the
caption of their Fig. 5, state that the Type I outbursts of Be XBs "happen
(quasi-)periodically at or near periastron", fed by the neutron star crossing
the decretion disc, and that "the Type I outburst periodicity, as measured from
X-ray or optical light curves, provides a way of determining the orbital
period" — that is, for these systems the epoch and the period are the same
measurement, and the epoch is a periastron passage by construction. The
conversion is therefore the identity, `T_per = T0`.

Whether it applies is decided by **F23's own `Class` column**, read from bytes
143-154 of `Input/HMXB/tablea.dat` (`Be`, `sg`, `SFXT`, `WR`), and not by any
judgement of ours. Of the 10 epochs of kind `outburst`, five belong to systems
F23 classes `Be` and are adopted: `IGR J01363+6610` (B1Ve), `IGR J11435-6109`
(B0.5Ve), `IGR J19294+1816` (B1Ve), `RX J0812.4-3114` (B0.2IVe) and
`SAX J2239.3+6116` (B0Ve). Each of their papers reports a *recurring* maximum,
not a single event: a sine fit with an explicit `+ n x P` term for
`IGR J01363+6610` and `IGR J11435-6109`, "8 maxima ... which coincide with the
times predicted by a 117 day period" for `IGR J19294+1816`, a maximum with an
assumed 81.3 d period for `RX J0812.4-3114`, and "five outbursts in total ...
with a regular interval time of 262 days" for `SAX J2239.3+6116`, read in
in 't Zand et al. (2000) itself.

The other five are **not** adopted, the quoted statement being about Be XBs and
not about them: `IGR J16465-4507`, `IGR J18462-0223` and `IGR J18483-0311` are
SFXTs and `IGR J19140+0951` and `XTE J1855-026` are supergiants. The outbursts
of eccentric SFXTs are also reported to cluster near periastron — the `status`
of `IGR J16465-4507` says so — but that is a different claim with a different
source and it is not applied here.

**The cost of this route.** Fornasini et al. describe Type I outbursts as
"usually short lived, lasting only for a small fraction of the orbit", and
"(quasi-)periodically at or near periastron" is not "at periastron". A
systematic of order the outburst width therefore rides on these five epochs.
It is **not** included in `e_MJD_Tper`, which carries only the published epoch
error, and the affected rows say so in `Notes`.

An optical maximum, an ellipsoidal minimum, a spectroscopic phase zero, a
folding origin, and an outburst maximum of a system that is not a Be XB have no
fixed relation to the orbit and are not converted.

**The frame of `omega` was checked system by system.** Equation (conj) needs
the compact object's `omega`. Pulse timing gives that; radial velocities of the
companion give the donor's, and the two differ by exactly 180 degrees, so
getting it wrong displaces the periastron by half a cycle. The diagnostic used
is which projected semi-major axis the paper quotes: `a_X sin i` in light
seconds means pulse timing and the compact object's frame, `K1` or
`a1 sin i` of the star means the donor's. Twenty of the 27 are in the compact
frame; seven are in the donor's — `1FGL J1018.6-5856`, `2S 0114+650`,
`4U 1700-377`, `Cir X-1`, `IGR J17544-2619`, `LS 5039`, `LS I+61 303` — and are
turned by 180 degrees first. The classification is a column,
`omega_frame`, and the evidence is in each row's `quote`. `4U 1700-377` is the
one that had to be chased: it has no pulsations, its `omega` reaches us through
Falanga et al. (2015) whose footnote `c` names Hammerschlag-Hensberge et al.
(2003), an optical study, and whose `a_X sin i` column gives the range `48-82`
rather than a measurement — the signature of a system without pulse timing.
The `Tpi/2` conversion needs no frame correction, `Tpi/2` and `omega` being
defined together in the same paper.

**What the footnotes of Falanga et al. (2015) carry.** That compilation is the
source of four of our values, and its footnotes are load-bearing, so they are
recorded here rather than left in the PDF:

- `m` — *"Mid-eclipse time, equivalent to time when mean longitude l equals
  pi/2 for a circular orbit."* This is the published statement of the identity
  the conversion rests on, and the one the implementation checks in the
  `e -> 0` limit.
- `b` (`Vela X-1`) — Bildsten et al. (1997); Raichur & Paul (2010). Both are
  pulse-timing analyses, which is what fixes its `omega` in the compact
  object's frame.
- `g` (`XTE J1855-026`) — Corbet & Mukai (2002), pulse timing, likewise the
  compact object's frame.
- `c` (`4U 1700-377`) — Rubin et al. (1996); Hammerschlag-Hensberge et al.
  (2003). The second is an optical radial-velocity study, which is what puts
  this system's `omega` in the donor's frame.
- `l` — *"The eccentricity for this source has been reported by
  Hammerschlag-Hensberge et al. (2003) to be 0.22(4) and later questioned by
  Clark (2000) and this work."* The source is `4U 1700-377`, and the caveat
  reaches our table twice over: the eccentricity we carry for it, 0.03, is
  F23's from Islam & Paul (2016) and is not the questioned 0.22, but the
  `omega = 49(11)` used to convert its mid-eclipse epoch comes from the same
  optical analysis whose eccentricity was questioned. Its `MJD_Tper` should be
  read with that in mind; the conversion happens to be insensitive to it,
  because at `e = 0.03` the shift is nearly the circular-orbit value.
- `f` — `SAX J1802.7-2017` is the same source as `IGR J18027-2016`, the name
  under which F23 and this table carry it. Falanga's row for it is how its
  period and mid-eclipse epoch can be found in that paper at all.

**One epoch not adopted.** Falanga et al. also publish a mid-eclipse epoch for
`XTE J1855-026`, MJD 51495.25 +/- 0.02. We keep Corbet's outburst maximum,
MJD 50289.1, because it is the epoch belonging to F23's own period reference;
but the Falanga epoch is convertible where the outburst maximum is not, so
adopting it would raise the count on the common convention from 9 conversions
to 10. The choice is recorded here rather than made silently.

**Validated on the one system that publishes both conventions.** Falanga et al.
give `Vela X-1` both `T_ecl = 42611.349(13)` and `Tpi/2 = 42611.1693(43)`, with
`e = 0.0898(12)` and `omega = 152.59(92)`. The formulae predict a separation of
**+0.2204 d**; the published values differ by **+0.1797 d**. The residual,
0.041 d, is 0.45% of the period and **1.5 sigma** once the uncertainty on
`omega` alone, which propagates to 0.023 d, is combined with the errors on the
two epochs. The two are independent fits, one to eclipse times and one to
pulse-timing epochs, so their difference carries the systematics of both and is
not constrained to reproduce the geometric relation exactly. The implementation
also reproduces `Tpi/2 = T_sup.conj.` in the limit `e -> 0` to 2.5e-12 d.

**Result.** 33 systems are on the common convention: 19 published as periastron
passages, 9 converted through the two geometric routes (6 from `Tpi/2`, 2 from
mid-eclipse, 1 from an inferior conjunction) and 5 adopted from an outburst
maximum of a Be XB. It was 28 before the Be route was added. Uncertainties come from 2e4 draws on the quoted errors
of `P`, `omega` and `e`, combined with the error on the published epoch; the
seed is fixed at 20260927 so the numbers are reproducible. The largest are
`1E 1145.1-6141` and `4U 1907+097` at +/- 0.51 d, whose `omega` carry +/- 8 and
+/- 20 degrees. `XTE J0421+560` converts with a shift of exactly zero, its
`omega` of about 270 degrees putting inferior conjunction at periastron.

**Deferred, at the user's direction, until the observation times are known.**
The cycle count from the published epoch to the epoch of interest, whose error
`n * sigma_P` usually dominates everything else and whose `n` is fixed only by
the observations to be folded; and apsidal motion, which makes `omega` a
function of epoch for the few systems where `omega-dot` is measured
(`4U 0115+634`, `4U 1538-522`). Both are properties of the epoch one converts
*to*, not of the conversion, and neither changes the values now in the table.

**Where it stands.** Of the **72** systems with a period, **57** carry a zero
point, **27** a longitude of periastron and **44** an eccentricity. The 15 without one are
`1A 0535+262`, `1H 1249-637`, `2E 1145.5-6155`, `3A 0656-072`, `3A 0726-260`,
`4U 1954+319`, `Cep X-4`, `IGR J11215-5952`, `IGR J16195-4945`,
`IGR J16207-5129`, `RX J0146.9+6121`, `RX J2030.5+4751`, `SAX J1818.6-1703`,
`SGR 0755-2933` and `XTE J1906+090`. In every one of those cases a paper was
read and found to state no epoch, not merely left unexamined: several report a
period from a periodogram and stop there (`IGR J16207-5129`, `Cep X-4`,
`3A 0726-260`); two rest on periods their own authors do not stand behind
(`SGR 0755-2933`, tentative; `XTE J1906+090`, whose two candidate periods the
paper calls of unclear origin); one has no closed orbit at all (`4U 1954+319`:
*the velocities have not yet closed*); and `IGR J11215-5952` folds its
outbursts on 164.6 d against the peak of the 2007 February 9 outburst without
ever writing that date as a number. `1A 0535+262` is the one real gap: its
periastron epoch is universally quoted as MJD 53613.0 from a 2006 contribution
by Finger et al. that has no eprint and could not be identified as a citable
publication, and its F23 period reference, Okazaki & Negueruela (2001), gives
no epoch.

**Nothing was entered from a search-engine summary.** Every value in the file
was read in the paper itself, and three candidates were rejected on that rule:
a phase zero of MJD 53671 for `SAX J1818.6-1703` attributed to Bird et al.
(2009), which that paper does not state — 53671 is the date of a flare in its
Table 1; the `1A 0535+262` epoch above; and an omega for `4U 2206+543`, whose
own table renders in the PDF as `omega (deg) 61 0.2±1` and cannot be read
unambiguously, so that row carries the epoch and no omega.

**A value whose cited source does not contain it.** F23 gives
`RX J2030.5+4751` an orbital period of 46.02 d and an eccentricity of 0.41,
citing Sidoli & Paizis (2018). Table 1 of that paper carries exactly those two
numbers in the row for `EXO 2030+375`, and the paper does not mention
`RX J2030.5+4751` anywhere — the strings `RX J2030` and `2030.5` do not occur
in it. Both values are left in the table and flagged in `Notes` through the
`MISATTRIBUTED` map in `build_orbital_table.py`, not discarded: unlike the
withdrawn preprint, the source is readable and anyone can check what it says.
It also explains why that system has no zero point — there is no orbital
solution of its own for one to come from.

**Papers fetched for this.** `f23_refs/` now also holds nine references that
F23 does not cite, added because they carry an epoch or an omega that F23's own
reference does not:

| bibcode | what it gave |
|---|---|
| `2015A&A...577A.130F` | Falanga et al., ten eclipsing HMXBs: mid-eclipse epochs for Cen X-3 and Vela X-1, omega for Vela X-1, XTE J1855-026 and 4U 1700-377 |
| `2009ApJ...698..514A` | Aragona et al.: periastron epochs and omega for LS 5039 and LS I +61 303 |
| `2010MNRAS.409.1220D` | Drave et al.: the folding origin of XTE J1739-302 |
| `2010MNRAS.406L..75C` | Clark et al.: the outburst-peak ephemeris of IGR J16465-4507 |
| `2010ApJ...709.1374K` | Kubota et al.: the primary-eclipse ephemeris of SS 433 |
| `2011A&A...527A...7S` | Staubert et al.: Tpi/2 for 1A 1118-615 |
| `2000ApJ...530L..33C` | Corbet & Peele: the outburst-maximum epoch of RX J0812.4-3114 |
| `2009MNRAS.393L..11B` | Bird et al., SAX J1818.6-1703: read, states no epoch |
| `2009ApJ...696.2068R` | Romano et al., IGR J11215-5952: read, states no numeric epoch |

The last two are kept although they yielded nothing, so that a reader checking
why those systems have no zero point can see the paper that was read.
`download_f23_refs.py` leaves all nine alone: it only fetches bibcodes that
appear in `r_Porb`.

**A withdrawn source, and what replaces it.** F23 cites
`2015arXiv150301087G` — arXiv:1503.01087, González-Galán (2015),
*Fundamental properties of High-Mass X-ray Binaries*, a thesis posted as a
preprint — for the orbital period, the eccentricity and the systemic radial
velocity of `AX J1841.0-0536` and `IGR J18450-0435`, and for nothing else in
the catalogue. The submission was withdrawn: ADS resolves the bibcode to the
arXiv abstract page, and arXiv answers the PDF request with *withdrawn and is
unavailable*. `download_f23_refs.py` therefore never obtained it, and the two
systems dropped silently out of the zero-point study. A value whose only
source cannot be read cannot be checked by anyone, so all six quantities are
discarded: the script matches on the reference column, not on the system name,
so any other quantity F23 ever attributed to that bibcode would go with them.
`Notes` of this table, and only of this table, records it in words for the two
systems concerned, naming the bibcode as F23 gives it.
The replacement comes from the two papers in `f23_refs/`:

- `AX J1841.0-0536` — Bozzo et al. (2024), MNRAS 528, 863, Sect. 2.1: *"The
  system orbital period remains so far elusive, while the tentative detection
  of a 4.7 s spin period by ASCA (Bamba et al. 2001) has been questioned by
  Bozzo et al. (2011) and never confirmed in other data sets."* Nothing
  replaces the discarded values; the system now has no period, eccentricity or
  systemic velocity in the table. Its `Pspin_s` of 4.7394 s is left in place,
  since F23 attributes it to Bamba et al. (2001) and not to the preprint, but
  `Notes` records that the same passage questions it.
- `IGR J18450-0435` — Goossens et al. (2013), MNRAS 434, 2182, measure
  `Porb = 5.7195 +/- 0.0007 d` for the same system under its other name
  `AX J1845.0-0433`, from a Lomb-Scargle analysis of INTEGRAL/IBIS 18-60 keV
  data with the error from bootstrap resampling. That period is adopted, with
  `r_Porb` set to `2013MNRAS.434.2182G`. It is not a refinement of F23's
  4.73983 d but a different value. The paper bounds the eccentricity at
  `< 0.37` from the Roche-lobe geometry rather than measuring it, and contains
  no radial velocities at all, so no eccentricity or systemic velocity is
  adopted and the bound is stated in `Notes`. Its ephemeris, phase 0 at MJD 52708.43297, is the start
  of the IBIS coverage and carries no physical meaning, so it is a folding
  origin and not a periastron or a conjunction.

Neither paper is cited by F23: `2013MNRAS.434.2182G` appears nowhere in
`tablea.dat`, so the published measurement existed and the catalogue took the
preprint instead. Both papers are in `f23_refs/` as the published articles,
`Goossens2013.pdf` and `Bozzo2024.pdf`, and both are in `references.bib`.
Because `r_Porb` now names the Goossens bibcode, `download_f23_refs.py` would
otherwise fetch the arXiv version of a paper already on disk, so it carries a
small `ALREADY` map from those two bibcodes to the filenames used here. Not checked: whether the withdrawn preprint's numbers appear
in any other publication by the same author.

**Why a redesignation is the commonest note, 39 of 109 systems.** A 4XMM or
5XMM designation encodes the position, quantised at 0.1 s in right ascension
— about 1.5 arcsec — and 1 arcsec in declination. 5XMM re-derives positions
from the stacked detection, so a shift of a few tenths of an arcsecond
renames the source. Measured: over those 39 systems the distance between the
position F23 quotes and the nearest of ours runs from 0.27 to 5.80 arcsec,
median 0.94; 21 agree within 1 arcsec, 36 within 2, and only one lies beyond
3 arcsec. That one is `AX J1700-419`, the single genuine disagreement of
X-06. So 38 of the 39 are renames of the same source, not competing
identifications.

**Checked against the notebook.** The `Gaia_F23` column agrees with the F23
column of `Output/gaia_id_comparison.csv`, which notebook 2 builds
independently, for all 89 systems.

**Not checked.** The table is assembled from the accepted tables of X-01 and
the 5 arcsec match of X-07; it does not rerun either. `Notes` reports
disagreement, not who is right.

## 5. What the X-ray catalogues report

### C-01 — CSC: the fitted column is a total, not a decomposition

- **Claim.** The column returned by the CSC spectral fits is the total
  column, and the model has a single absorbing screen.
- **Source.** `Papers/Evans2024.pdf`, Sect. 3.13 *Spectral Model Fits*,
  and Table 4.
- **How checked.** read. Sect. 3.13 gives the absorbed power-law model as
  `f(E) = e^{-N_H σ_E} A E^{-Γ}` — one exponential, one column. Table 4
  describes `powlaw_nh` as "**Total** neutral Hydrogen column density,
  N_H, of the best fitting absorbed power-law model spectrum to the
  source region aperture"; `bb_nh`, `brems_nh` and `apec_nh` carry the
  same wording for the other three models.
- **Also read from the same section.** The fit threshold is 150 net
  counts in the broad (0.5–7.0 keV) band; the statistic is χ² with data
  variance σ²ᵢ = N_{i,S} + (A_S/A_B)² N_{i,B}; the cross-sections are
  Bałucińska-Church & McCammon and the metal abundances Anders &
  Grevesse, "assumed for all models, unless otherwise noted".

### C-02 — CSC: the Galactic column exists, but only freezes the flux models

- **Claim.** `nh_gal` never enters a fit as a separate component; its
  only role is to freeze the canonical models used to turn aperture
  photometry into an energy flux. It is integrated through the whole
  Galaxy and carries no distance information.
- **Source.** `Papers/Evans2024.pdf`, Sect. 3.14 *Spectral Model Energy
  Fluxes*, and Table 4.
- **How checked.** read. Sect. 3.14: "For all of the canonical spectral
  models, we define the absorption model component total neutral
  hydrogen column density to be fixed with N_H = N_H(Gal), the measured
  Galactic absorption column density from [an H I survey]." The same
  section states that for the canonical models "all of the model
  parameters are frozen except for the normalization". Table 4 describes
  `nh_gal` as "Galactic neutral Hydrogen column density, N_H(Gal), in the
  direction of the source", with no distance entering.
- **Note.** The underlying H I survey is named in `Evans2024.pdf`. It is
  not cited in the manuscript because that paper is not in `Papers/`
  (`Agreement2.md`: never cite a paper that is not in there); the
  manuscript attributes it to `\citet{evans2024}` instead.

### C-03 — CSC: canonical model parameters

- **Claim.** Power law Γ = 2.0; black body kT = 0.75 keV; bremsstrahlung
  kT = 3.5 keV; APEC kT = 6.5 keV.
- **Source.** `Papers/Evans2024.pdf`, Sect. 3.14.
- **How checked.** read. The section derives each value from the
  distribution of CSC 1.1 fits with reduced χ² ≤ 1.25 and then states
  the adopted canonical value.

### C-04 — CSC: energy bands

- **Claim.** ACIS ultrasoft 0.2–0.5, soft 0.5–1.2, medium 1.2–2.0, hard
  2.0–7.0, broad 0.5–7.0 keV; HRC wide 0.1–10 keV.
- **Source.** `Papers/Evans2024.pdf`, Table 1 *CSC Energy Bands*.
- **How checked.** read.

### C-05 — CSC: hardness ratio definition

- **Claim.** H_xy = (F_x − F_y)/(F_x + F_y), with F the **photon** flux in
  the PSF 90% ECF aperture and x always the harder band; changed from
  release 1, where the denominator was the broad-band flux.
- **Source.** `Papers/Evans2024.pdf`, Sect. 3.15 *Spectral Hardness
  Ratios*.
- **How checked.** read.

### C-06 — 5XMM: the fitted column is explicitly not Galactic

- **Claim.** The column of the 5XMM spectral fit is a single total
  column, and the catalogue says so.
- **Source.** `Papers/5XMMdraft.pdf`, Sect. 6 *Spectral fitting*.
- **How checked.** read. The model is `cflux * phabs * powerlw` — one
  `phabs`. The three free parameters are described as "the flux
  (`SPEC_FLUX_PL` in the catalogue, **observed flux not corrected for
  absorption**), the column density (`SPEC_NH_PL`, **not constrained by
  the column density of our Galaxy in the direction of the source**) and
  the spectral slope (`SPEC_GAMMA_PL`)".
- **Also read from the same section.** BXA/UltraNest sampling; a
  log-uniform (Jeffreys) prior on N_H over [0.001, 1000] × 10²² cm⁻², i.e.
  10¹⁹–10²⁵ cm⁻²; a uniform prior on Γ over [1.0, 3.0]; a uniform prior on
  log₁₀(flux) over [−15, −9]; Cash statistic; one pn and one MOS spectrum
  per stacked source, the highest-S/N detection of each.

### C-07 — 5XMM: the detection-level fit is also a single absorber

- **Claim.** `STACK_NH` and `STACK_GAMMA` come from an absorbed power law
  fitted to the five band count rates during source detection, on an ECF
  grid covering 10¹⁹–10²³ cm⁻² and 0 ≤ Γ ≤ 5.
- **Source.** `Papers/5XMMdraft.pdf`, Sect. 4 *Stacking and source
  detection*.
- **How checked.** read. Sect. 6 itself calls this fit "restricted to
  five energy points … and therefore limited".

### C-08 — 5XMM: no separate Galactic column identified in the paper

- **Claim.** Unlike the CSC, 5XMM tabulates no Galactic column.
- **Source.** `Papers/5XMMdraft.pdf`, whole text.
- **How checked.** counted. Every occurrence of `N_H` / `nh` /
  "hydrogen column" / "Galactic column" in the paper was listed and
  inspected. They are: `STACK_NH` (Sect. 4), `SPEC_NH_PL` (Sect. 6), the
  ECF grid limits (Sect. 4), and the `SPEC_FLAG_PL` values 8–11 that flag
  a pegged N_H. None of them is a foreground column. **Not checked:** the
  online column list in the 5XMM-DR15 User Guide was not gone through
  column by column; the claim rests on the catalogue paper.

### C-09 — absorption-corrected fluxes were not identified in the papers

- **Claim.** The words "unabsorbed" and "intrinsic absorption" occur
  nowhere in either catalogue paper, and 5XMM says its flux is "not
  corrected for absorption".
- **Source.** `Papers/Evans2024.pdf`, `Papers/5XMMdraft.pdf`.
- **How checked.** counted, case-insensitively, over the full extracted
  text of both papers:

  | string | Evans2024 | 5XMMdraft |
  |---|---|---|
  | `unabsorb` | 0 | 0 |
  | `de-absorb` | 0 | 0 |
  | `intrinsic absorption` | 0 | 0 |
  | `corrected for absorption` | 0 | 1 (the negation, Sect. 6) |
  | `luminos` | 0 | 3 |
  | `absorb` | 95 | 8 |
  | `flux` | 254 | 28 |

  The last two rows are the control: the text of both papers is being
  read correctly, so the zeros above are real absences and not a broken
  extraction.
- **Not checked.** This is a statement about the catalogue *papers*. It
  is not a statement about every column of the released catalogue
  products.

### C-10 — 5XMM: hardness ratio definition

- **Claim.** HR_i = (r_{i+1} − r_i)/(r_{i+1} + r_i) between the **count
  rates** r in adjacent bands i and i+1.
- **Source.** `Papers/Traulsen2019.pdf`, the column description table
  (Cols. 140–171, `EP_HRi`): "Equivalent all-EPIC hardness ratios
  (r_{i+1} − r_i)/(r_{i+1} + r_i) between the count rates r in energy
  bands i and i + 1." Carried over into `Papers/Traulsen2020.pdf`
  (Cols. 142–149) and used unchanged in 5XMM-DR15.
- **How checked.** read.

### C-11 — 5XMM: bands and catalogue size

- **Claim.** Bands 1–5 = 0.2–0.5, 0.5–1.0, 1.0–2.0, 2.0–4.5,
  4.5–12.0 keV; EP_8 = 0.2–12.0 keV. 818 656 unique sources from
  2 578 752 detections.
- **Source.** `Papers/5XMMdraft.pdf`, Sect. 3 (bands) and the abstract
  (counts).
- **How checked.** read.

### C-12 — 5XMM runs a classifier that consumes these quantities

- **Claim.** 5XMM classifies its X-ray sources with an adapted CLAXBOI
  (Tranin et al. 2022) using hardness ratios, spectral-fit parameters and
  X-ray luminosity as features, and labels 26 100 sources as Galactic
  X-ray binaries.
- **Source.** `Papers/5XMMdraft.pdf`, Sect. 9 *Classification*.
- **How checked.** read.

---

### C-13 — how each catalogue builds a source-level value

**How checked.** read, on 26 September, and written into Sects. 3.1-3.3 of
the manuscript.

- **CSC.** Bayesian Blocks over the contributing observations, the fitness
  of a block being the product of the per-observation flux MPDFs, per-band
  blockings intersected, run once ordered by epoch and once by flux. The
  published source value comes from the flux-ordered block with the longest
  exposure; the other observations do not enter it. Spectral fits are made
  per block with the detections fitted simultaneously.
  `../Papers/Evans2024.pdf` Sects. 3.11.2, 3.13 and the Table 5 column
  descriptions.
- **5XMM.** Observations overlapping by at least 1 arcmin are stacked and
  searched simultaneously, assuming "that the flux of each source remains
  constant over all exposures" and an absorbed power law; the free
  parameters are position, one mean flux, N_H and Gamma, reported in the
  `STACK_` columns. After detection the assumption is dropped and forced
  photometry gives per-image rates. The fitted spectra combine nothing: one
  pn and one MOS detection, the highest signal-to-noise of each.
  `../Papers/5XMMdraft.pdf` Sects. 4 and 6.

**Consequence recorded in the manuscript.** Neither source-level value is
defined for a variable source, so this work stays at the detection level.
The argument given is that the interstellar column at a fixed distance is
time-invariant while the intrinsic one is not, so the epoch-to-epoch
variation is what separates the two absorbers and aggregation removes it.

**Not used.** No Bayesian Blocks analysis is performed anywhere in this
work. Stage 1 refits per-detection products only: CSC observation-level
fits and 5XMM `SPEC_` fits.

**Pending.** A figure comparing the per-detection values of the
best-sampled sources with the single number each catalogue publishes;
candidates are `2CXO J130247.6-635008` with 21 CSC fits and 5XMM
`3009282010100002` with 20 detections carrying spectra.

## 6. Observations and their reduction

### X-10 — the per-system observation tables, and the orbital phase of each observation

`build_image_tables.py` lays out `../Observations/Images/`, one directory per
system of the working sample that has data on disk, named after the identifier
with runs of whitespace replaced by a single underscore (`LS I+61 303` becomes
`LS_I+61_303`); the identifier itself is the first line of the table inside, so
nothing is lost. Each directory holds `Chandra_images/`, `XMM-Newton_images/`
or both, according to which missions observed the system, and
`observations.txt`.

**The table.** One row per observation, with `ChandraObsID`,
`XMM-NewtonObsID`, `MJD` and `orbital_phase`. An observation belongs to one
mission, so exactly one identifier column is filled per row. Rows are ordered
by MJD.

**Which observation is whose.** Chandra: the kept rows of
`Output/hmxb_csc_match_5arcsec.csv` give the 2CXO names of each system, and the
5 arcsec CSCview pull gives the observations each of those names was detected
in. XMM-Newton: `Output/hmxb_5xmm_crossmatch.csv` gives the `SRCID` of each
system and the 5XMM-DR15 detection table gives the observations of that
`SRCID`.

**MJD** is the mid-point of the observation, `MJDREF + (TSTART + TSTOP)/2/86400`
read from the reprocessed event files themselves and not from a catalogue. For
XMM-Newton the span is taken over all imaging event files of the observation,
since one observation has several exposures and cameras.

**Checked.** 89 systems, 61 with Chandra data and 60 with XMM-Newton. 466 rows,
containing all **193** Chandra and all **192** XMM-Newton observations of the
lists in X-05 and nothing else. The row count exceeds the sum because 81
Chandra observations contain more than one of our systems and so appear in more
than one table.

**The orbital phase, and the cycle-count term that X-09 deferred.** The phase
is counted from the periastron passage, `phase = ((MJD - MJD_Tper) / Porb) mod
1`, using the epochs brought to the common convention in X-09. Where a system
has no such epoch the column is empty: a phase measured from an outburst
maximum or an arbitrary folding origin is not the same quantity and is not
written as though it were.

Now that the observation times are known, the term X-09 left aside can be
evaluated. Carrying an epoch across the `n` cycles that separate it from an
observation costs `n * sigma_Porb`, and it dominates: the phase uncertainty
per observation is `hypot(sigma_Tper, n * sigma_Porb) / Porb`. The phase is
written when that is below **0.5 cycles** and left empty above it. Phase is a
circular quantity, so the measure that matters is the concentration of the
wrapped normal distribution, `R = exp(-2 pi^2 sigma^2)`: it is 0.82 at
`sigma = 0.1`, 0.45 at 0.2, 0.17 at 0.3 and **0.007 at 0.5**, where the
distribution is uniform to under one per cent and the value excludes no phase.
A `+/- 1 sigma` interval of half a cycle also spans the whole cycle, which is
the same statement made crudely. The threshold is a choice, and the result is
insensitive to it: of the 145 observations belonging to a system with an epoch,
129 survive a cut at 0.5, 126 at 0.3, 121 at 0.2 and 113 at 0.1, while the two
systems rejected — at 6.7 and 1.8 cycles — lie far outside any of them. Two systems lose
their phase entirely on this test:

- `SAX J0635.2+0533` — `Porb = 11.2 +/- 0.5 d`, a 4.5 per cent period, and its
  observations sit 151 cycles from the epoch: **6.7 cycles** of uncertainty.
- `Cir X-1` — `Porb = 16.68 +/- 0.15 d` at 195 cycles: **1.8 cycles**.

Four more are usable but should be read with care, their uncertainty being a
large fraction of a cycle: `V 0332+53` (0.31), `GS 0834-430` (0.27),
`XTE J0421+560` (0.20) and `IGR J17544-2619` (0.15). At the other end
`PSR B1259-63`, `H 1417-624`, `EXO 2030+375` and `4U 1901+03` are good to a few
thousandths. Each table states its own number in the header, with the cycle
count it was computed at.

Of the 466 rows, **134** carry a phase: the systems with an epoch on the common
convention, minus the observations where the cycle count has destroyed it. It
was 129 before the Be outburst route of X-09 was added, which put five more
systems on the convention and four of them into these tables; the fifth,
`RX J0812.4-3114`, loses its single observation to a cycle count of 65 on a
period of 80.39 +/- 3.0 d, 3.7 per cent, which is 2.4 cycles of phase
uncertainty.

**Not done.** Apsidal motion, the third correction, is still not applied; it
matters for `4U 0115+634` and `4U 1538-522`, whose phase uncertainties here are
0.004 and 0.042 cycles, so it is not the limiting term for either.

**A limitation of these lists: they are not pointings.** The observations of a
system are those in which its matched CSC or 5XMM source has a record, and a
catalogue record exists for an observation that merely *covers* the position,
not only for one aimed at it. Nothing in `build_image_tables.py` filters on
where the source falls in the field. Measured from `RA_NOM`/`DEC_NOM` of the
image headers over all 242 (system, Chandra observation) pairs:

| off-axis | pairs |
|---|---|
| 0-1 arcmin | 85 |
| 1-3 | 53 |
| 3-6 | 29 |
| 6-10 | 41 |
| beyond 10 | 34 |

So **75 pairs sit beyond 6 arcmin and 34 beyond 10**, where the Chandra PSF has
degraded past any use for this. They concentrate in a few systems whose
positions fall in fields observed for something else: the SNR G21.5-0.9
pointings supply most of them, and the three HRC observations of
\object{GS 0834-430} are pointings at the Vela remnant with the system 14
arcmin off-axis. That is why nothing is visible at the source in those tiles:
counts in a 10 arcsec box at the position are 134, 318 and 357 against expected
backgrounds of 151, 294 and 334, which is background and nothing else, while
its two ACIS observations, which *are* pointings at the system at 0.05 and 0.3
arcmin, give 33 against 5.3 and 4 against 0.7.

**No off-axis cut is applied, and none is intended.** These observations are
part of the sample and are studied like the rest, at the user's direction. What
a large off-axis angle changes is the PSF, the effective area and the
background under the source, so it belongs in the extraction and in the error
budget of each measurement, never in a criterion for discarding observations.
The angle itself is worth having as information — it is not yet a column of
`observations.txt` — but as a quantity to fold into the analysis, not as a
threshold.

### X-11 — the image and source-detection scripts

`../Observations/Images/chandra_images.sh` and
`../Observations/Images/xmm-newton_images.sh` fill the directories X-10 lays
out. They are adapted from the `# Images` and source-detection sections of
`Extra/Scripts/Chandra.sh` and `Extra/Scripts/EPIC.sh`. **Both have now been
run to completion** — see *What the run produced* at the end of this entry for
what is on disk and what failed.

**What changed from the originals.** The originals walked a hard-coded pair of
observations, made one image in a single band, and ran the detection
interactively behind a `y/n` prompt. These read the observations from the
per-system tables of X-10, write into each system's `Chandra_images/` or
`XMM-Newton_images/`, make one broad-band image, run twelve jobs at a time,
give each job its own log, and skip any product already present.

**One broad band per observation, since 28 September.** An ACIS observation is
imaged in `csc-b`, 0.5-7 keV, the band CSC 2.1 itself detects in; an HRC
observation in `csc-w`, its 0.1-10 keV wide band, unfiltered, HRC having no
energy column; and each EPIC camera and exposure in `EP8`, 0.2-12 keV, the full
band of 5XMM-DR15, which lies inside the 200-12000 eV the events were already
filtered to.

**What this replaced, and why the sub-bands were deleted.** Until then a
Chandra observation was imaged in nine bands — the five ACIS bands of CSC 2.1
plus the Watson et al. (2009) EPIC bands that ACIS can supply, with band 5
truncated at 7 keV as `xmm-EP5alt` and 0.2-0.5 keV written once as
`csc-u_xmm-EP1`, the two catalogues defining it identically — and an
XMM-Newton exposure in all six EPIC bands. Nothing in the analysis read them:
the spectral work is per-detection fits to spectra, not to images, and the
mosaics show the broad band only. At the user's direction the 1480 Chandra and
3260 XMM-Newton sub-band files were deleted on 28 September, 190 GB, taking
`Observations/Images/` from 418 to 230 GB, and both scripts now make the broad
band alone. They are reproducible from the same scripts if a sub-band is ever
wanted: the band definitions are one line each.

**Where the counts per band are, for each mission.** They are not symmetric,
and it matters for any comparison built on them:

- **Chandra** — the 5 arcsec CSCview pull already carries them per detection
  per observation: `src_cnts_aper_{u,s,m,h,b}` are the net, aperture-corrected
  source counts in each ACIS band, with `cnts_aper_*` the gross counts in the
  aperture and `cnts_aperbkg_*` the background in it. Nothing has to be
  measured to have the catalogue's own numbers.
- **XMM-Newton** — 5XMM-DR15 publishes **no per-band counts**. It has
  `EP_CTS` and the per-camera `PN_CTS`, `M1_CTS`, `M2_CTS`, all for the total
  band only. Per band it gives rates, `PN_1_RATE` ... `M2_5_RATE`, and fluxes,
  `EP_1_FLUX` ... `EP_5_FLUX`. Counts per band therefore have to be
  reconstructed as rate times the camera's exposure, `PN_ONTIME`, `M1_ONTIME`,
  `M2_ONTIME`, and the result is not a raw count: the catalogue rates are
  vignetting- and PSF-corrected, so it is an equivalent on-axis count.

To measure counts in our own apertures instead, which is what makes the two
missions comparable in the same region rather than in each catalogue's own:
`dmextract` on the filtered events per band with a source and a background
region for Chandra, aperture-corrected with `arfcorr` or done whole by
`srcflux`; `eregionanalyse` per band and camera for XMM-Newton, which returns
the source counts, the background and the encircled-energy correction. A plain
pixel sum over one of the images made here is **not** that number: it is gross,
un-aperture-corrected and un-background-subtracted. And counts alone are not
comparable between the missions even in an identical band, the effective areas
and exposures differing; the comparable quantity is a flux, which needs the
response and belongs to Stage 2.

**Source detection.** Chandra: `mkpsfmap` then `wavdetect`, once per
observation, on the `csc-b` (0.5-7 keV) image, which is the band CSC 2.1 itself
detects in; running it in more bands would multiply the cost by as many and
would not be what the catalogue did. XMM-Newton: `edetect_chain` once per
observation over the three cameras together, on their `EP8` (0.2-12 keV)
images, with the attitude file of that observation. Both write the resulting
list as FITS and as a CSV anyone can read without CIAO or SAS.

**Two passes, and the scratch directory.** Each script runs the images first and
the detection second, as the originals do through two separate prompts, because
the two cost very different things and need different job counts: 12 observations
imaged at a time, 4 detected. The first run of `chandra_images.sh` filled `/tmp`
and every concurrent `wavdetect` then failed with *No space left on device*.
The cause: CIAO writes its working files in `ASCDS_WORK_PATH`, which `ciaoinit`
sets to `/tmp`, and `/tmp` on this machine is a **32 GB tmpfs**; `wavdetect`
works in double precision, so on an unblocked 8192 x 8192 ACIS field each of its
five normalised-background intermediates is 537 MB, 2.7 GB per observation, and
twelve at once is 32 GB. Three changes follow:

- `ASCDS_WORK_PATH` and `SAS_TMPDIR` are set per job to
  `Observations/Images/scratch/<obsid>_{img,det}`, on the project disk, and
  removed when the job ends. Nothing is written to `/tmp` any more.
- each script refuses to start unless `min_free_gb` (60) is free on that
  filesystem, so a full disk is a refusal up front rather than a discovery
  halfway through.
- the Chandra detection runs on a copy of the `csc-b` image blocked by
  `detect_block` (4), which divides both the intermediates and the three
  required outputs by sixteen, to 34 MB and 17 MB. The cost is resolution,
  1.97 arcsec per pixel against the 0.492 of the unblocked image, coarser than
  the on-axis PSF. That is the right trade here: this detection is for finding
  and marking the sources in the field, while the CSC positions this work
  actually uses come from the catalogue. `detect_block=1` restores full
  resolution at 2.7 GB of scratch and 800 MB of output per observation. The
  images the mosaics use are unblocked either way.

**HRC observations get one band, not nine.** 36 of the 194 reprocessed Chandra
observations are HRC, and an HRC event list has no `energy` column at all, only
`pha` and `pi`; the first run failed on all 36 with *Could not find identifier
energy*. HRC also has essentially no energy resolution, which is why CSC 2.1
gives HRC detections one band and one only, the wide band 0.1-10 keV. So an
HRC observation now gets a single image, `csc-w`, made **without** an energy
filter: the band is the whole HRC bandpass, and selecting on `pi` would imply a
precision on the energies that HRC does not have. Detection runs on that image
instead of on `csc-b`. Three systems have Chandra data that is HRC and nothing
else — `4U 1954+319`, `IGR J00370+6122` and `IGR J18483-0311` — and without
this they would have had no Chandra image at all. This also means the HRC wide
band, which Table of bands in the manuscript said was not used, *is* used, for
those 36 observations and only for them.

**What the failed first run left behind.** The 216 *Error in parameter file*
messages in the logs were the same ENOSPC in another guise: CIAO also puts its
per-process parameter directory in `/tmp`
(`cp: error writing '/tmp/71730_param/wavdetect.par': No space left on
device`), which `ASCDS_WORK_PATH` fixes along with the rest. Of the 596 files
the run had produced, two were truncated by the kill,
`1722_xmm-EP4_img.fits` and `1723_csc-b_img.fits`, both missing their END card;
they were deleted so the rerun remakes them, and the other 594 open cleanly. No
source list was completed, so nothing has to be checked there. The skip logic
honours whatever is on disk, so this check has to be repeated after any
interrupted run: `fits.open` on every product, delete what fails.

**Parallelism.** The unit of work is a (system, observation) pair, twelve at a
time, each in a subshell with its own working directory — and, for CIAO, its
own `PFILES`, since concurrent tools writing the same parameter file corrupt
each other, and for SAS its own `SAS_CCF` and `SAS_ODF`. 81 observations
contain more than one of our systems, so there are 242 Chandra and 224
XMM-Newton jobs behind 193 and 192 distinct observations.

**The duplication that caused, and how it is resolved.** Writing every product
once per system made 74.8 GB of redundant copies, a third of everything under
`Observations/Images/`: 441 Chandra files over 49 observations and 674
XMM-Newton files over 32. Almost all of it was HRC, whose 2.0 GB `csc-w` images
alone accounted for 42 GB; the SNR G21.5-0.9 pointings, which fall in both
`ATO J278.3657-10.5901` and `SNR 021.5-00.9`, carried 3.25 GB each.

**They were the same products.** The detection does *not* work on a region
around the system: `wavdetect` and `edetect_chain` run on the whole field, so
two copies of a source list are the same source list. Checked on obsid 144,
whose two copies have byte-identical `source_list.csv` and `src.reg` and whose
FITS files differ only in the `DATE` and `CHECKSUM` cards recording when each
was written — every data unit compares equal.

**Nothing was deleted.** `work/dedup_image_products.py` compares the data units
of every duplicate group and, where they agree, replaces the copies with hard
links to one inode. Every path still exists and every system directory still
lists exactly what it did; only the bytes behind them are shared. Run with
`--apply` on 28 September: **1115 files linked, 74.8 GB reclaimed**, no group
left where the data disagree, and `Observations/Images/` fell from 229 to
155 GB with the file counts unchanged (484 Chandra images, 652 EPIC images, 242
and 224 source lists). NaN has to be compared as equal to itself for this, or
every float column of a detection table reads as different.

**And it does not recur.** Both scripts now look for an existing copy of a
product under another system's directory before making it, and hard-link it
when they find one; the detection pass does the same for the whole product set
of an observation, which also saves 49 redundant `wavdetect` runs and 32
`edetect_chain` runs. Two jobs starting the same observation at the same moment
can still both make it, so the dedup script remains the sweep to run afterwards.

**Checked.** `bash -n` passes on both. The `awk` that
reads the observation tables returns the right identifiers from a real one.
The globs match real files on disk: `*repro_evt2.fits` under
`Chandra_repro/<obsid>` and `flag-pattern_*<camera>*ImagingEvts.fits` under
`XMM-Newton_repro/<obsid>`. The embedded python that converts a source list to
CSV runs on a real FITS table. Both have since been run to completion and the
products counted, which is what *What the run produced* below records; the one
thing that only running could show, that the CIAO and SAS tools accept these
parameter combinations, it did show, with the single exception of the stale PSF
maps described there.

**What the run produced**, counted on disk on 28 September. 421 GB in
`Observations/Images/`.

- **Chandra images: complete.** One image per (system, observation) pair over
  242 pairs, 185 ACIS in `csc-b` and 57 HRC in `csc-w`, each with the blocked
  copy its detection runs on. Walking the per-system tables and asking for each
  expected file by name returns **no missing image**, HRC included. Before the
  sub-bands were deleted this was 1722 files.
- **The HRC images are real and correct.** All 36 HRC observations took the
  unfiltered `csc-w` branch — the branch says so in its own log — and the 57
  files it wrote (36 observations, some shared between two systems) all open,
  all carry `INSTRUME = HRC`, and all are 32768 x 32768 16-bit images with
  `CDELT1 = 3.6611e-05` deg, that is **0.1318 arcsec per pixel**. Each is
  2.1 GB and they are 115 GB of the 421. No ACIS-band image exists for any HRC
  observation, so nothing was written before the branch was added.
- **XMM-Newton: complete.** 652 `EP8` images, one per camera and exposure, and
  224 `emllist` source lists over 224 pairs, that is every one. Before the
  sub-bands were deleted this was 3912 images.
- **Chandra detection: 242 of 242**, after a rerun. It was 228 at first, and
  the 14 that failed had one cause.
  Those 14 are ACIS observations that had been detected once before
  `detect_block` existed, when the PSF map was built from the unblocked
  8192 x 8192 image. `mkpsfmap` runs with `clobber=no`, so on the rerun it
  refused to replace that stale map, said so, and returned a passing status
  anyway; `wavdetect` then stopped with `psf image size (8192, 8192) does not
  match input image size (2048, 2048)`, and the failure only surfaced at the
  end as a Python traceback about a missing source list. The 8192 x 8192 PSF
  maps of the HRC observations are *not* stale: an HRC image blocked by 4 is
  8192 x 8192, so for those the size is the right one.

**The three changes this caused in `chandra_images.sh`.** The detection pass
now (i) deletes the products of an attempt that did not finish before trying
again — it is reached only when the source list is absent, and the source list
is written last, so anything else present is from a failed attempt; (ii) stops
with a message when `mkpsfmap` writes no PSF map, instead of carrying on;
(iii) stops with a message when `wavdetect` writes no source list, instead of
handing the absent file to Python. Rerunning the script is therefore enough to
recover the 14, and everything already made is still skipped. **Checked** on
`8242` (1A 0535+262) before the rerun: the stale map removed, `mkpsfmap` remade
it at 2048 x 2048 and `wavdetect` completed. The rerun has since been done and
the products counted on disk: **242 source lists over 242 pairs**, a `.csv`
beside every one, no PSF map whose size disagrees with the image its detection
ran on, and no scratch or parameter directory left behind.

### X-12 — the phase mosaic, and what it is showing

`../Observations/Images/ds9_phase_mosaic.sh` opens one SAOImage ds9 per system
and tiles every observation of that system in order of orbital phase, Chandra
and XMM-Newton in the same window, each tile centred on the source. It is a
looking tool, not a measurement: nothing it draws enters a number anywhere, and
the decisions below are about legibility, not about photometry. They are
recorded because a picture that is read wrongly is worse than no picture.

`work/build_mosaic_list.py` writes the two lists of what can be drawn, each
ordered by how many tiles the mosaic will hold, with the directory name in the
form the script takes as its argument:

- `Output/mosaic_systems.txt` — the **31 systems** with at least one phased
  observation, which are the phase mosaics, 134 tiles in all.
- `Output/mosaic_systems_no_phase.txt` — the **58 systems** with observations
  and no phase on any of them, **326 observations**, of which **31 systems have
  more than one** and are worth a mosaic even so. Its `why` column says where
  each one stops: no orbital period; a period and no epoch published; an epoch
  of a kind that converts, `Tpi/2`, mid-eclipse or a conjunction, but with no
  longitude of periastron ever published; an epoch that bears no fixed relation
  to the orbit at all, an outburst maximum, a folding origin, an optical
  maximum, an ellipsoidal minimum, a spectroscopic phase zero; or a phase lost
  to the period error. That distinction is worth keeping: a missing $\omega$
  can be supplied by one paper, while a folding origin never can. An outburst
  maximum is no longer in the second group by default — for a Be XB it *is* the
  periastron passage (X-09), and only the outburst epochs of the supergiant and
  SFXT systems remain unconvertible. Checked against
  `Output/hmxb_orbital_parameters.csv`: every unconverted epoch of a
  convertible kind is missing $\omega$ and nothing else, and the one system
  with an $\omega$ and no conversion, XTE J1855-026, is blocked by the kind of
  its epoch instead.

**Time order for the systems with no phase.** `--by-mjd` tiles every
observation by the mid-point of its exposure and writes the MJD on each tile in
place of the phase, the label always being the quantity the mosaic is ordered
by. A system with no phased observation falls back to it on its own, with a
message, rather than refusing to draw observations it does have. In phase order
the unphased rows are still skipped, since they cannot be placed in that
sequence. Checked on SS 433: 14 observations, 6 Chandra and 8 XMM-Newton, in
time order from MJD 52731.66531, five across in three rows.

**Which observations appear, and which do not.** The script reads
`Observations/Images/<system>/observations.txt` (X-10) and takes only the rows
that carry an orbital phase; a row without one is skipped and counted in the
line it prints at the end. So the mosaic shows the observations whose phase
survives X-09 and the half-cycle test of X-10, not every observation of the
system. Cir X-1 is the visible case: 2 of its 8 observations are drawn.

**The position the tiles are centred on** is our adopted Gaia counterpart, from
`Output/hmxb_sample_positions.vot` — the best position for the *system*, not an
X-ray detection. The red circle of 10 arcsec radius marks it. That circle is
deliberately larger than any positional error involved; it is there to say
*here is where the star is*, and, since its radius is known, to serve as the
scale of the tile. A labelled 10 arcsec scale bar did the same job and is
commented out rather than deleted, the two being the same 10 arcsec twice.

**Not every tile is a pointing at the system.** The observation lists include
observations that merely cover the position — 75 of the 242 Chandra pairs put
the source beyond 6 arcmin off-axis and 34 beyond 10, the numbers being in
X-10. An empty tile therefore means *not detected here*, which is not the same
as *not detected*, and for the far off-axis ones it is very often a statement
about the pointing rather than about the source. They are kept and studied all
the same; X-10 says why the off-axis angle is not a cut.

**The band of each tile.** Chandra `csc-b`, 0.5-7 keV, the band CSC 2.1 itself
detects in; XMM-Newton `EP8`, 0.2-12 keV, the full EPIC band. They are the
broad band of each mission and not the same energies, which matters when
comparing a Chandra tile with an XMM-Newton one by eye. An HRC observation has
no energy column at all and shows `csc-w`, the 0.1-10 keV wide band (X-11).

**Pan and zoom, never crop.** ds9's `crop` does not frame an image, it discards
the rest of it from the display: a cropped frame cannot be panned or zoomed
outside the box, and ds9 recomputes the scale limits over what is left, which
on a sparse EPIC image flattens the tile to one colour. Panning to the source
and setting a zoom leaves the whole image loaded, so the mosaic opens on the
field and can then be zoomed out past the edges of the detector and back.

**Three plate scales, three zooms, one field.** An unblocked ACIS image is
0.492 arcsec per pixel, our EPIC images are binned to 4 arcsec (X-11), and an
HRC image is 0.1318 arcsec, read off `CDELT1` of the images themselves. Given
one zoom these would show fields differing by factors of 3.7 and 8.1.
`zoom_chandra` is the only knob; `zoom_xmm` and the HRC zoom are derived from
it, so that **an HRC tile covers the same sky as a Chandra tile** and an
XMM-Newton tile covers `field_ratio_xmm` times as much. The XMM-Newton field is
deliberately *not* equal: the EPIC PSF is about six arcsec against half an
arcsec on axis for ACIS, so at equal field the EPIC tile is a handful of
blown-up squares. `field_ratio_xmm = 3` is the value in use and is of the order
of the ratio of the two point spread functions. At the default
`zoom_chandra = 3` that is 0.164 arcsec per screen pixel on Chandra and 0.492
on XMM-Newton.

**HRC is shown blocked.** At full resolution an HRC image is 32768 x 32768 —
a thousand million pixels, 2.1 GB. `chandra_images.sh` already writes the
blocked copy its detection runs on, `csc-w_block4`, 8192 x 8192 and 129 MB,
whose 4 x 0.1318 = 0.527 arcsec pixel is within seven per cent of the ACIS one,
so it is what belongs beside an ACIS tile on its own merits. The block factor
is read from the file name and the zoom follows it. This was also a hard
failure and not only a preference: with smoothing on, three full-resolution HRC
frames put ds9 at 25 GB of memory and 98 per cent of a core for seven minutes
with nothing drawn. Smoothing is now never applied to a full-resolution HRC
image, whatever the knob says.

**XMM-Newton as one RGB frame.** An observation with all three cameras is drawn
as a single RGB frame, EMOS1 red, EMOS2 green, EPN blue, rather than as three
tiles. **161 of the 192** observations have all three. The other **31 carry no
pn** and fall back to one grey frame from the first camera they do have, EPN,
EMOS1, EMOS2 in that order; it is all three channels or none, because a
two-channel RGB puts a colour on the source that means *this camera was
missing*. **Known and not yet decided:** 8 observations have all three cameras
but more than one exposure of one of them — `0112430103`, `0137150301`,
`0139760101`, `0206380701`, `0672050201`, `0728371001`, `0729560601`,
`0823990901` — and the script takes the first exposure alphabetically, which
picks `EMOS1_U002` for `0604700101` where other frames use `S001`.

**The stretch, and why tiles must not be compared by brightness.** `zscale` is
built for sky-background-dominated optical frames and collapses on sparse
Poisson counts, leaving the tile flat. Each frame therefore uses a square-root
stretch to the 99.5th percentile **of its own data**, `scale scope local`, and
`lock scalelimits` and `lock colorbar` are both off so that the contrast of one
frame can be adjusted without touching the others. The consequence has to be
stated plainly: **the tiles are not on a common intensity scale**, the images
are raw counts with no exposure or vignetting correction, and the exposure
times differ between observations, so a source looking brighter in one tile
than in another says nothing about flux. The colour bar under each tile carries
its own limits and is the only honest reading of its levels. RGB channels use
`minmax` instead of the percentile, because a high percentile on smoothed
sparse data sits just under the peak, all three channels saturate together and
the tile washes out to white.

**Two colormaps.** HRC is given a colormap of its own, `cmap_hrc`, distinct
from the `cmap` that Chandra ACIS and any single-camera XMM-Newton frame use.
An HRC field holds 1.1 to 2.6 million counts where an ACIS image of the same
sky holds a few thousand, so an HRC tile is a full frame of background speckle
rather than a dark field with a source in it, and a different colormap says at
a glance that it is a different kind of picture. The choice is a preference and
the two knobs are meant to be changed; the header the script prints reports
whatever is in force. RGB frames take no colormap at all.

**Smoothing.** EPIC images binned to 4 arcsec have very few counts per pixel, so
they get a Gaussian of radius 2 to make a source read as a source rather than as
three isolated pixels. Chandra is not smoothed, its PSF being narrower than a
display pixel, and blocked HRC is not smoothed either, blocking by 4 having
already summed sixteen pixels. Smoothing is a display setting a new frame
**inherits**, so it is turned off explicitly on every frame that should not have
it; without that, one smoothed XMM-Newton frame smooths every frame after it.

**The coordinate grid.** Type `analysis`, with axes and numerics drawn inside
the image. Type `publication` puts them in a margin that a tile in a mosaic does
not have, so the numerics are simply dropped and the grid says nothing. The
gaps adapt to the zoom either way — measured on an ACIS field, 1 s of RA at
zoom 4, 0.5 s at zoom 8, 0.1 s at zoom 64 — and `analysis` is what makes that
visible. The cost is that with no margin ds9 never draws the frame title, which
is why the phase is written as a text region instead.

**The phase label.** The number above each source is the orbital phase of that
observation, nothing else. It is placed 70 **screen pixels** above the source,
not a fixed angle: a region is placed on the sky, so a fixed angular offset
grows on screen with the zoom, and at zoom 16 on Chandra the 15 arcsec used
before was 487 pixels and the label had left the tile. In pixels it sits the
same distance from the centre of every tile at every zoom.

**The window is the picture.** ds9 has no headless mode, and `saveimage`
photographs the window as it stands on screen: a small window gives a small
picture whatever the frames contain. The window is therefore sized to one
monitor, less `win_margin` for the decorations. The desktop here spans two
monitors, 5952 pixels in all, so the size is taken from the monitor `xrandr`
marks primary and not from the whole desktop. The saved PNG is 3072 x 1213.

**The tiling.** Every possible number of rows is costed and the cheapest wins.
The cost of a layout is how far one tile is from square, `|log((w/c)/(h/r))|`,
plus `empty_cell_cost` for each cell left empty in the last row. Neither term
works alone: squareness by itself puts five frames in one row of five, which
leaves no cell empty and no row uneven, and evenness by itself does the same.
Together they give 3 and 2. Checked over 1 to 20 frames and at 25, 30 and 41 in
a 3012 x 1668 window: no layout leaves the last row more than two short of a
full one — four come out 2 and 2, five 3 and 2, six 3 and 3, seventeen six
across in three rows, and the 41 of PSR B1259-63 seven across in six.
`tile_cols` forces a number and ignores both terms.

**Frames are not locked.** `lock frame wcs` forces every frame onto the sky
scale of the current one, which throws away the per-mission zoom the moment it
is applied. The script prints the line to send when they *should* be tied.

**Two ds9 behaviours that had to be worked around.**

- *Edit mode.* ds9 8.6 comes up in mode `none`, checked by launching a bare
  ds9 and asking it: the mouse is then inert, a tile cannot be clicked to become
  the current frame, and the wheel only ever zooms whichever frame was current
  already, which reads as the zoom having stopped working. The script sends
  `mode pointer`, which ds9 reports back as `region`.
- *XPA registration.* ds9 announces itself to the name server `xpans` on port
  14285, and everything here is driven through that. Measured on this machine:
  a ds9 that does register can take **over two minutes** to open the connection,
  on an unloaded machine, and one that fails to open it at all never recovers —
  a ds9 was seen alive for four minutes with no connection to 14285 at all. The
  script waits 180 s and then replaces that ds9 rather than giving up, twice. It
  also tests the name server by asking it a question rather than by `pgrep`,
  since a wedged server holds the port and answers nothing, and a `pgrep` guard
  would skip the very restart that fixes it.

**Checked.** Run end to end on `1FGL_J1018.6-5856` (3 Chandra, 3 XMM-Newton),
on `4U_1901+03` (2 ACIS and 1 HRC, which is how the HRC zoom was verified: the
three tiles carry the same grid) and on `GS_0834-430` (2 ACIS and 3 HRC, 21 s
with the blocked images against seven minutes hung without them). The saved
PNGs were read back and inspected. **Not checked:** any system with more than
17 frames; PSR B1259-63 at 41 frames has not been run.

### X-05 — observation lists

The current lists are `../Observations/Chandra/obsids_chandra.txt` (187) and
`../Observations/XMM-Newton/obsids_xmm.txt` (192). The original notebook wrote
lists in Output; those are not their present locations. Chandra obsids came
from accepted detection rows; XMM obsids from all detection rows sharing an
accepted SRCID. These reduction lists are distinct from Stage 1 target lists.

### R-02 — inventory and integrity limits

Checked now by directory/file names: 187 raw Chandra directories, 192 raw XMM
directories, 187 Chandra repro directories with 188 evt2 filenames, and 192
XMM repro directories with 586 imaging and 65 timing event filenames.
This establishes presence, not successful science processing of every file.
The earlier session recorded repair of nine truncated XMM archives and of
Chandra obsid 1556 (HTML error pages), and removal of three empty evt1 files
from interrupted reprocessing. These are historical repairs; no files were
deleted or downloaded since. Earlier raw-data checks covered names,
archive readability and magic bytes, not a byte-for-byte archive audit.
The old contradictory claim that Chandra reprocessing had not run is withdrawn.

### R-01 — current XMM reprocessing recipe

Actual script: `../Observations/XMM-Newton_repro/xmm-newton_repro.sh`, read on
26 September. It runs cifbuild, odfingest, emproc and epproc when event files
are absent. It creates per-exposure imaging light curves and PDF plots with
100-second bins, PATTERN==0, MOS PI>10000, pn 10000<PI<12000, using #XMMEA_EM
or #XMMEA_EP. Thresholds 0.35 (MOS) and 0.40 (pn) counts/s are reference lines
only. Axis labels/ticks use fontsize 14. Names include camera/exposure/mode;
plots are `rate_<tag>.pdf`. Existing event files/plots are skipped and logs
appended. Presence-based skipping is not a full completeness validation of a
partially processed observation.

### R-03 — FLAG/PATTERN, no solar-flare GTIs

Actual script: `../Observations/XMM-Newton_repro/xmm-newton_flag-pattern.sh`,
read on 26 September. It processes only *ImagingEvts.ds with an open filter;
Closed and CalClosed are skipped. MOS keeps #XMMEA_EM, 200<=PI<=12000,
PATTERN<=12; pn keeps #XMMEA_EP and FLAG==0, PATTERN==0 at 200–500 eV,
PATTERN<=4 above 500 eV through 12 keV. The lower edge was 300 eV until
26 September, when the user set it to 200 eV; the 561 filtered event lists
were deleted and rebuilt with the new limit, checked in the logs
(`PI>=200`). Outputs are
`flag-pattern_<full-original-event-stem>.fits`, with per-obsid logs.
There is **no tabgtigen call, no gti(...) selection and no solar-flare cut**.
The old bkgcorr/clean/GTI recipe is superseded and exists only in the archive.

Checked now by names: 561 flag-pattern files, 651 rate PDFs, zero clean_*
and zero gti_* FITS products. Existing timing PDFs can predate the imaging-only
rule. This is not a content audit of those 561 event lists. The prior header
audit recorded 25 closed-filter imaging files, and 36 open-filter MOS outer-CCD
imaging files with central CCD in timing mode. ONTIME/LIVETIME=0 in their
central-CCD keywords does not mean the outer CCDs have no exposure; per-CCD
EXPOSU/GTI extensions are the relevant evidence. Their suitability for an HMXB
is decided at extraction, not by discarding all such files.

Sources: `../Papers/watson2009.pdf` §4.3 for PATTERN; `../Papers/SASUSG.pdf`
Issue 19.0 §4.1.2 and §4.4.4 for modes and flare screening. The Closed/CalClosed
definitions were previously checked at the SOC watchout URL recorded in the
archive. Hard bright sources can contaminate whole-field >10 keV curves;
the source-exclusion question belongs to the deferred flare-screening step.

## 7. Repeating the catalogues' own fits (Stage 1)

### R-04 — inputs and target selection

`select_stage1_targets.py` rebuilds the CSC target list from
`../hmxb_5arcsec_csc.tsv` by (obsid, name), dropping the header line the file
repeats among its rows and keeping the rows with a published `powlaw_nh`.
This preserves the user's decision to keep the additional sources of the
updated cone-search pull; it is not restricted back to the 243-row
crossmatch. It gives 115 detections, 108 observations, 47 sources. No fresh
count cut and no positional match is made. Checked: the fitted rows are
exactly those above the catalogue's 150 net broad-band counts, minimum 157.7
and median 815.7, and none falls below it.

The 11 fitted 5XMM SRCIDs come from `Output/hmxb_5xmm_crossmatch.fits`; the
released catalogue is joined on SRCID and `SPECTRA` to give 60 detections
over 50 observations. That join is necessary because the `SPEC_` results sit
on the stacked-source rows while `SPECTRA` marks the detection rows that
carry a spectrum, and over the whole catalogue the two sets are disjoint.

`csc_products/`, `xmm_pps/` and `epic_rmf/` are read-only inputs, fetched by
`download_csc_products.py`, `download_xmm_pps.sh` and `download_epic_rmf.py`.
Checked: all 115 Chandra detections have all four of `pha3`, `arf3`, `rmf3`
and `reg3`; all 60 XMM detections have their EPIC spectra, 123 files; the 58
canned response matrices are readable and carry a MATRIX extension. The
XMM-Newton pipeline distributes no RMF, so each spectrum names a canned one
in `RESPFILE`. The pn matrices were taken as `v22.0`, the version the
repository offers, while the header names none: their agreement with the
catalogue pipeline calibration is **not** established.

### R-05 — Chandra fits

`refit_csc_stage1.py` reads each original `pha3` with its background, ARF and
RMF, groups to 16 counts, subtracts the background and fits
`xsphabs*powlaw1d` with angr/bcmc and `chi2datavar` over 0.5-7.0 keV,
following `../Papers/Evans2024.pdf` §3.13. A fit worse than 1.2 times the
published reduced statistic triggers a seeded `moncar`, the better solution
being kept; this catalogue-dependent trigger is an optimisation diagnostic,
not an independent success criterion. It fired on 18 of the 115. Our N_H is
converted once from the 10^22 cm^-2 of `xsphabs` to the 10^20 cm^-2 of the
catalogue's columns.

Result: 105 of 115 reach the published statistic. Among those, N_H falls
inside the published 68% limits for 92 and Gamma for 98, and the median
relative flux difference is +2.3%. Over all 115 the counts are 100 and 104 of
the 112 detections that carry limits. On obsid 88, `2CXO J112115.1-603725`,
our reduced chi-square is 1.436 against their 1.438. Reproducing a statistic
is evidence that the grouping, subtraction and variance match; it is not by
itself proof that the whole configuration was reproduced.

Checked: the same script run serially and 12 at a time gives identical
reduced chi-squares for all 115 and parameters agreeing to better than 1e-6
for 114; the one exception differs by 1.3e-4 and goes through the stochastic
`moncar`. Wall time falls from about 25 minutes to 38 seconds.

**Not checked:** the ten detections that never reach the published statistic;
and no confidence limits are computed on our own parameters, so the
comparison is our best fit against their interval, not interval against
interval.

### R-06 — XMM-Newton fits

`refit_xmm_stage1.py` follows the **local draft §6**: positive source,
background and area-scaled net counts in 0.2-12 keV; grouping to at least one
count; `constant*cflux*phabs*powerlaw` with the Cash statistic; priors
log-uniform in N_H over 1e19-1e25 cm^-2, uniform in Gamma over [1,3] and
uniform in log10(observed flux) over [-15,-9]; the pn constant fixed at unity
and the MOS constant free in a joint fit. It runs in `env_bxa`, a venv on
HEASoft's python with `--system-site-packages` so that it sees `xspec`
without anything being installed into HEASoft: XSPEC 12.15.1, UltraNest 4.5.2,
BXA 5.1.1.

**Which camera the catalogue used is recoverable, and is no longer assumed.**
The catalogue does not name the instrument, but `SPEC_DOF_PL` equals the
number of fitted bins minus the free parameters, three for one instrument and
four for two. Counting our own bins in each combination reproduces the
published degrees of freedom **exactly for all eleven sources**: pn only for
3009282010100002, 3012270010100009, 3060482030100046 and 3072796120100002;
MOS only for 3012270010100011, 3014687010100002 and 3050528010100001; pn and
MOS for 3020120020100001, 3050332030100003, 3050619010100001 and
3060385010100002 — which are precisely the four that carry a published
`SPEC_IIN_PL`. The script now selects on this match instead of on the
highest signal-to-noise rule, which disagreed with it for four sources.

Two things follow from the same match, and both were previously recorded as
assumptions. The **fit band is 0.3-10 keV**, since our bin counts in that band
are the ones that reproduce the published degrees of freedom. And our
**grouping reproduces theirs bin for bin**, since the counts agree exactly
rather than approximately.

**The flux is quoted over the fitted band, not the detection band.** Taking
our fitted model and rescaling the flux from 0.2-12 keV to 0.3-10 keV
reproduces `SPEC_FLUX_PL` to better than 1% on every source tested:
3060385010100002 0.992, 3050332030100003 0.997, 3020120020100001 1.007,
3050619010100001 0.992. This removed a systematic 16% excess that had been
left unexplained. `cflux` now integrates over the fitted band.

Remaining explicit assumptions: the abundance table and cross-sections, which
the draft does not state. An attached Poisson background with XSPEC `cstat`
uses W-stat and is not Gaussian background subtraction; the draft's wording
about Cash and background subtraction is insufficient to establish exact
equivalence. Published fitted sources are retained without re-applying a
preliminary p-value screen, and the catalogue's permutation p-value is not
recomputed.

The online 5XMM User Guide §3.6, read on 26 September, materially conflicts
with the local draft: it describes merged spectra, tightened preliminary-fit
priors, an empirical background model and 5/95 percentiles, whereas the local
draft describes highest-S/N single spectra, fixed priors and 16/84
percentiles. This work follows the requested local paper and records the
conflict; it does not mix the two recipes or claim the discrepancy resolved.
Source: https://xmmssc.irap.omp.eu/Catalogue/5XMM-DR15/5XMM-DR15_Catalogue_User_Guide.html
XSPEC statistical convention:
https://heasarc.gsfc.nasa.gov/docs/software/xspec/manual/XSappendixStatistics.html

### R-07 — two errors in the comparison, corrected

Both were in the comparison rather than in the fit, and both are corrected in
`refit_xmm_stage1.py`.

BXA's log-uniform transform stores **log10(nH)** in sampling space; its
aftertransform is `10**x`, which is why the posterior column is named
`log(nH)`. Reading that column as linear and scaling it to cm^-2 made N_H
come out about a factor of 13 low, and negative for five sources, which a
log-uniform prior over 1e19-1e25 cm^-2 cannot produce. Source code:
`env_bxa/lib/python3.13/site-packages/bxa/xspec/priors.py`;
https://johannesbuchner.github.io/BXA/_modules/bxa/xspec/priors.html

5XMM `SPEC_*_ERR_LO/UP_PL` are **percentile endpoints**, not offsets to
subtract and add: for 3060385010100002 they read 1.02e22 and 1.56e22 around
a median of 1.27e22, and the Gamma pair brackets its median in the same way.
They had been subtracted, which would have made every published interval
wrong.

A hypothesis was tested and rejected along the way. The draft states 0.3-10
keV only for its Levenberg-Marquardt pre-screen and never states the band of
the Bayesian fit, so the band was the first suspect for the N_H difference.
Refitting 3060385010100002 over 0.2-12 keV moved N_H by 1% and left Gamma,
the flux and the MOS normalisation unchanged. The band was not the cause; the
decoding above was.

### R-08 — current status

The Chandra side of Stage 1 is complete and its numbers are in R-05 and in
`Output/stage1_csc_refit.tsv`. The XMM-Newton side is being refitted from
scratch with the camera selection of R-06 and the corrected flux band; the
earlier XMM numbers were produced before both corrections and are not
carried forward. No XMM success count is asserted until that run finishes
and its per-source results are checked. Stage 2, extraction from our own
reprocessed event files, has not been started for either mission.

## 8. The manuscript

### M-01 — manuscript synchronization

The manuscript's reduction section and build describe the earlier drafting
state and are now out of date in one specific way: the draft box on page 1
states that the *Chandra* reprocessing and the flare, FLAG and PATTERN
screening "have not been run yet", and both have since run — 188 Chandra
reprocessing directories with 188 `evt2` files, and 561 `flag-pattern_*`
files with no `clean_*` or `gti_*`, checked by name on 26 September. That
sentence must be corrected when the manuscript is next revised, and must not
be used to override the inventory above. Stage 1 results and remaining
configuration limitations below are the current record; manuscript revision
and Stage 2 extraction remain separate work.
