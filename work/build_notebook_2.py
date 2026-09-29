"""Emit 2_XrayCrossmatch.ipynb.

The X-ray half of the sample definition: the Chandra (CSC 2.1.1) and
XMM-Newton (5XMM-DR15) counterparts of the 109 HMXBs that have a Gaia
DR3 counterpart and a Bailer-Jones distance.

The recipe is the one of ~/Doctorado/2025/XRays_ML/M101/1_Crossmatch.ipynb,
cells 5 and 29: rescale each catalogue's error region to a covariance
matrix with `matrix_rotation`, sum the two covariances, and accept on the
Mahalanobis distance.  Two differences with that notebook are deliberate
and are stated in the notebook itself: the Omega rotation is applied
inside `crossmatch.Mahalanobis_distance_squared` here, and the position
of the HMXB is the Gaia one rather than the F23 one.
"""

import json

cells = []


def _lines(src):
    return src.strip("\n").splitlines(keepends=True)


def md(src):
    cells.append({"cell_type": "markdown", "metadata": {},
                  "source": _lines(src)})


def code(src):
    cells.append({
        "cell_type": "code", "metadata": {}, "execution_count": None,
        "outputs": [], "source": _lines(src),
    })


# ---------------------------------------------------------------- header
md(r"""
# *Chandra* and *XMM-Newton* counterparts of the HMXB sample

Input: the 109 systems of `Output/hmxb_sample_positions.vot`, i.e. the
F23 sources with an accepted Gaia DR3 counterpart *and* a Bailer-Jones
distance. Output: for each of the two X-ray catalogues, the accepted
counterparts with **every column the catalogue publishes**, plus the
Gaia identifier, the `Notes` column and the statistic of the Gaia
cross-match.

The statistic is the one of Pineau et al. (2011), Appendix A, as
implemented in `crossmatch.py`, and the recipe is the one of
`~/Doctorado/2025/XRays_ML/M101/1_Crossmatch.ipynb` (cells 5 and 29).
Two departures from that notebook are deliberate:

1. **The Omega rotation is not added to the position angle.** In the
   M101 notebook `Omega_function` is called explicitly and its result is
   added to the ellipse position angle
   (`beta_chandra = err_ellipse_ang + Omega`), because the
   `Mahalanobis_distance_squared` of that notebook does not rotate
   anything itself. The `crossmatch.py` used here computes both Omega
   angles internally and rotates *both* covariance matrices. Adding
   Omega to `beta` as well would apply it twice, so `beta` is the
   catalogue position angle alone.
2. **The HMXB position is the Gaia one.** The cone search that produced
   the *Chandra* table was made around the Gaia positions, and the Gaia
   position of the adopted counterpart is the best position we have for
   the system (sub-mas, against the arcsecond-scale `ePos` of F23). The
   match with the F23 position and its 90% error, which is what the M101
   notebook does, is computed as well and kept in the column
   `dM2_xray_F23`, so the alternative acceptance can be read off
   directly.
""")

code(r"""
import gzip
import os
import shutil
import tempfile
import time

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.table import Table

import crossmatch

pd.set_option('display.width', 200)
pd.set_option('display.max_columns', 40)

# Acceptance threshold: 99% point of chi-squared with two degrees of
# freedom, the same one used for the Gaia cross-match.
ACCEPT = 9.21

# Rescaling factors, one per error region.  See the section below.
CV_F23 = 4.61      # F23 ePos: 90% confidence, 2 dof
CV_CSC = 5.99      # CSC err_ellipse_r0/r1: 95% confidence, 2 dof
CV_XMM_ELL = 2.30  # 5XMM ELL_MAJOR/MINOR: Delta C = 2.3, 68% for 2 dof
SIG_SYS_XMM = 0.88 # arcsec, 5XMM systematic astrometric error, per axis

XMM_GZ = "/home/marina/Doctorado/2026/HMXB_project/5XMM-DR15/5XMM_DR15.fits.gz"
XMM_FITS = os.environ.get(
    "XMM_FITS_CACHE", os.path.join(tempfile.gettempdir(), "5XMM_DR15.fits"))
KEEP_XMM_FITS = False  # the decompressed catalogue is 8.5 GB
""")

# ------------------------------------------------- what the catalogues say
md(r"""
## 1. What each catalogue publishes as a position error

This had to be established before any rescaling, by reading the
catalogue papers and the column metadata. `matrix_rotation(a, b, cv, pa)`
divides both semi-axes by `sqrt(cv)`, so `cv` must be the chi-squared
value the published region corresponds to, and `cv = 1` means the
published number is already a standard deviation.

**F23 (`ePos`).** Circular, 90% confidence, two degrees of freedom:
`cv = 4.61`. This is the value used for the Gaia cross-match and it is
used again here for the secondary statistic.

**CSC 2.1.1 (`err_ellipse_r0`, `err_ellipse_r1`, `err_ellipse_ang`).**
The `#Column` header of the cone-search file reads *"Major radius of the
95% confidence level position error ellipse"*, *"Minor radius of the 95%
confidence level position error ellipse"* and *"Position angle (ref.
local true north) of the major axis of the 95% confidence level error
ellipse"*. Two degrees of freedom at 95% is `cv = 5.99`, the value used
in the M101 notebook.

**5XMM-DR15.** This is the case that had to be checked, and the answer
is that the two position-error columns need *different* treatment.

* `RADEC_ERR` is documented in the FITS header as *"Quadratic sum of
  1-sigma errors in RA and DEC"*, and Sect. 4.1.2 of the catalogue paper
  writes the astrometric standard deviation per axis as
  $\sigma = \sqrt{\sigma_\mathrm{stat}^2 + \sigma_\mathrm{sys}^2}$ with
  $\sigma_\mathrm{stat} = \texttt{RADEC\_ERR}/\sqrt{2}$. So it is **not**
  a confidence region and must **not** be rescaled by a chi-squared
  critical value; it must be divided by $\sqrt 2$ to become a per-axis
  standard deviation.
* The same section gives the systematic term:
  $\sigma_\mathrm{sys} = 0\farcs88 \pm 0\farcs01$, inferred against Gaia
  and unWISE QSOs, *"larger than the one derived for 4XMM-DR9s
  (0.227")... because the 5XMM data were not rectified astrometrically
  when producing DR15"*. It is added in quadrature, per axis.
* `ELL_MAJOR`, `ELL_MINOR`, `ELL_PA` *are* a confidence region. Sect. 3
  of the paper: *"The parameters of an ideal ellipse containing 68% of
  the positions could be calculated from any three points on this
  ellipse. We thus determine $C_\mathrm{min}+2.3$ (1$\sigma$ for two
  parameters) along the x-axis, the y-axis, and their diagonal, and
  derive the semi-axes and orientation of an ideal error ellipse from
  them."* $\Delta C = 2.3$ with two degrees of freedom is the 68.3%
  contour, so `cv = 2.30`.

**Not established by reading.** The reference direction of `ELL_PA` is
not stated either in the FITS header (*"Rotation angle of the positional
error ellipse"*) or anywhere in the catalogue paper, which never names
the three `ELL_*` columns. The elliptical statistic below therefore
assumes the same convention as the CSC, position angle from local true
north, and it is reported as a secondary quantity only: the accepted set
is decided on the circular `RADEC_ERR` route, which needs no angle and
is the one the paper writes down.
""")

# --------------------------------------------------------------- the sample
md("## 2. The sample")

code(r"""
_S = {'ID': str, 'GaiaDR3_adopted': str, 'Notes': str}
sample = pd.read_csv("Output/hmxb_gaia_adopted.csv", dtype=_S).fillna('')
gaia = pd.read_csv("Output/gaia_data_Tables/hmxbs_gaia_data.csv",
                   dtype={'DR3Name': str, 'ID': str})
f23 = pd.read_csv("Output/f23_flagged.csv", dtype={'ID': str})

# The 109 with a distance are the sample; the other four accepted Gaia
# counterparts have no Bailer-Jones entry and are not carried here.
sample = sample.rename(columns={'dM2': 'dM2_gaia_match'})
sample = sample.merge(gaia[['ID', 'DR3Name', 'RA_ICRS', 'DE_ICRS',
                            'e_RA_ICRS', 'e_DE_ICRS', 'RA_DEC_CORR']],
                      on='ID', how='inner')
sample = sample.merge(f23[['ID', 'RAdeg', 'DEdeg', 'ePos']], on='ID', how='left')
assert sample['ID'].is_unique
assert (sample['DR3Name'] == sample['GaiaDR3_adopted']).all()
print("systems in the sample:", len(sample))


def gaia_cov(row):
    '''Gaia covariance in arcsec^2.  The catalogue errors are 1-sigma,
    so there is nothing to rescale (cv = 1), and RA_DEC_CORR is the
    correlation between the two coordinates.'''
    sa = row['e_RA_ICRS'] / 1000.0
    sd = row['e_DE_ICRS'] / 1000.0
    c = row['RA_DEC_CORR'] * sa * sd
    return sa**2, c, c, sd**2


def f23_cov(row):
    '''F23 ePos is a circular 90% region, two degrees of freedom.'''
    return crossmatch.matrix_rotation(row['ePos'], row['ePos'], CV_F23, 0.0)
""")

# ------------------------------------------------------------------ Chandra
md(r"""
## 3. *Chandra*: CSC 2.1.1

`../hmxb_5arcsec_csc.tsv` is a 5" cone search of CSC 2.1.1 around
each of the 109 positions of `Output/hmxb_sample_positions.vot`, run
through CSCview. It is kept exactly as delivered, and it carries
everything the match needs: the `master_source` position, the
`master_source` error ellipse and the observation-level rows, one per
`obsid`, so each master source appears several times.

The position comes as sexagesimal strings, so it is converted with
`sexagesimal_to_decimal`, the function of cell 5 of
`~/Doctorado/2025/XRays_ML/M101/1_Crossmatch.ipynb`.
""")

code(r"""
csc = pd.read_csv("../hmxb_5arcsec_csc.tsv", sep='\t', comment='#',
                  dtype=str)
csc['usrid'] = csc['usrid'].str.strip()
csc['name'] = csc['name'].str.strip()
print("rows:", len(csc), "| columns:", csc.shape[1],
      "| HMXBs with a CSC source inside 2'':", csc['usrid'].nunique(),
      "| distinct 2CXO sources:", csc['name'].nunique())
assert set(csc['usrid']) <= set(sample['ID'])
""")

code(r"""
def sexagesimal_to_decimal(mu_alpha_sexagesimal, mu_delta_sexagesimal):
    # VERBATIM from cell 5 of M101/1_Crossmatch.ipynb.
    coords = SkyCoord(ra=mu_alpha_sexagesimal, dec=mu_delta_sexagesimal,
                      unit=('hourangle', 'deg'))
    mu_alpha_decimal = coords.ra.deg
    mu_delta_decimal = coords.dec.deg
    return mu_alpha_decimal, mu_delta_decimal


csc['ra_decimal'], csc['dec_decimal'] = sexagesimal_to_decimal(
    csc['ra'].values, csc['dec'].values)
for c in ('err_ellipse_r0', 'err_ellipse_r1', 'err_ellipse_ang'):
    csc[c] = csc[c].astype(float)

# One master source per name: position and ellipse must not depend on
# which observation row they were read from.
ell = csc[['name', 'ra_decimal', 'dec_decimal', 'err_ellipse_r0',
           'err_ellipse_r1', 'err_ellipse_ang']].drop_duplicates()
assert len(ell) == csc['name'].nunique()
assert ell.notna().all().all()
print("master sources:", len(ell),
      "| ellipse semi-major (arcsec): median %.2f  max %.2f"
      % (ell['err_ellipse_r0'].median(), ell['err_ellipse_r0'].max()))
print(ell.head(3).to_string(index=False))
""")

md(r"""
The match itself. `beta` is `err_ellipse_ang` alone, for the reason given
at the top of the notebook.
""")

code(r"""
pairs = csc[['usrid', 'name']].drop_duplicates().merge(ell, on='name', how='left')
assert pairs[['err_ellipse_r0', 'err_ellipse_r1', 'err_ellipse_ang']].notna().all().all()

g_by_id = sample.set_index('ID')
rows = []
for _, p in pairs.iterrows():
    s = g_by_id.loc[p['usrid']]
    X = crossmatch.matrix_rotation(p['err_ellipse_r0'], p['err_ellipse_r1'],
                                   CV_CSC, p['err_ellipse_ang'])
    dM2 = crossmatch.Mahalanobis_distance_squared(
        *gaia_cov(s), *X, s['RA_ICRS'], s['DE_ICRS'],
        p['ra_decimal'], p['dec_decimal'])
    dM2_f23 = crossmatch.Mahalanobis_distance_squared(
        *f23_cov(s), *X, s['RAdeg'], s['DEdeg'],
        p['ra_decimal'], p['dec_decimal'])
    rows.append({'usrid': p['usrid'], 'name': p['name'],
                 'sep_gaia': crossmatch.haversine_distance(
                     s['RA_ICRS'], s['DE_ICRS'],
                     p['ra_decimal'], p['dec_decimal']),
                 'dM2_xray': dM2, 'dM2_xray_F23': dM2_f23})

csc_dm2 = pd.DataFrame(rows)
csc_dm2['accepted'] = csc_dm2['dM2_xray'] < ACCEPT
print("candidate (HMXB, 2CXO) pairs      :", len(csc_dm2))
print("accepted at dM2 < 9.21            :", int(csc_dm2['accepted'].sum()))
print("HMXBs with >= 1 accepted source   :",
      csc_dm2.loc[csc_dm2['accepted'], 'usrid'].nunique(), "of", len(sample))
print("same acceptance with F23 as origin:",
      bool(((csc_dm2['dM2_xray_F23'] < ACCEPT) == csc_dm2['accepted']).all()))
csc_dm2.sort_values('dM2_xray').head(10)
""")

md(r"""
### The *Chandra* table

Every column of the cone-search file is kept, together with the decimal
positions derived from it, the Gaia identifier, `Notes`, the statistic
of the Gaia cross-match and the two statistics of this one. One row per
(observation, source), as delivered.
""")

code(r"""
keep = csc_dm2.loc[csc_dm2['accepted'], ['usrid', 'name', 'sep_gaia',
                                         'dM2_xray', 'dM2_xray_F23']]
csc_out = (csc.merge(keep, on=['usrid', 'name'], how='inner')
              .merge(sample[['ID', 'GaiaDR3_adopted', 'Notes', 'dM2_gaia_match']],
                     left_on='usrid', right_on='ID', how='left')
              .drop(columns='ID'))
assert csc_out['GaiaDR3_adopted'].notna().all()
csc_out.to_csv("Output/hmxb_csc_crossmatch.csv", index=False)
print("rows:", len(csc_out), "| columns:", csc_out.shape[1],
      "| HMXBs:", csc_out['usrid'].nunique(),
      "| 2CXO sources:", csc_out['name'].nunique(),
      "| HMXBs with more than one:",
      int((csc_out.groupby('usrid')['name'].nunique() > 1).sum()))
csc_out[['usrid', 'name', 'obsid', 'GaiaDR3_adopted', 'dM2_gaia_match',
         'dM2_xray', 'Notes']].head(6)
""")

# ---------------------------------------------------------------------- XMM
md(r"""
## 4. *XMM-Newton*: 5XMM-DR15

No cone search was run against a server here: the whole catalogue is on
disk (`5XMM_DR15.fits.gz`, 3 397 248 rows, 421 columns). Of those rows,
818 656 are the stacked *sources* (`OBS_ID` empty, one per `SRCID`) and
the rest are the individual detections that contribute to them. The
match is made against the source rows, which is the level that
corresponds to CSC `master_source`.

The efficient part is the pre-filter. Sorting the catalogue by right
ascension once and taking, for each HMXB, the slice
`searchsorted(ra, ra0 -+ R/cos(dec))` followed by a declination cut
reduces 818 656 rows to a few hundred candidates in about a second; the
Mahalanobis distance is then computed on those only. The radius `R` is
not chosen by hand: with the Gaia error negligible,
$d_M^2 \geq (\Delta/\sigma_\mathrm{max})^2$, so no pair with
$d_M^2 < 9.21$ can lie farther than
$\sqrt{9.21}\,\sigma_\mathrm{max}$, and `R` is set to that. The
pre-filter therefore cannot discard an acceptable pair.
""")

code(r"""
if not os.path.exists(XMM_FITS):
    t0 = time.time()
    with gzip.open(XMM_GZ, 'rb') as src, open(XMM_FITS, 'wb') as dst:
        shutil.copyfileobj(src, dst, length=32 * 1024 * 1024)
    print(f"decompressed in {time.time() - t0:.0f} s")
hdul = fits.open(XMM_FITS, memmap=True)
xmm = hdul[1].data
print("rows:", len(xmm), "| columns:", len(hdul[1].columns))

is_source = np.char.strip(np.asarray(xmm['OBS_ID'])) == ''
src_idx = np.nonzero(is_source)[0]
x_ra = np.asarray(xmm['RA'], float)[src_idx]
x_de = np.asarray(xmm['DEC'], float)[src_idx]
x_re = np.asarray(xmm['RADEC_ERR'], float)[src_idx]
print("stacked sources:", len(src_idx), "| detections:", len(xmm) - len(src_idx))

# Per-axis standard deviation, catalogue paper Eq. (1).
x_sig = np.sqrt((x_re / np.sqrt(2.0))**2 + SIG_SYS_XMM**2)
print("sigma per axis (arcsec): median %.2f  99th %.2f  max %.2f  undefined %d"
      % (np.nanmedian(x_sig), np.nanpercentile(x_sig, 99), np.nanmax(x_sig),
         int(np.isnan(x_sig).sum())))
RADIUS = np.sqrt(ACCEPT) * np.nanmax(x_sig)
print("pre-filter radius (arcsec): %.1f" % RADIUS)
""")

code(r"""
order = np.argsort(x_ra)
ra_sorted = x_ra[order]


def ra_windows(ra0, dec0, radius_arcsec):
    '''Indices into the RA-sorted array whose RA is within the window,
    wrapping at 0/360.'''
    dra = radius_arcsec / 3600.0 / max(np.cos(np.deg2rad(dec0)), 1e-6)
    lo, hi = ra0 - dra, ra0 + dra
    spans = [(lo, hi)]
    if lo < 0:
        spans = [(0.0, hi), (lo + 360.0, 360.0)]
    elif hi > 360.0:
        spans = [(lo, 360.0), (0.0, hi - 360.0)]
    out = []
    for a, b in spans:
        i, j = np.searchsorted(ra_sorted, [a, b])
        out.append(order[i:j])
    return np.concatenate(out)


t0 = time.time()
cand = []
for _, s in sample.iterrows():
    w = ra_windows(s['RA_ICRS'], s['DE_ICRS'], RADIUS)
    w = w[np.abs(x_de[w] - s['DE_ICRS']) < RADIUS / 3600.0]
    for k in w:
        dh = crossmatch.haversine_distance(s['RA_ICRS'], s['DE_ICRS'],
                                           x_ra[k], x_de[k])
        if dh <= RADIUS:
            cand.append((s['ID'], int(src_idx[k]), dh))
print(f"candidates within {RADIUS:.0f}'': {len(cand)}  ({time.time() - t0:.1f} s)")
""")

code(r"""
rows = []
for hid, row, dh in cand:
    s = g_by_id.loc[hid]
    ra_x = float(xmm['RA'][row])
    de_x = float(xmm['DEC'][row])
    sig = float(np.sqrt((xmm['RADEC_ERR'][row] / np.sqrt(2.0))**2 + SIG_SYS_XMM**2))

    # Circular route: RADEC_ERR is already 1-sigma, so cv = 1.
    Xc = crossmatch.matrix_rotation(sig, sig, 1.0, 0.0)
    dM2 = crossmatch.Mahalanobis_distance_squared(
        *gaia_cov(s), *Xc, s['RA_ICRS'], s['DE_ICRS'], ra_x, de_x)
    dM2_f23 = crossmatch.Mahalanobis_distance_squared(
        *f23_cov(s), *Xc, s['RAdeg'], s['DEdeg'], ra_x, de_x)

    # Elliptical route: 68% for two parameters, plus the systematic term,
    # which is isotropic and so commutes with the rotation.
    a = float(xmm['ELL_MAJOR'][row])
    b = float(xmm['ELL_MINOR'][row])
    pa = float(xmm['ELL_PA'][row])
    if np.isfinite([a, b, pa]).all() and a > 0:
        e11, e12, e21, e22 = crossmatch.matrix_rotation(a, b, CV_XMM_ELL, pa)
        dM2_ell = crossmatch.Mahalanobis_distance_squared(
            *gaia_cov(s), e11 + SIG_SYS_XMM**2, e12, e21, e22 + SIG_SYS_XMM**2,
            s['RA_ICRS'], s['DE_ICRS'], ra_x, de_x)
    else:
        dM2_ell = np.nan

    rows.append({'ID': hid, 'row': row,
                 'IAUNAME': str(xmm['IAUNAME'][row]).strip(),
                 'SRCID': int(xmm['SRCID'][row]),
                 'xmm_ra': ra_x, 'xmm_de': de_x, 'sep_gaia': dh,
                 'sigma_axis': sig, 'dM2_xray': dM2,
                 'dM2_xray_F23': dM2_f23, 'dM2_xray_ellipse': dM2_ell})

xmm_dm2 = pd.DataFrame(rows)
xmm_dm2['accepted'] = xmm_dm2['dM2_xray'] < ACCEPT
print("candidates                        :", len(xmm_dm2))
print("  with RADEC_ERR undefined        :", int(xmm_dm2['sigma_axis'].isna().sum()))
print("  with the ellipse undefined      :", int(xmm_dm2['dM2_xray_ellipse'].isna().sum()))
print("accepted at dM2 < 9.21 (circular) :", int(xmm_dm2['accepted'].sum()))
print("HMXBs with >= 1 accepted source   :",
      xmm_dm2.loc[xmm_dm2['accepted'], 'ID'].nunique(), "of", len(sample))
print("accepted on the ellipse instead   :",
      int((xmm_dm2['dM2_xray_ellipse'] < ACCEPT).sum()),
      "| decisions that differ:",
      int(((xmm_dm2['dM2_xray_ellipse'] < ACCEPT) != xmm_dm2['accepted']).sum()))
print("same acceptance with F23 as origin:",
      bool(((xmm_dm2['dM2_xray_F23'] < ACCEPT) == xmm_dm2['accepted']).all()))
""")

md(r"""
### The *XMM-Newton* table

All 421 catalogue columns of the accepted source rows, plus the same
five added columns as the *Chandra* table. Written as FITS, which keeps
the catalogue's own types, and as CSV.
""")

code(r"""
acc = xmm_dm2[xmm_dm2['accepted']].reset_index(drop=True)
sub = Table(xmm[np.asarray(acc['row'])])
xmm_out = sub.to_pandas()
for col in ('sep_gaia', 'sigma_axis', 'dM2_xray', 'dM2_xray_F23',
            'dM2_xray_ellipse'):
    xmm_out[col] = acc[col].values
xmm_out['ID'] = acc['ID'].values
xmm_out = xmm_out.merge(
    sample[['ID', 'GaiaDR3_adopted', 'Notes', 'dM2_gaia_match']],
    on='ID', how='left')
assert xmm_out['GaiaDR3_adopted'].notna().all()

xmm_out.to_csv("Output/hmxb_5xmm_crossmatch.csv", index=False)
Table.from_pandas(xmm_out).write("Output/hmxb_5xmm_crossmatch.fits",
                                 overwrite=True)
print("rows:", len(xmm_out), "| columns:", xmm_out.shape[1],
      "| HMXBs:", xmm_out['ID'].nunique(),
      "| HMXBs with more than one:",
      int((xmm_out.groupby('ID')['IAUNAME'].nunique() > 1).sum()))
xmm_out[['ID', 'IAUNAME', 'SRCID', 'GaiaDR3_adopted', 'dM2_gaia_match',
         'dM2_xray', 'Notes']].head(6)
""")

md(r"""
### The observations to download

One observation identifier per line, for each telescope: the `obsid`
column of the accepted *Chandra* rows, and, for *XMM-Newton*, the
`OBS_ID` of every detection that contributes to an accepted stacked
source, which is where the identifiers live since the source rows
themselves carry none.
""")

code(r"""
cxo_obsids = sorted({o.strip() for o in csc_out['obsid'] if o.strip()})
with open("../Observations/Chandra/obsids_chandra.txt", "w") as fh:
    fh.write("\n".join(cxo_obsids) + "\n")

acc_srcid = set(acc['SRCID'])
srcid_all = np.asarray(xmm['SRCID'])
obs_all = np.char.strip(np.asarray(xmm['OBS_ID']))
in_acc = np.isin(srcid_all, list(acc_srcid)) & (obs_all != '')
xmm_obsids = sorted(set(obs_all[in_acc].tolist()))
with open("../Observations/XMM-Newton/obsids_xmm.txt", "w") as fh:
    fh.write("\n".join(xmm_obsids) + "\n")

print("Chandra observations:", len(cxo_obsids))
print("XMM-Newton observations:", len(xmm_obsids))
assert all(o and o.strip() == o for o in cxo_obsids + xmm_obsids)
""")

code(r"""
hdul.close()
if not KEEP_XMM_FITS and os.path.exists(XMM_FITS):
    os.remove(XMM_FITS)
    print("removed the decompressed catalogue:", XMM_FITS)
""")

# --------------------------------------------------------------- comparison
md(r"""
## 5. The identifiers, against F23 and SIMBAD

Same comparison as for the Gaia identifier. The F23 columns are read as
strings out of the fixed-width record, at the byte ranges the CDS ReadMe
declares: `1271-1291 A21 XMM` (*"Identifier in 4XMM DR11"*) and
`1341-1362 A22 Chandra` (*"Identifier in Chandra CSC2"*). The SIMBAD
identifiers are the ones collected in notebook 1.

Two things make a plain string comparison insufficient on the
*XMM-Newton* side. F23 quotes **4XMM-DR11** and we match against
**5XMM-DR15**: the prefix differs by construction, and the
`Jhhmmss.s+ddmmss` root is built from the position, which moved between
the two releases. The same holds, more weakly, between CSC 2 and CSC
2.1.1. The comparison is by name. It is made on the root, so that the
prefix alone never counts as a difference, and a root that differs only
because the source was redesignated in the newer release does not count
as a difference either: a designation encodes a position, quantised at
0.1 s in right ascension and 1" in declination, so two designations of
one source can differ by about 2" from the quantisation alone, and any
two roots whose encoded positions agree within 3" are taken to name the
same source. Nothing but the names enters this. The identifier kept is
the 5XMM-DR15 one, or the CSC 2.1.1 one.

An entry from another mission is not an identifier of the catalogue
being compared, and counts as no identifier rather than as a
difference.
""")

code(r"""
records = [l for l in open("Input/HMXB/tablea.dat").read().split('\n') if l.strip()]


def _f(lo, hi):
    return [l[lo - 1:hi].strip() for l in records]


def _num(lo, hi):
    return pd.to_numeric(pd.Series(_f(lo, hi)), errors='coerce')


f23_x = pd.DataFrame({
    'ID': _f(1, 23),
    'F23_XMM': _f(1271, 1291), 'F23_XMM_RA': _num(1293, 1303),
    'F23_XMM_DE': _num(1305, 1314),
    'F23_CXO': _f(1341, 1362), 'F23_CXO_RA': _num(1364, 1383),
    'F23_CXO_DE': _num(1385, 1404),
})
# The positions are read only to report, further down, how far apart the
# two catalogues place a source whose name differs.  They take no part
# in the comparison, which is by name.

_D = {'ID': str, 'Chandra_simbad': str, 'XMM_simbad': str}
simbad = pd.read_csv("Output/simbad_ids.csv", dtype=_D).fillna('')
print("F23 rows with an XMM identifier    :", int((f23_x['F23_XMM'] != '').sum()))
print("F23 rows with a Chandra identifier :", int((f23_x['F23_CXO'] != '').sum()))
print("SIMBAD, in the sample: Chandra", int((simbad['Chandra_simbad'] != '').sum()),
      "| XMM", int((simbad['XMM_simbad'] != '').sum()))
""")

code(r"""
import re

ROOT = re.compile(r'J\d{6}\.\d[+-]\d{6}')


def root(name):
    '''The Jhhmmss.s+ddmmss part of an X-ray source name, without the
    catalogue prefix and without the CSC ambiguity suffix.'''
    m = ROOT.search(str(name))
    return m.group(0) if m else ''


def roots(field):
    '''SIMBAD returns a ';'-joined list; keep every root it contains.'''
    return {root(x) for x in str(field).split(';') if root(x)}


CXO_FAMILY = re.compile(r'^(2CXO|CXO|CXOU|CXOGSG)\b')
XMM_FAMILY = re.compile(r'^([1-5]XMM|XMMU|XMMSL)\b')

# Two designations of the same source, issued by different releases of
# the same catalogue, differ because the name is recomputed from an
# updated position.  A name is quantised at 0.1 s in right ascension,
# which is 1.5'', and at 1'' in declination, so two designations of one
# source can differ by about 2'' from the quantisation alone.
RENAME_TOL = 3.0  # arcsec


def encoded_position(r):
    '''The position a Jhhmmss.s+ddmmss root encodes, in degrees.'''
    ra = (int(r[1:3]) + int(r[3:5]) / 60 + float(r[5:9]) / 3600) * 15.0
    sign = 1.0 if r[9] == '+' else -1.0
    de = sign * (int(r[10:12]) + int(r[12:14]) / 60 + int(r[14:16]) / 3600)
    return ra, de


def same_designation(a, b):
    '''True when two roots name the same source: identical, or the same
    source redesignated between releases.'''
    if a == b:
        return True
    ra1, de1 = encoded_position(a)
    ra2, de2 = encoded_position(b)
    return crossmatch.haversine_distance(ra1, de1, ra2, de2) < RENAME_TOL


def compare_set(ours_roots, theirs_roots):
    '''Compare by name.  Every source of a system that passes the
    cross-match is kept, so it is enough that one of them carries the
    identifier.  Returns 'yes' when a name is identical, 'renamed' when
    it is the same source under a redesignation -- which is not a
    difference -- 'NO' otherwise, '-' when there is nothing to compare.'''
    if not theirs_roots:
        return '-'
    if ours_roots & theirs_roots:
        return 'yes'
    if any(same_designation(a, b) for a in ours_roots for b in theirs_roots):
        return 'renamed'
    return 'NO'


def agreed(v):
    '''A redesignation is a coincidence, not a difference.'''
    return 'yes' if v == 'renamed' else v
""")

md("### *Chandra*")

code(r"""
best = (csc_dm2[csc_dm2['accepted']]
        .sort_values('dM2_xray').groupby('usrid').first().reset_index()
        .rename(columns={'usrid': 'ID', 'name': 'XM_CXO'}))
cxo = (best[['ID', 'XM_CXO', 'dM2_xray']]
       .merge(f23_x[['ID', 'F23_CXO', 'F23_CXO_RA', 'F23_CXO_DE']], on='ID', how='left')
       .merge(simbad[['ID', 'Chandra_simbad']], on='ID', how='left')
       .merge(sample[['ID', 'RA_ICRS', 'DE_ICRS', 'Notes']], on='ID', how='left'))
for c in ('F23_CXO', 'Chandra_simbad', 'Notes'):
    cxo[c] = cxo[c].fillna('')
cxo = cxo.merge(ell[['name', 'ra_decimal', 'dec_decimal']],
                left_on='XM_CXO', right_on='name', how='left').drop(columns='name')

acc_roots = (csc_dm2[csc_dm2['accepted']].groupby('usrid')['name']
             .apply(lambda g: {root(x) for x in g}).to_dict())


def roots_of(field, family):
    '''Only designations of the same mission are comparable: an entry
    from another mission is not an identifier of this catalogue.'''
    return {root(x) for x in str(field).split(';')
            if root(x) and family.match(x.strip())}


cxo['XM=F23'] = [compare_set(acc_roots[i], roots_of(b, CXO_FAMILY))
                 for i, b in zip(cxo['ID'], cxo['F23_CXO'])]
cxo['XM=SIMBAD'] = [compare_set(acc_roots[i], roots_of(b, CXO_FAMILY))
                    for i, b in zip(cxo['ID'], cxo['Chandra_simbad'])]
cxo['sep_F23_XM'] = [
    crossmatch.haversine_distance(r1, d1, r2, d2)
    if pd.notna(r1) and pd.notna(r2) else np.nan
    for r1, d1, r2, d2 in zip(cxo['F23_CXO_RA'], cxo['F23_CXO_DE'],
                              cxo['ra_decimal'], cxo['dec_decimal'])]

for col in ('XM=F23', 'XM=SIMBAD'):
    v = cxo[col].value_counts()
    print(f"{col:11s} same name {v.get('yes', 0):3d}  redesignated "
          f"{v.get('renamed', 0):3d}  differ {v.get('NO', 0):3d}"
          f"  no identifier {v.get('-', 0):3d}")
_nz = [v for v in f23_x['F23_CXO'] if v]
print("\nF23 'Chandra' entries that are not CSC names:",
      sum(1 for v in _nz if not v.startswith('2CXO')), "of", len(_nz),
      sorted({v.split()[0] for v in _nz if not v.startswith('2CXO')}))
cxo.to_csv("Output/xray_id_comparison_chandra.csv", index=False)
cxo.loc[(cxo['XM=F23'] == 'NO') | (cxo['XM=SIMBAD'] == 'NO'),
        ['ID', 'XM_CXO', 'F23_CXO', 'Chandra_simbad', 'XM=F23', 'XM=SIMBAD',
         'sep_F23_XM']]
""")

md("### *XMM-Newton*")

code(r"""
bestx = (xmm_dm2[xmm_dm2['accepted']]
         .sort_values('dM2_xray').groupby('ID').first().reset_index()
         .rename(columns={'IAUNAME': 'XM_XMM'}))

xm = (bestx[['ID', 'XM_XMM', 'SRCID', 'dM2_xray', 'xmm_ra', 'xmm_de']]
      .merge(f23_x[['ID', 'F23_XMM', 'F23_XMM_RA', 'F23_XMM_DE']], on='ID', how='left')
      .merge(simbad[['ID', 'XMM_simbad']], on='ID', how='left')
      .merge(sample[['ID', 'Notes']], on='ID', how='left'))
for c in ('F23_XMM', 'XMM_simbad', 'Notes'):
    xm[c] = xm[c].fillna('')

acc_roots_x = (xmm_dm2[xmm_dm2['accepted']].groupby('ID')['IAUNAME']
               .apply(lambda g: {root(x) for x in g}).to_dict())

xm['XM=F23'] = [compare_set(acc_roots_x[i], roots_of(b, XMM_FAMILY))
                for i, b in zip(xm['ID'], xm['F23_XMM'])]
xm['XM=SIMBAD'] = [compare_set(acc_roots_x[i], roots_of(b, XMM_FAMILY))
                   for i, b in zip(xm['ID'], xm['XMM_simbad'])]
xm['sep_F23_XM'] = [
    crossmatch.haversine_distance(r1, d1, r2, d2) if pd.notna(r1) else np.nan
    for r1, d1, r2, d2 in zip(xm['F23_XMM_RA'], xm['F23_XMM_DE'],
                              xm['xmm_ra'], xm['xmm_de'])]

for col in ('XM=F23', 'XM=SIMBAD'):
    v = xm[col].value_counts()
    print(f"{col:11s} same name {v.get('yes', 0):3d}  redesignated "
          f"{v.get('renamed', 0):3d}  differ {v.get('NO', 0):3d}"
          f"  no identifier {v.get('-', 0):3d}")
_r = xm['XM=F23'] == 'renamed'
print("\nredesignated against F23:", int(_r.sum()),
      "| separation of the two catalogue positions: median %.2f'' max %.2f''"
      % (xm.loc[_r, 'sep_F23_XM'].median(), xm.loc[_r, 'sep_F23_XM'].max()))
xm.to_csv("Output/xray_id_comparison_xmm.csv", index=False)
xm.loc[(xm['XM=F23'] == 'NO') | (xm['XM=SIMBAD'] == 'NO'),
       ['ID', 'XM_XMM', 'F23_XMM', 'XMM_simbad', 'XM=F23', 'XM=SIMBAD',
        'sep_F23_XM']]
""")

md(r"""
## 6. The X-ray identifiers in `Notes`

The same tag the Gaia identifier carries, once per X-ray catalogue, so
that one string says where each counterpart came from and what
corroborates it. `CXO(n)` and `XMM(n)` give the number of accepted
sources when a system has more than one: every source that passes the
test is kept, and the comparison asks whether *any* of them carries the
identifier F23 or SIMBAD gives. A system with no accepted source reads
`none`. The comparison is by name, as in the cells above, and a
redesignation counts as a coincidence.
""")

code(r"""
def xray_tag(prefix, n_sources, vs_f23, vs_simbad):
    if n_sources == 0:
        return f'{prefix}: none'
    head = prefix if n_sources == 1 else f'{prefix}({n_sources})'
    f = {'yes': '=', 'NO': ' differs', '-': ' absent'}[vs_f23]
    s = {'yes': '=', 'NO': ' differs', '-': ' absent'}[vs_simbad]
    if f == '=' and s == '=':
        return f'{head}: XM=F23=SIMBAD'
    if f == '=':
        return f'{head}: XM=F23, SIMBAD' + s
    if s == '=':
        return f'{head}: XM=SIMBAD, F23' + f
    if f == s:
        return f'{head}: XM only, F23+SIMBAD ' + ('differ' if f == ' differs' else 'absent')
    return f'{head}: XM only, F23' + f + ', SIMBAD' + s


n_cxo = csc_dm2[csc_dm2['accepted']].groupby('usrid')['name'].nunique().to_dict()
n_xmm = xmm_dm2[xmm_dm2['accepted']].groupby('ID')['IAUNAME'].nunique().to_dict()
cxo_cmp = cxo.set_index('ID')[['XM=F23', 'XM=SIMBAD']].to_dict('index')
xmm_cmp = xm.set_index('ID')[['XM=F23', 'XM=SIMBAD']].to_dict('index')

notes = {}
for _, r in sample.iterrows():
    i = r['ID']
    c = cxo_cmp.get(i, {'XM=F23': '-', 'XM=SIMBAD': '-'})
    x = xmm_cmp.get(i, {'XM=F23': '-', 'XM=SIMBAD': '-'})
    notes[i] = '; '.join([
        r['Notes'],
        xray_tag('CXO', n_cxo.get(i, 0), agreed(c['XM=F23']), agreed(c['XM=SIMBAD'])),
        xray_tag('XMM', n_xmm.get(i, 0), agreed(x['XM=F23']), agreed(x['XM=SIMBAD'])),
    ])

sample['Notes'] = sample['ID'].map(notes)
sample[['ID', 'GaiaDR3_adopted', 'dM2_gaia_match', 'Notes']] \
    .to_csv("Output/hmxb_sample_notes.csv", index=False)

# The two cross-match tables carry the same string.
csc_out['Notes'] = csc_out['usrid'].map(notes)
xmm_out['Notes'] = xmm_out['ID'].map(notes)
assert csc_out['Notes'].notna().all() and xmm_out['Notes'].notna().all()
csc_out.to_csv("Output/hmxb_csc_crossmatch.csv", index=False)
xmm_out.to_csv("Output/hmxb_5xmm_crossmatch.csv", index=False)
Table.from_pandas(xmm_out).write("Output/hmxb_5xmm_crossmatch.fits", overwrite=True)

print("Chandra table:", len(csc_out), "rows x", csc_out.shape[1], "columns")
print("XMM table    :", len(xmm_out), "rows x", xmm_out.shape[1], "columns")
print()
print(sample['Notes'].str.split('; ').str[-2].value_counts().to_string())
print()
print(sample['Notes'].str.split('; ').str[-1].value_counts().to_string())
""")

code(r"""
sample.loc[~sample['Notes'].str.contains('CXO: none') | ~sample['Notes'].str.contains('XMM: none'),
           ['ID', 'Notes']].head(20)
""")

code(r"""
print("written:", sorted(f for f in os.listdir("Output")
                         if 'csc' in f or 'xmm' in f or 'xray' in f))
""")

nb = {
    "cells": cells,
    "metadata": {
        "kernelspec": {"display_name": "hmxb", "language": "python", "name": "python3"},
        "language_info": {"name": "python", "version": "3.12"},
    },
    "nbformat": 4,
    "nbformat_minor": 5,
}

with open("2_XrayCrossmatch.ipynb", "w") as fh:
    json.dump(nb, fh, indent=1)
print(f"wrote 2_XrayCrossmatch.ipynb with {len(cells)} cells")
