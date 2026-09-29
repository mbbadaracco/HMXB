"""Emit 1_Crossmatch.ipynb.

The notebook is generated from this file so that the cell contents are
diffable and so that it is obvious which cells were copied verbatim from
~/Doctorado/2026/FotometriaSintetica/0_GaiaPhotometry+.ipynb and which
are new.  Cells marked "VERBATIM (cell N)" are byte-identical to cell N
of that notebook except for the module rename AppendixB -> crossmatch
and the output path.
"""

import json

cells = []


def _lines(src):
    """ipynb stores source as a list of lines, each keeping its newline."""
    return src.strip("\n").splitlines(keepends=True)


def md(src):
    cells.append({"cell_type": "markdown", "metadata": {},
                  "source": _lines(src)})


def code(src):
    cells.append({
        "cell_type": "code", "metadata": {}, "execution_count": None,
        "outputs": [], "source": _lines(src),
    })


def code_from(path):
    """Verbatim code cell read from a file, so nested quotes survive."""
    code(open(path).read())


# ---------------------------------------------------------------- header
md(r"""
# Cross-match of the Fortin et al. (2023) HMXB catalogue

Recovered from `~/Doctorado/2026/FotometriaSintetica/0_GaiaPhotometry+.ipynb`.
Cells flagged **VERBATIM** reproduce that notebook unchanged, apart from
renaming the module `AppendixB` to `crossmatch` and pointing the output
path at `Output/`.

Three things are done here:

1. the duplicate search between Fortin et al. (2023, F23, HMXBs) and
   Fortin et al. (2024, F24, LMXBs), positionally and on identifiers;
2. the search for the Gaia DR3 counterpart of every F23 source, which is
   what fixes the distance and therefore what the 3D `N_H` correction
   needs;
3. a comparison between the Gaia DR3 identifier that F23 reports and the
   one SIMBAD lists for the same object, together with the *Chandra* and
   *XMM-Newton* identifiers SIMBAD carries.

**Difference with respect to the original notebook.** There, every
source found in both F23 and F24 was removed from the HMXB sample. Here
only the cases *confirmed* not to be HMXBs are removed; the ambiguous
ones are kept and flagged, and the flag is carried through the rest of
the analysis.

The cross-match statistic is the one of Pineau et al. (2011),
Appendix A, implemented in `crossmatch.py`.
""")

code(r"""
import os
import time
from collections import defaultdict

import numpy as np
import pandas as pd
from astropy.table import Table
from astroquery.gaia import Gaia
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier

import crossmatch

v = Vizier(columns=["**"], row_limit=-1)

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)
pd.set_option('display.max_colwidth', None)

os.makedirs("Output/gaia_data_Tables", exist_ok=True)

# gaia_data() below sends a LEFT JOIN against gaiadr3.astrophysical_parameters
# for a few hundred source_ids.  The Gaia *synchronous* TAP endpoint has a
# short server-side budget and answers that query with
# "Error 408: Job timeout/aborted"; the asynchronous endpoint is the one
# meant for it.  Rather than edit the copied cells, redirect Gaia.launch_job
# to the async endpoint, with a retry for genuinely transient failures.
# This changes no result: same query, same service, longer budget.
_launch_job_sync = Gaia.launch_job
_launch_job_async = Gaia.launch_job_async

def _launch_job_robust(*args, **kwargs):
    last = None
    for attempt in range(4):
        try:
            return _launch_job_async(*args, **kwargs)
        except Exception as exc:
            last = exc
            print(f"  Gaia TAP async attempt {attempt + 1}/4 failed "
                  f"({type(exc).__name__}: {exc}); retrying")
            time.sleep(15 * (attempt + 1))
    raise last

Gaia.launch_job = _launch_job_robust
""")

md(r"""
## 0. Helper functions

VERBATIM (cell 4 of the original notebook), keeping only the two
functions this notebook uses.
""")

code_from("_helpers_cell.py")

md("## 1. Read F23, F24 data\n\nVERBATIM (cell 8).")

code(r"""
hmxb_data = Table.read("Input/HMXB/tablea.dat",readme="Input/HMXB/ReadMe",format="ascii.cds").to_pandas()
lmxb_data = Table.read("Input/LMXB/lmxbcat.dat",readme="Input/LMXB/ReadMe",format="ascii.cds").to_pandas()

hmxb_data['ID'] = hmxb_data['ID'].str.strip()
hmxb_data['ePos'] = pd.to_numeric(hmxb_data['ePos'], errors='coerce') * 3600

lmxb_data['MainID'] = lmxb_data['MainID'].str.strip()
lmxb_data = lmxb_data.drop([5, 187])
lmxb_data = lmxb_data.reset_index(drop=True)
lmxb_data['PosErr'] = pd.to_numeric(lmxb_data['PosErr'], errors='coerce') * 3600

lmxb_data = lmxb_data.rename(columns={
    'GaiaDist': 'Dist',
    'e_GaiaDist': 'e_Dist',
    'E_GaiaDist': 'E_Dist',
    'Period': 'Per',
    'e_Period': 'e_Per',
    'PosErr': 'ePos',
    'MainID': 'ID'
})
""")

code(r"""
print("F23 (HMXB):", len(hmxb_data), "sources")
print("F24 (LMXB):", len(lmxb_data), "sources after dropping rows 5 and 187")
print("F23 sources carrying a Gaia DR3 identifier:",
      int(hmxb_data['GaiaDR3'].astype(str).str.strip().str.fullmatch(r'\d+').sum()))
""")

md("## 2. Crossmatch between F23 and F24\n\nVERBATIM (cell 10).")

code(r"""
XRB_incommon = []
critical_value_90 = 4.61
beta = 0
umbral = 1e-10
for i in range(len(hmxb_data)):
    h_row = hmxb_data.iloc[i]
    for j in range(len(lmxb_data)):
        l_row = lmxb_data.iloc[j]
        h_dist = crossmatch.haversine_distance(h_row['RAdeg'], h_row['DEdeg'], l_row['RAdeg'], l_row['DEdeg'])
        if not (h_dist < 120):
            continue
        dRA = abs(h_row['RAdeg'] - l_row['RAdeg'])
        dDE = abs(h_row['DEdeg'] - l_row['DEdeg'])
        dErr = abs(h_row['ePos']/3600 - l_row['ePos']/3600)
        dM2 = np.nan
        if dRA > umbral and dDE > umbral and dErr > umbral:
            O11, O12, O21, O22 = crossmatch.matrix_rotation(h_row['ePos'], h_row['ePos'], critical_value_90, beta)
            X11, X12, X21, X22 = crossmatch.matrix_rotation(l_row['ePos'], l_row['ePos'], critical_value_90, beta)
            dM2 = crossmatch.Mahalanobis_distance_squared(O11, O12, O21, O22, X11, X12, X21, X22, h_row['RAdeg'], h_row['DEdeg'], l_row['RAdeg'], l_row['DEdeg'])
        XRB_incommon.append({'ID_Fortin23': h_row['ID'], 'ID_Fortin24': l_row['ID'], 'h_D': h_dist, 'dM2': dM2})
df_incommon = pd.DataFrame(XRB_incommon)
""")

code(r"""
df_incommon.sort_values('h_D')
""")

md(r"""
### Check the IDs

VERBATIM (cell 13). The positional search above only sees pairs closer
than 120''. The identifier search below does not depend on the
coordinates at all, which is what lets it find the two pairs whose F23
and F24 positions disagree by more than a thousand arcseconds.
""")

code(r"""
def clean_columns(df, columns):
    for col in columns:
        df[col] = df[col].astype(str).str.strip().replace({"nan": "", "<NA>": "", "None": ""})

hmxb_data_columns = ["ID", "AGILE", "HEAO", "UHURU4", "ARIEL3", "IGR", "ROSAT", "ROSATF", "FERMI", "SAX", "Swift", "XMM", "Chandra", "2MASS", "GaiaDR3"]
lmxb_data_columns = ["ID", "PopularID", "bestID", "HEAO", "UHURU4", "ARIEL3", "IGR", "2E", "ROSAT", "ROSATF", "FERMI", "SAX", "Swift", "XMM", "Chandra", "2MASS", "GaiaDR3"]

clean_columns(hmxb_data, hmxb_data_columns)
clean_columns(lmxb_data, lmxb_data_columns)

for i, row_h in hmxb_data.iterrows():
    ra_h = row_h["RAdeg"]
    dec_h = row_h["DEdeg"]
    id_h = row_h["ID"]
    for j, row_l in lmxb_data.iterrows():
        ra_l = row_l["RAdeg"]
        dec_l = row_l["DEdeg"]
        id_l = row_l["ID"]
        match_found = False
        for col_h in hmxb_data_columns:
            for col_l in lmxb_data_columns:
                val_h = row_h[col_h]
                val_l = row_l[col_l]
                if val_h and val_h == val_l and val_h != "nan":
                    if not match_found:
                        ang_d = crossmatch.haversine_distance(ra_h, dec_h, ra_l, dec_l)
                        print(f"hmxb_data index {i} (ID: {id_h}) and lmxb_data index {j} (ID: {id_l})")
                        print(f"Angular distance: {ang_d}")
                        print("Matches of IDs:")
                        match_found = True
                    print(f"Match: {val_h} in columns {col_h} (Fortin+23) and {col_l} (Fortin+24)")
        if match_found:
            print('----------------------------------------------------------------------------------')
""")

md("### Summary of the two searches\n\nNEW. Collects what the two loops above printed, so it can be tabulated.")

code(r"""
id_pairs = []
for i, row_h in hmxb_data.iterrows():
    for j, row_l in lmxb_data.iterrows():
        n = sum(1 for ch in hmxb_data_columns for cl in lmxb_data_columns
                if row_h[ch] and row_h[ch] == row_l[cl] and row_h[ch] != "nan")
        if n:
            h_dist = crossmatch.haversine_distance(row_h['RAdeg'], row_h['DEdeg'],
                                                   row_l['RAdeg'], row_l['DEdeg'])
            # Same guard as the positional search: identical positions or
            # identical errors leave dM2 undefined.
            dM2 = np.nan
            if (abs(row_h['RAdeg'] - row_l['RAdeg']) > umbral
                    and abs(row_h['DEdeg'] - row_l['DEdeg']) > umbral
                    and abs(row_h['ePos'] / 3600 - row_l['ePos'] / 3600) > umbral):
                O11, O12, O21, O22 = crossmatch.matrix_rotation(row_h['ePos'], row_h['ePos'],
                                                                critical_value_90, beta)
                X11, X12, X21, X22 = crossmatch.matrix_rotation(row_l['ePos'], row_l['ePos'],
                                                                critical_value_90, beta)
                dM2 = crossmatch.Mahalanobis_distance_squared(
                    O11, O12, O21, O22, X11, X12, X21, X22,
                    row_h['RAdeg'], row_h['DEdeg'], row_l['RAdeg'], row_l['DEdeg'])
            id_pairs.append({
                'ID_Fortin23': row_h['ID'], 'ID_Fortin24': row_l['ID'], 'n_ID_matches': n,
                'h_D': h_dist, 'dM2': dM2})
df_idmatch = pd.DataFrame(id_pairs).sort_values('h_D').reset_index(drop=True)

pos_pairs = set(map(tuple, df_incommon[['ID_Fortin23', 'ID_Fortin24']].values))
df_idmatch['found_positionally'] = [tuple(x) in pos_pairs
                                    for x in df_idmatch[['ID_Fortin23', 'ID_Fortin24']].values]
print("pairs from the positional search (< 120''):", len(df_incommon))
print("pairs from the identifier search          :", len(df_idmatch))
df_idmatch
""")

md(r"""
### The two pairs the positional search cannot reach

VERBATIM (cells 14 and 15). The Mahalanobis distance is enormous in both
cases: the two catalogues place the same object more than 1000'' apart,
with positional errors well below an arcsecond, so no positional
criterion can associate them.
""")

code(r"""
critical_value_90 = 4.61
beta = 0
h = hmxb_data[hmxb_data['ID'] == 'SAX J1324.4-6200'].iloc[0]
l = lmxb_data[lmxb_data['ID'] == '4U 1323-62'].iloc[0]
print(crossmatch.haversine_distance(h['RAdeg'], h['DEdeg'], l['RAdeg'], l['DEdeg']))
O11, O12, O21, O22 = crossmatch.matrix_rotation(h['ePos'], h['ePos'], critical_value_90, beta)
X11, X12, X21, X22 = crossmatch.matrix_rotation(l['ePos'], l['ePos'], critical_value_90, beta)
print(crossmatch.Mahalanobis_distance_squared(O11, O12, O21, O22, X11, X12, X21, X22, h['RAdeg'], h['DEdeg'], l['RAdeg'], l['DEdeg']))
""")

code(r"""
critical_value_90 = 4.61
beta = 0
h = hmxb_data[hmxb_data['ID'] == 'IGR J21347+4737'].iloc[0]
l = lmxb_data[lmxb_data['ID'] == 'V* V1727 Cyg'].iloc[0]
print(crossmatch.haversine_distance(h['RAdeg'], h['DEdeg'], l['RAdeg'], l['DEdeg']))
O11, O12, O21, O22 = crossmatch.matrix_rotation(h['ePos'], h['ePos'], critical_value_90, beta)
X11, X12, X21, X22 = crossmatch.matrix_rotation(l['ePos'], l['ePos'], critical_value_90, beta)
print(crossmatch.Mahalanobis_distance_squared(O11, O12, O21, O22, X11, X12, X21, X22, h['RAdeg'], h['DEdeg'], l['RAdeg'], l['DEdeg']))
""")

md(r"""
## 3. Confirmed duplicates removed, ambiguous ones kept

VERBATIM (cell 17): the two F23 entries confirmed through the literature
not to be HMXBs are dropped, and the two F24 entries that duplicate a
confirmed HMXB are dropped from the LMXB side.
""")

code(r"""
hmxb_data = hmxb_data[(hmxb_data['ID'] != '1E 1740.7-2942') & (hmxb_data['ID'] != 'GRS 1758-258')].reset_index(drop=True)
lmxb_data = lmxb_data[(lmxb_data['ID'] != '2XMM J180112.4-254436') & (lmxb_data['ID'] != 'EXO 1846-031')].reset_index(drop=True)
""")

md(r"""
NEW, and this is where we depart from the original notebook: the
remaining sources that appear in both catalogues are **kept** in the
HMXB sample and flagged, instead of being removed. The flag is a free
text column `Notes` carrying `Ambiguous with LMXB: <F24 ID>`, and it
travels with the source through the rest of the analysis. The same
column later receives the provenance of the Gaia counterpart.
""")

code(r"""
# Every F23 source still present that also appears in F24, by either
# search.  A pair whose F24 side has just been removed is no longer an
# ambiguity: it was resolved in favour of the HMXB classification, so it
# is excluded here.
in_common = pd.concat([
    df_incommon[['ID_Fortin23', 'ID_Fortin24']],
    df_idmatch[['ID_Fortin23', 'ID_Fortin24']],
]).drop_duplicates()
in_common = in_common[in_common['ID_Fortin23'].isin(hmxb_data['ID'])
                      & in_common['ID_Fortin24'].isin(lmxb_data['ID'])]

hmxb_data['ID_Fortin24'] = hmxb_data['ID'].map(
    in_common.groupby('ID_Fortin23')['ID_Fortin24'].apply(lambda s: '; '.join(sorted(set(s)))))
hmxb_data['ID_Fortin24'] = hmxb_data['ID_Fortin24'].fillna('')

# The flag is the text itself, so that it stays readable in the output
# tables without a legend.
hmxb_data['Notes'] = np.where(hmxb_data['ID_Fortin24'] != '',
                              'Ambiguous with LMXB: ' + hmxb_data['ID_Fortin24'], '')

print("F23 sources kept:", len(hmxb_data))
print("of which flagged as ambiguous with F24:", int((hmxb_data['Notes'] != '').sum()))
hmxb_data.loc[hmxb_data['Notes'] != '', ['ID', 'Notes', 'RAdeg', 'DEdeg', 'ePos', 'Class']]
""")

md(r"""
## 4. Search for the Gaia DR3 counterpart

VERBATIM (cells 74 and 89). `Input/OnlyCoords/HMXB_output.tsv` is the
VizieR answer to a 2'' cone search around every F23 position in Gaia DR3
(`I/355`); `Input/OnlyCoords/readme.txt` records the search radius and
the request URL is preserved in the header of the `.tsv` itself. Each
candidate is then accepted if the Mahalanobis distance between the F23
position and the Gaia position satisfies `dM2 < 9.21`, the 99% point of
a chi-squared distribution with two degrees of freedom.
""")

code(r"""
# HMXBs from Fortin et al. 2023
gaia_blocks_hmxb = GAIA_BLOCKS_f(hmxb_data, "RAdeg", "HMXB_output.tsv")
gaia_corr_hmxb = gaia_data(gaia_blocks_hmxb["DR3Name"].drop_duplicates().tolist())[["DR3Name", "RA_DEC_CORR"]]
gaia_blocks_hmxb = gaia_blocks_hmxb.astype({"DR3Name": "string"}).merge(gaia_corr_hmxb.astype({"DR3Name": "string"}), on="DR3Name", how="left")
print("Gaia DR3 candidates inside 2'' of an F23 position:", len(gaia_blocks_hmxb))
""")

code(r"""
if os.path.exists("Output/gaia_data_Tables/hmxbs_gaia_data.csv"):
    # DR3Name likewise: read it as a string, never as a float.
    hmxbs_gaia = pd.read_csv("Output/gaia_data_Tables/hmxbs_gaia_data.csv",
                             dtype={'DR3Name': str, 'ID': str})
else:
    hmxbs_crossmatch = []
    for i, hmxb_row in hmxb_data.iterrows():
        block = gaia_blocks_hmxb[gaia_blocks_hmxb["index"] == i]
        for _, gaia_row in block.iterrows():
            O11, O12, O21, O22 = crossmatch.matrix_rotation(hmxb_row['ePos'], hmxb_row['ePos'], 4.61, 0)
            G11, G12, G21, G22 = ((gaia_row["e_RA_ICRS"] / 1000)**2, gaia_row["RA_DEC_CORR"] * (gaia_row["e_RA_ICRS"] / 1000) * (gaia_row["e_DE_ICRS"] / 1000), gaia_row["RA_DEC_CORR"] * (gaia_row["e_RA_ICRS"] / 1000) * (gaia_row["e_DE_ICRS"] / 1000), (gaia_row["e_DE_ICRS"] / 1000)**2)
            dM2 = crossmatch.Mahalanobis_distance_squared(O11, O12, O21, O22, G11, G12, G21, G22, hmxb_row["RAdeg"],
                                               hmxb_row["DEdeg"], gaia_row["RA_ICRS"], gaia_row["DE_ICRS"])
            if dM2 < 9.21:
                hmxbs_crossmatch.append((int(gaia_row["DR3Name"]), i))
    hmxbs_gaia = gaia_data([dr3 for dr3, _ in hmxbs_crossmatch])
    hmxbs_gaia["ID"] = hmxbs_gaia["DR3Name"].map({dr3: hmxb_data.loc[idx, "ID"] for dr3, idx in hmxbs_crossmatch})
    sources = hmxbs_gaia['DR3Name'].tolist()
    batch_size = 50
    dfs = []
    for i in range(0, len(sources), batch_size):
        batch = sources[i:i+batch_size]
        results = v.query_constraints(catalog="I/352", Source=",".join(str(s) for s in batch))
        if results:
            dfs.append(results[0].to_pandas())
    df_vizier = pd.concat(dfs, ignore_index=True)
    df_vizier = df_vizier.rename(columns={'Source': 'DR3Name'})
    df_vizier = df_vizier[['DR3Name', 'rgeo', 'b_rgeo', 'B_rgeo', 'rpgeo', 'b_rpgeo', 'B_rpgeo']]
    hmxbs_gaia = hmxbs_gaia.merge(df_vizier, on='DR3Name', how='inner')
    hmxbs_gaia.to_csv("Output/gaia_data_Tables/hmxbs_gaia_data.csv", index=False)
    del hmxbs_crossmatch, sources, dfs, df_vizier
""")

md(r"""
### The cross-match result on its own

The acceptance loop of the cell above, repeated here so that the
cross-match result is available independently of the merge with the
Bailer-Jones distances, which is a separate requirement.
""")

code(r"""
xm_rows = []
for i, hmxb_row in hmxb_data.iterrows():
    block = gaia_blocks_hmxb[gaia_blocks_hmxb["index"] == i]
    for _, gaia_row in block.iterrows():
        O11, O12, O21, O22 = crossmatch.matrix_rotation(hmxb_row['ePos'], hmxb_row['ePos'], 4.61, 0)
        G11, G12, G21, G22 = ((gaia_row["e_RA_ICRS"] / 1000)**2, gaia_row["RA_DEC_CORR"] * (gaia_row["e_RA_ICRS"] / 1000) * (gaia_row["e_DE_ICRS"] / 1000), gaia_row["RA_DEC_CORR"] * (gaia_row["e_RA_ICRS"] / 1000) * (gaia_row["e_DE_ICRS"] / 1000), (gaia_row["e_DE_ICRS"] / 1000)**2)
        dM2 = crossmatch.Mahalanobis_distance_squared(O11, O12, O21, O22, G11, G12, G21, G22, hmxb_row["RAdeg"],
                                           hmxb_row["DEdeg"], gaia_row["RA_ICRS"], gaia_row["DE_ICRS"])
        if dM2 < 9.21:
            xm_rows.append({'ID': hmxb_row['ID'], 'XM': str(gaia_row['DR3Name']), 'dM2': dM2,
                            'Notes': hmxb_row['Notes']})

xm = pd.DataFrame(xm_rows)
print("F23 sources searched      :", len(hmxb_data))
print("accepted at dM2 < 9.21    :", xm['ID'].nunique(), f"({len(xm)} Gaia sources)")
print("with a Bailer-Jones dist. :", hmxbs_gaia['ID'].nunique(), "(separate requirement)")
""")

md("## 5. SIMBAD identifiers")

code(r"""
import re
CXO_RE = re.compile(r'^(2CXO|CXO|CXOU|CXOGSG)\b')
XMM_RE = re.compile(r'^([1-4]XMM|XMMU|XMMSL)\b')
_D = {'ID': str, 'SIMBAD': str, 'SIMBAD_DR2': str, 'Chandra_simbad': str, 'XMM_simbad': str}

if os.path.exists("Output/simbad_ids.csv"):
    simbad_ids = pd.read_csv("Output/simbad_ids.csv", dtype=_D).fillna('')
else:
    rows = []
    for name in sorted(xm['ID'].unique()):
        try:
            tab = Simbad.query_objectids(name)
        except Exception as exc:
            tab = None
            print(f"  SIMBAD failed for {name!r}: {type(exc).__name__}")
        idents = [str(x).strip() for x in tab['id']] if tab is not None else []
        g3 = [i.split()[-1] for i in idents if i.startswith('Gaia DR3')]
        g2 = [i.split()[-1] for i in idents if i.startswith('Gaia DR2')]
        rows.append({'ID': name,
                     'SIMBAD': g3[0] if g3 else '',
                     'SIMBAD_DR2': g2[0] if g2 else '',
                     'Chandra_simbad': '; '.join(i for i in idents if CXO_RE.match(i)),
                     'XMM_simbad': '; '.join(i for i in idents if XMM_RE.match(i))})
    simbad_ids = pd.DataFrame(rows)
    simbad_ids.to_csv("Output/simbad_ids.csv", index=False)

print("queried:", len(simbad_ids),
      "| with Gaia DR3:", int(simbad_ids['SIMBAD'].str.len().gt(0).sum()),
      "| with Chandra:", int(simbad_ids['Chandra_simbad'].str.len().gt(0).sum()),
      "| with XMM:", int(simbad_ids['XMM_simbad'].str.len().gt(0).sum()))
""")

md(r"""
## 6. Cross-match against F23 and SIMBAD

The `GaiaDR3` column of F23 is `I19`; a Gaia `source_id` exceeds
$2^{53}$, so it is read as a string straight from the record. `0` is the
ReadMe's "no value".
""")

code(r"""
records = [l for l in open("Input/HMXB/tablea.dat").read().split('\n') if l.strip()]
f23_ids = pd.DataFrame({'ID':  [l[:23].strip()       for l in records],
                        'F23': [l[1521:1540].strip() for l in records]})
f23_ids['F23'] = f23_ids['F23'].replace('0', '')

tab = (xm[['ID', 'XM', 'dM2', 'Notes']]
       .merge(f23_ids, on='ID', how='left')
       .merge(simbad_ids[['ID', 'SIMBAD']], on='ID', how='left')
       .fillna(''))

tab['XM=F23']    = np.where(tab['F23'] == '',    '-', np.where(tab['XM'] == tab['F23'],    'yes', 'NO'))
tab['XM=SIMBAD'] = np.where(tab['SIMBAD'] == '', '-', np.where(tab['XM'] == tab['SIMBAD'], 'yes', 'NO'))
tab['F23=SIMBAD'] = np.where((tab['F23'] == '') | (tab['SIMBAD'] == ''), '-',
                             np.where(tab['F23'] == tab['SIMBAD'], 'yes', 'NO'))

for col in ['XM=F23', 'XM=SIMBAD', 'F23=SIMBAD']:
    v = tab[col].value_counts()
    print(f"{col:11s}  agree {v.get('yes', 0):3d}   differ {v.get('NO', 0):3d}   not comparable {v.get('-', 0):3d}")
tab.to_csv("Output/gaia_id_comparison.csv", index=False)
""")

md("### The rows that disagree")

code(r"""
bad = tab[(tab['XM=F23'] == 'NO') | (tab['XM=SIMBAD'] == 'NO') | (tab['F23=SIMBAD'] == 'NO')]
bad[['ID', 'XM', 'F23', 'SIMBAD', 'XM=F23', 'XM=SIMBAD', 'F23=SIMBAD']]
""")

md("### Sources where one of the three has no identifier")

code(r"""
tab.loc[(tab['F23'] == '') | (tab['SIMBAD'] == ''),
        ['ID', 'XM', 'F23', 'SIMBAD']]
""")

md(r"""
### Adopted identifier

One identifier per source is adopted: the one the cross-match returns
where it is corroborated by SIMBAD, the one the cross-match returns
where it is corroborated by F23, and the cross-match value for
3A 0656-072, where the three sources disagree. All three rules select
the same column, so the adopted identifier is the cross-match one
throughout; what differs between sources is what corroborates it. That
is appended to `Notes` as a short tag, `XM` standing for the
cross-match of Sect. 4.
""")

code(r"""
tab['GaiaDR3_adopted'] = tab['XM']

def _provenance(r):
    # 'differs' and 'absent' are kept apart: they are different statements.
    f23 = {'yes': '=', 'NO': ' differs', '-': ' absent'}[r['XM=F23']]
    sim = {'yes': '=', 'NO': ' differs', '-': ' absent'}[r['XM=SIMBAD']]
    if f23 == '=' and sim == '=':
        return 'Gaia: XM=F23=SIMBAD'
    if f23 == '=':
        return 'Gaia: XM=F23, SIMBAD' + sim
    if sim == '=':
        return 'Gaia: XM=SIMBAD, F23' + f23
    if f23 == sim:
        return 'Gaia: XM only, F23+SIMBAD ' + ('differ' if f23 == ' differs' else 'absent')
    return 'Gaia: XM only, F23' + f23 + ', SIMBAD' + sim

tab['Gaia_prov'] = tab.apply(_provenance, axis=1)
tab['Notes'] = np.where(tab['Notes'] == '', tab['Gaia_prov'],
                        tab['Notes'] + '; ' + tab['Gaia_prov'])
tab = tab.drop(columns='Gaia_prov')

assert (tab['GaiaDR3_adopted'] != '').all()
assert (tab['GaiaDR3_adopted'] == tab['XM']).all()
assert set(tab.loc[tab['Notes'].str.contains('XM only, F23+SIMBAD differ', regex=False), 'ID']) == {'3A 0656-072'}
assert tab['Notes'].str.len().max() < 80

print(tab['Notes'].str.split('; ').str[-1].value_counts().to_string())
print("\nadopted for every source:", len(tab))
tab.loc[~tab['Notes'].str.endswith('XM=F23=SIMBAD') | tab['Notes'].str.startswith('Ambiguous'),
        ['ID', 'GaiaDR3_adopted', 'F23', 'SIMBAD', 'Notes']]
""")

md(r"""
### The final sample, as a VOTable

The 109 systems that have both an accepted Gaia counterpart and a
Bailer-Jones distance, with the position of the adopted Gaia DR3 source
and its uncertainties. This is the table the X-ray catalogue searches
start from, so it is written in VOTable form as well.
""")

code(r"""
from astropy.io.votable import from_table, writeto

pos = hmxbs_gaia[['ID', 'RA_ICRS', 'DE_ICRS', 'e_RA_ICRS', 'e_DE_ICRS']] \
        .sort_values('ID').reset_index(drop=True)
assert pos['ID'].is_unique
assert pos.notna().all().all()

vot_table = Table.from_pandas(pos)
for col, unit in (('RA_ICRS', 'deg'), ('DE_ICRS', 'deg'),
                  ('e_RA_ICRS', 'mas'), ('e_DE_ICRS', 'mas')):
    vot_table[col].unit = unit

votable = from_table(vot_table)
field = votable.get_first_table()
for col, ucd in (('ID', 'meta.id;meta.main'),
                 ('RA_ICRS', 'pos.eq.ra;meta.main'),
                 ('DE_ICRS', 'pos.eq.dec;meta.main'),
                 ('e_RA_ICRS', 'stat.error;pos.eq.ra'),
                 ('e_DE_ICRS', 'stat.error;pos.eq.dec')):
    field.get_field_by_id_or_name(col).ucd = ucd

writeto(votable, "Output/hmxb_sample_positions.vot")
print("rows written:", len(pos))
pos.head()
""")

md("## 7. Save")

code(r"""
hmxb_data.to_csv("Output/f23_flagged.csv", index=False)
tab.to_csv("Output/gaia_id_comparison.csv", index=False)
df_incommon.to_csv("Output/f23_f24_positional.csv", index=False)
df_idmatch.to_csv("Output/f23_f24_identifiers.csv", index=False)
tab.to_csv("Output/gaia_id_comparison.csv", index=False)
tab[['ID', 'GaiaDR3_adopted', 'dM2', 'Notes']] \
   .to_csv("Output/hmxb_gaia_adopted.csv", index=False)
print("written:", sorted(os.listdir("Output")))
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

with open("1_Crossmatch.ipynb", "w") as fh:
    json.dump(nb, fh, indent=1)
print(f"wrote 1_Crossmatch.ipynb with {len(cells)} cells")
