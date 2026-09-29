#!/usr/bin/env python
"""Emit Output/hmxb_orbital_parameters.csv, the system parameters F23 reports.

    env/bin/python build_orbital_table.py

One row per HMXB of the working sample, that is per row of
Output/hmxb_crossmatch.csv, which is therefore built first. For each of the
seven quantities F23 measures -- the mass of the compact object, the mass of
the companion, the orbital period, the superorbital period, the eccentricity,
the spin period and the systemic radial velocity -- the table carries the
value, its error and the reference F23 gives for it. The catalogued distance
is deliberately absent: it is not an orbital parameter, and the distance this
work uses is Bailer-Jones (2021).

They live apart from the cross-match table because they answer a different
question. That table is about which source in each catalogue is this system;
this one is about what is known of the orbit. Joining them on ID costs one
line and keeps either readable on its own.

A Notes column records the one reference in the set that cannot be consulted.
No network is needed.
"""

import os

import numpy as np
import pandas as pd

base = "/home/marina/Doctorado/2026/HMXB_project"
work = os.path.join(base, "work")
out = os.path.join(work, "Output", "hmxb_orbital_parameters.csv")

f23_table = os.path.join(work, "Input", "HMXB", "tablea.dat")
sample = os.path.join(work, "Output", "hmxb_crossmatch.csv")
ephem = os.path.join(work, "Input", "Ephemerides", "porb_ephemerides.tsv")

# Byte ranges as the F23 ReadMe declares them, 1-based and both ends included.
# The order is the catalogue's own.
F23_PARAMS = [
    # value column,  error column,    reference column, value, error, reference
    ("Mx_Msun",      "e_Mx_Msun",     "r_Mx",      (156, 161), (163, 167), (169, 198)),
    ("Mopt_Msun",    "e_Mopt_Msun",   "r_Mopt",    (200, 205), (207, 211), (213, 242)),
    ("Porb_d",       "e_Porb_d",      "r_Porb",    (244, 265), (267, 274), (276, 294)),
    ("Psuporb_d",    "e_Psuporb_d",   "r_Psuporb", (296, 304), (306, 313), (315, 333)),
    ("ecc",          "e_ecc",         "r_ecc",     (335, 344), (346, 351), (353, 371)),
    ("Pspin_s",      "e_Pspin_s",     "r_Pspin",   (373, 387), (389, 396), (398, 416)),
    ("RV_kms",       "e_RV_kms",      "r_RV",      (418, 425), (427, 433), (435, 453)),
]
PARAM_COLUMNS = [c for row in F23_PARAMS for c in row[:3]]
PARAM_TRIPLE = {row[0]: (row[1], row[2]) for row in F23_PARAMS}

# arXiv:1503.01087, the source F23 cites for the orbital period, eccentricity
# and systemic radial velocity of two of our systems, was withdrawn from arXiv
# and is unavailable: ADS resolves the bibcode to the abstract page and arXiv
# answers the PDF request with "withdrawn and is unavailable". A value whose
# only source cannot be read is not a measurement anyone can check, so every
# quantity F23 attributes to it is dropped here. The match is on the reference
# column and not on the system name, so anything else ever attributed to that
# bibcode would go with it.
WITHDRAWN = "2015arXiv150301087G"

# What the readable literature says instead, from the two papers in f23_refs/.
# Goossens et al. (2013) measure the orbital period of IGR J18450-0435 under
# its other name AX J1845.0-0433; Bozzo et al. (2024) state that the orbital
# period of AX J1841.0-0536 is still undetermined, so nothing replaces it.
SUPERSEDED = {
    "IGR J18450-0435": {"Porb_d": ("5.7195", "0.0007", "2013MNRAS.434.2182G")},
}

# A value whose cited source does not contain it. Found while looking for an
# epoch: Sidoli & Paizis (2018) Table 1 carries a row for EXO 2030+375 reading
# 46.02 and 0.41, and the paper does not mention RX J2030.5+4751 anywhere. F23
# attributes those two numbers to that paper for RX J2030.5+4751. The values
# are left in place and flagged, not discarded: unlike the withdrawn preprint,
# the source here is readable, and what it says can be checked by anyone.
MISATTRIBUTED = {
    "RX J2030.5+4751":
        "the orbital period and eccentricity F23 gives for this system, 46.02 d "
        "and 0.41, are the values its cited source (Sidoli & Paizis 2018, "
        "Table 1) tabulates for EXO 2030+375; that paper does not mention "
        "RX J2030.5+4751 at all, so both are treated as unverified here",
}

# How each discarded quantity is named in a sentence.
PARAM_WORDS = {"Mx_Msun": "compact-object mass", "Mopt_Msun": "companion mass",
               "Porb_d": "orbital period", "Psuporb_d": "superorbital period",
               "ecc": "eccentricity", "Pspin_s": "spin period",
               "RV_kms": "systemic radial velocity"}

WITHDRAWN_NOTE = {
    "AX J1841.0-0536":
        "Bozzo et al. (2024) report that the orbital period of this system "
        "remains undetermined, and that the 4.7 s spin period of Bamba et al. "
        "(2001) that F23 also tabulates has been questioned and never "
        "confirmed, so no value replaces them",
    "IGR J18450-0435":
        "Goossens et al. (2013) measure Porb = 5.7195 +/- 0.0007 d for the same "
        "system under its other name AX J1845.0-0433, from INTEGRAL/IBIS 18-60 "
        "keV data, and that period is adopted here; they bound the eccentricity "
        "at < 0.37 from the Roche-lobe geometry without measuring it, so no "
        "eccentricity is adopted",
}


def cut(line, span):
    """A byte range as the ReadMe writes it: 1-based, both ends included."""
    return line[span[0] - 1:span[1]].strip()


def wordlist(items):
    """"a", "a and b", "a, b and c"."""
    items = list(items)
    return items[0] if len(items) == 1 else \
        ", ".join(items[:-1]) + " and " + items[-1]


f23 = {}
for line in open(f23_table):
    if not line.strip():
        continue
    row = {}
    for value, error, reference, vb, eb, rb in F23_PARAMS:
        row[value] = cut(line, vb) or None
        row[error] = cut(line, eb) or None
        row[reference] = cut(line, rb) or None
        # The HMXB class, bytes 143-154 of the ReadMe: Be, sg, SFXT or WR. It is not
    # a parameter of the orbit and is not written out, but it decides whether an
    # outburst maximum can be read as a periastron passage; see below.
    row["Class"] = cut(line, (143, 154)) or None
    f23[cut(line, (1, 23))] = row

# The discard, and the replacement where one exists.
discarded = {}
for key, row in f23.items():
    for value, error, reference, _, _, _ in F23_PARAMS:
        if row[reference] != WITHDRAWN:
            continue
        discarded.setdefault(key, []).append(value)
        row[value] = row[error] = row[reference] = None
    for value, (v, e, r) in SUPERSEDED.get(key, {}).items():
        error, reference = PARAM_TRIPLE[value]
        row[value], row[error], row[reference] = v, e, r

t = pd.DataFrame({"ID": pd.read_csv(sample, dtype=str)["ID"].str.strip()})
for column in PARAM_COLUMNS:
    t[column] = t["ID"].map(lambda i, c=column: f23.get(i, {}).get(c))
    if not column.startswith("r_"):
        t[column] = pd.to_numeric(t[column], errors="coerce")

# ---------------------------------------------- the phase zero point and omega
# Neither is in F23: the catalogue gives periods without epochs, and no
# argument of periastron at all. Both are read by hand from the literature into
# Input/Ephemerides/porb_ephemerides.tsv, which carries the quote each value
# was read from, and are merged here on the system name.
eph = pd.read_csv(ephem, sep="\t", comment="#", dtype=str)
eph["ID"] = eph["ID"].str.strip()
EPHEM_COLUMNS = ["MJD_T0", "e_MJD_T0", "T0_kind", "r_T0",
                 "omega_deg", "e_omega_deg", "r_omega", "omega_frame"]
# Read but not emitted: the eccentricity found in the literature fills the F23
# column where F23 has none, so r_ecc always names the source of the value that
# is there; status is a sentence for Notes.
READ_ONLY = ["ecc_lit", "e_ecc_lit", "r_ecc_lit", "status"]
for column in EPHEM_COLUMNS + READ_ONLY:
    t[column] = t["ID"].map(dict(zip(eph["ID"], eph[column])))
    if column not in ("T0_kind", "r_T0", "r_omega", "omega_frame",
                      "r_ecc_lit", "status"):
        t[column] = pd.to_numeric(t[column], errors="coerce")

# An eccentricity read from the literature fills the column only where F23 has
# none; F23's own value is never overwritten.
from_lit = t["ecc"].isna() & t["ecc_lit"].notna()
for a, b in (("ecc", "ecc_lit"), ("e_ecc", "e_ecc_lit"), ("r_ecc", "r_ecc_lit")):
    t.loc[from_lit, a] = t.loc[from_lit, b]

def note(r):
    """Everything about this row a reader would otherwise have to go and find."""
    parts = []
    if r.ID in discarded:
        parts.append("the original reference for the %s of this system is %s "
                     "(arXiv:1503.01087), which is no longer available, having "
                     "been withdrawn from arXiv, so %s discarded here; %s"
                     % (wordlist(PARAM_WORDS[c] for c in discarded[r.ID]),
                        WITHDRAWN,
                        "they are" if len(discarded[r.ID]) > 1 else "it is",
                        WITHDRAWN_NOTE[r.ID]))
    if r.ID in MISATTRIBUTED:
        parts.append(MISATTRIBUTED[r.ID])

    # The epoch, its kind, its reference, the longitude of periastron and its
    # reference are columns of this table, so Notes does not repeat them. It
    # carries only what no column can: why a value is absent, where one came
    # from when the column cannot say, and what is wrong with one that is there.
    if pd.notna(r.ecc) and r.ID in ecc_from_lit:
        parts.append("the eccentricity is not F23's, which reports none for "
                     "this system")
    if r.ID in be_outburst:
        parts.append("the epoch is an outburst maximum and is adopted as the "
                     "periastron passage: F23 classes this system as a Be XB, "
                     "and the Type I outbursts of Be XBs recur at or near "
                     "periastron (Fornasini et al. 2023, Sect. 5.1.2), the "
                     "recurrence being what the orbital period is measured "
                     "from. The identification is only as sharp as the outburst "
                     "is narrow, and that systematic is not included in "
                     "e_MJD_Tper")
    if isinstance(r.status, str) and r.status:
        parts.append(r.status)
    return "; ".join(parts)


# ------------------------------------- the epoch brought to a common convention
# The published epochs are of nine kinds (T0_kind). Periastron passage is the
# one they can all be referred to that is defined for every eccentric orbit and
# is a single instant shared by both bodies, so it is the target.
#
#   Tpi/2      the mean longitude is l = M + omega and M advances linearly, so
#              T_per = T0 - (pi/2 - omega) P / 2pi. Exact, and it needs neither
#              the eccentricity nor a Kepler solve. It uses the omega of the
#              paper that defined Tpi/2, in that paper's own frame, because the
#              two are defined together.
#   conjunction  a conjunction is fixed in true anomaly, not in mean anomaly,
#              so the eccentricity enters. Superior conjunction of the compact
#              object is at nu = pi/2 - omega_compact, inferior at
#              nu = -pi/2 - omega_compact; then E from nu, M = E - e sin E and
#              T_per = T0 - M P / 2pi. Here omega must be the compact object's,
#              so a donor-frame omega is turned by 180 degrees first.
#   periastron identity.
#
#   outburst   only for a Be X-ray binary, and then it is an identity. Fornasini
#              et al. (2023), Sect. 5.1.2 and the caption of their Fig. 5: the
#              Type I outbursts of Be XBs "happen (quasi-)periodically at or
#              near periastron", being fed by the neutron star crossing the
#              decretion disc, and their recurrence is what the orbital period
#              of such a system is measured from in the first place. An epoch of
#              maximum flux for a Be system is therefore a periastron passage.
#              The identification is only as sharp as the outburst is narrow:
#              Fornasini et al. describe Type I outbursts as "usually short
#              lived, lasting only for a small fraction of the orbit", so a
#              systematic of that order rides on these epochs and is not
#              quantified here. It is not applied to the supergiant or SFXT
#              systems, whose outbursts the quoted statement does not cover.
#
# An optical maximum, an ellipsoidal minimum, a spectroscopic phase zero, a
# folding origin, and an outburst maximum of a system that is not a Be XB have
# no fixed relation to the orbit and are not converted.
#
# Two corrections are deliberately left out until the observation times are
# known: the cycle count from the published epoch to the epoch of interest,
# whose error n * sigma_Porb usually dominates everything else, and apsidal
# motion, which makes omega a function of epoch for the few systems where it is
# measured. Both are properties of the epoch one converts *to*, not of the
# conversion, and neither changes the numbers below.
SUPERIOR = ("mid-eclipse", "superior conjunction")
INFERIOR = ("inferior conjunction",)
DRAWS = 20000
rng = np.random.default_rng(20260927)


def mean_anomaly(nu, e):
    """Mean anomaly at true anomaly nu. Direct, no iteration."""
    E = 2 * np.arctan2(np.sqrt(1 - e) * np.sin(nu / 2),
                       np.sqrt(1 + e) * np.cos(nu / 2))
    return E - e * np.sin(E)


def to_periastron(row, n=1):
    """Offset T_per - T0 in days, sampled n times from the quoted errors."""
    P = row.Porb_d + (rng.normal(0, row.e_Porb_d, n)
                      if n > 1 and pd.notna(row.e_Porb_d) else 0)
    w = np.radians(row.omega_deg + (rng.normal(0, row.e_omega_deg, n)
                                    if n > 1 and pd.notna(row.e_omega_deg) else 0))
    if row.T0_kind == "Tpi/2":
        return -(np.pi / 2 - w) * P / (2 * np.pi)
    e = row.ecc + (rng.normal(0, row.e_ecc, n)
                   if n > 1 and pd.notna(row.e_ecc) else 0)
    e = np.clip(e, 0, 0.999)
    if row.omega_frame == "donor":
        w = w + np.pi
    nu = (np.pi / 2 if row.T0_kind in SUPERIOR else -np.pi / 2) - w
    return -mean_anomaly(nu, e) * P / (2 * np.pi)


f23_class = {k: (v.get("Class") or "") for k, v in f23.items()}
be_outburst = set()

tper, e_tper, how = [], [], []
for row in t.itertuples():
    if pd.isna(row.MJD_T0) or pd.isna(row.T0_kind):
        tper.append(np.nan); e_tper.append(np.nan); how.append(None); continue
    if row.T0_kind == "periastron":
        tper.append(row.MJD_T0); e_tper.append(row.e_MJD_T0)
        how.append("as published"); continue
    if row.T0_kind == "outburst" and f23_class.get(row.ID, "").startswith("Be"):
        be_outburst.add(row.ID)
        tper.append(row.MJD_T0); e_tper.append(row.e_MJD_T0)
        how.append("outburst maximum of a Be XB"); continue
    convertible = (row.T0_kind == "Tpi/2" or row.T0_kind in SUPERIOR + INFERIOR)
    if not (convertible and pd.notna(row.omega_deg)
            and (row.T0_kind == "Tpi/2" or pd.notna(row.ecc))):
        tper.append(np.nan); e_tper.append(np.nan); how.append(None); continue
    # The offset is an angle, so it is only defined modulo one orbit; take the
    # periastron nearest the published epoch, |shift| <= P/2. The spread is
    # computed on the unwrapped draws, which cluster around the same cycle.
    raw = float(to_periastron(row, 1))
    shift = (raw + row.Porb_d / 2) % row.Porb_d - row.Porb_d / 2
    spread = float(np.std(to_periastron(row, DRAWS)))
    tper.append(row.MJD_T0 + shift)
    e_tper.append(float(np.hypot(row.e_MJD_T0 if pd.notna(row.e_MJD_T0) else 0,
                                 spread)))
    how.append("from %s" % row.T0_kind)

t["MJD_Tper"] = tper
t["e_MJD_Tper"] = e_tper
t["Tper_from"] = how

ecc_from_lit = set(t.loc[from_lit, "ID"])
t["Notes"] = [note(r) for r in t.itertuples()]

t = t[["ID"] + PARAM_COLUMNS + EPHEM_COLUMNS
      + ["MJD_Tper", "e_MJD_Tper", "Tper_from", "Notes"]].sort_values("ID")
t.to_csv(out, index=False)

print("%d systems -> %s" % (len(t), out))
for column, _, _, _, _, _ in F23_PARAMS:
    print("  %-12s: %d" % (column, t[column].notna().sum()))
print("  with a phase zero point   : %d" % t.MJD_T0.notna().sum())
print("  with a longitude of periastron: %d" % t.omega_deg.notna().sum())
print("  on the common convention (T_per): %d, of which converted: %d"
      % (t.MJD_Tper.notna().sum(),
         int((t.Tper_from.notna() & (t.Tper_from != "as published")).sum())))
print("    outburst maxima adopted as periastron (Be XBs): %d  %s"
      % (len(be_outburst), ", ".join(sorted(be_outburst))))
print("    outburst maxima not adopted, the system not being a Be XB: %s"
      % ", ".join(sorted(r.ID for r in t.itertuples()
                         if r.T0_kind == "outburst" and r.ID not in be_outburst)))
print("  periods still without a zero point: %d"
      % int((t.Porb_d.notna() & t.MJD_T0.isna()).sum()))
print("  quantities discarded with the withdrawn reference:")
for key in sorted(discarded):
    if key in set(t["ID"]):
        print("    %-18s %s" % (key, ", ".join(discarded[key])))
