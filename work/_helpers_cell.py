def GAIA_BLOCKS_f(df, ra_string, simbad_output):
    '''
    Reads output from https://vizier.cds.unistra.fr/viz-bin/VizieR-3 where a search around 
    2'' for each input source have been performed, and builds blocks for each one.
    '''
    tol = 1 / 3600
    records = []
    current_constraint_ra = None
    current_index = None
    with open(f"Input/OnlyCoords/{simbad_output}") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#Constraint"):
                value = line.split()[1]
                if "," in value:
                    ra_str = value.split(",")[0]
                    current_constraint_ra = float(ra_str)
                else:
                    current_constraint_ra = float(value)
                diffs = np.abs(df[ra_string].values - current_constraint_ra)
                best = int(np.argmin(diffs))
                current_index = best if diffs[best] < tol else None
                continue
            if line.startswith("#"):
                continue
            if line.startswith("Gaia DR3") and current_index is not None:
                parts = line.split()
                records.append({
                    "index": current_index,
                    "Constraint_RA": current_constraint_ra,
                    "DR3Name": parts[2],
                    "RA_ICRS": float(parts[3]),
                    "DE_ICRS": float(parts[4]),
                    "e_RA_ICRS": float(parts[5]),
                    "e_DE_ICRS": float(parts[6])
                })
    return pd.DataFrame(records)

def gaia_data(list_, chunk_size=1000):
    dfs = []

    for i in range(0, len(list_), chunk_size):
        chunk = list_[i:i + chunk_size]
        ids = tuple(int(x) for x in chunk)

        query = f"""
            SELECT
                g.source_id AS DR3Name,
                g.ra AS RA_ICRS,
                g.dec AS DE_ICRS,
                g.ra_error AS e_RA_ICRS,
                g.dec_error AS e_DE_ICRS,
                g.ra_dec_corr AS RA_DEC_CORR,
                g.distance_gspphot AS Dist,
                g.distance_gspphot_lower AS b_Dist,
                g.distance_gspphot_upper AS B_Dist,
                g.parallax AS parallax,
                g.parallax_error AS parallax_error,
                g.phot_g_mean_mag AS phot_g_mean_mag,
                g.bp_rp AS bp_rp,
                g.pseudocolour AS pseudocolour,
                ap.ew_espels_halpha AS ew_Halpha
            FROM gaiadr3.gaia_source AS g
            LEFT JOIN gaiadr3.astrophysical_parameters AS ap
                ON g.source_id = ap.source_id
            WHERE g.source_id IN {ids}
        """

        job = Gaia.launch_job(query)
        results = job.get_results()
        dfs.append(results.to_pandas())

    return pd.concat(dfs, ignore_index=True)
