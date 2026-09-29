#!/bin/bash

base_path="/home/marina/Doctorado/2026/HMXB_project/Observations/XMM-Newton_repro"
cd "$base_path" # where the reprocessed files and the plots are stored
mkdir -p logs # one log per observation, since the jobs run at the same time

# Initialize HEASoft and SAS
heainit
sasinit

sas_ccfpath="/home/marina/Software/SAS/CCF"
odf_path="/home/marina/Doctorado/2026/HMXB_project/Observations/XMM-Newton"
obsid_file="${odf_path}/obsids_xmm.txt"
python_bin="/home/marina/Doctorado/2026/HMXB_project/work/env/bin/python"
export SAS_CCFPATH=${sas_ccfpath}

cameras=("EMOS1" "EMOS2" "EPN")
timebinsize=100 # seconds, the time resolution of the light curve
njobs=12 # observations processed at the same time

# Default thresholds in counts/s. They are only drawn on the plots, as a
# reference to look at: no filtering is done here.
thr_mos=0.35
thr_pn=0.40

# Single events (PATTERN==0) above 10 keV, one light curve per exposure
expression_mos='#XMMEA_EM && (PI>10000) && (PATTERN==0)'
expression_pn='#XMMEA_EP && (PI>10000 && PI<12000) && (PATTERN==0)'

echo "######################################################"
echo " "
echo "# Default SAS_CCFPATH is ${sas_ccfpath} #"
echo "# The observations to reprocess are read from ${obsid_file} #"
echo "# Their ODFs are read from ${odf_path}/<obsid> #"
echo "# The reprocessed files and the plots are written to ${base_path}/<obsid> #"
echo "# An observation that already has its event files is not reprocessed again #"
echo "# One light curve and one plot are made for each camera and each imaging exposure #"
echo "# The background light curves are plotted with ${python_bin} #"
echo "# The output of each observation is appended to ${base_path}/logs/<obsid>_repro.log #"
echo "# ${njobs} observations are processed at the same time #"
echo "# Change the script for customization. #"
echo " "
echo "######################################################"

obsids=$(cat "$obsid_file")

echo "Beginning reprocessing of XMM-Newton data..."
for obsid in $obsids; do
	# One subshell per observation: it gets its own working directory and
	# its own SAS_ODF and SAS_CCF, which the parallel jobs would otherwise
	# overwrite for each other. The light curves of an observation are
	# extracted and plotted as soon as it is reprocessed.
	(
	# Reprocessing
	mkdir -p "$obsid"
	cd "$obsid" || exit
	if ls *Evts.ds > /dev/null 2>&1; then
		echo "Skipping the reprocessing of ${obsid}, its event files are already there"
	else
		export SAS_ODF="${odf_path}/${obsid}"
		cifbuild # Identifies the necessary calibration files among all CCF for this particular ODF
		export SAS_CCF="$base_path/${obsid}/ccf.cif"
		odfingest # Updates the summary file containing all the information of the ODF
		target_sas=$(ls *SUM.SAS)
		export SAS_ODF="$base_path/${obsid}/$target_sas"
		emproc # Meta-tasks. Default values works for majority of cases.
		epproc
	fi
	# The light curves are extracted with the calibration of the observation,
	# whether it was reprocessed now or in an earlier run
	export SAS_CCF="$base_path/${obsid}/ccf.cif"
	target_sas=$(ls *SUM.SAS)
	export SAS_ODF="$base_path/${obsid}/$target_sas"

	# Background light curves (solar flares), to be inspected by eye
	for cam in "${cameras[@]}"; do
		for infile in *"${cam}"*ImagingEvts.ds; do # only imaging mode
			[ -e "$infile" ] || continue # the camera was not used in this observation
			if [[ "$cam" == "EPN" ]]; then
				expression="$expression_pn"
				threshold="$thr_pn"
			else
				expression="$expression_mos"
				threshold="$thr_mos"
			fi
			# <camera>_<exposure>_<mode>Evts, so that an observation with
			# several exposures of the same camera does not overwrite itself
			tag="${infile#*_*_}"
			tag="${tag%.ds}"
			if [ -f "rate_${tag}.pdf" ]; then
				echo "Skipping ${tag}, its light curve is already plotted"
				continue
			fi
			evselect table="$infile" withrateset=Y rateset="rate_${tag}.fits" maketimecolumn=Y timebinsize=$timebinsize makeratecolumn=Y expression="$expression"
			$python_bin - "rate_${tag}.fits" "rate_${tag}.pdf" "$obsid" "$tag" "$threshold" <<'EOF'
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.io import fits

rateset, plotfile, obsid, tag, threshold = sys.argv[1:6]
with fits.open(rateset) as hdu:
    time = hdu[1].data["TIME"]
    rate = hdu[1].data["RATE"]

plt.figure(figsize=(6,6))
plt.plot(time - time[0], rate, lw=0.8, color="black")
plt.axhline(float(threshold), color="red", ls="--", lw=0.8,
            label=f"default threshold {threshold} counts/s")
plt.xlabel("Time since the beginning of the observation [s]", fontsize=14)
plt.ylabel("Rate [counts/s]", fontsize=14)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.title(f"{obsid} {tag}: single events above 10 keV")
plt.legend()
plt.tight_layout()
plt.savefig(plotfile, dpi=300)
EOF
		done
	done
	) >> "logs/${obsid}_repro.log" 2>&1 &
	# Wait for a slot before launching the next observation
	while [[ $(jobs -rp | wc -l) -ge $njobs ]]; do wait -n; done
done
wait
echo "Done with reprocessing XMM-Newton data and extracting the background light curves."
