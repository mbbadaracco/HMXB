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
energy_band=(200 12000)
intermezzo=500 # middle value in the energy_band where the condition of PATTERN during the filtering changes
njobs=12 # observations processed at the same time

echo "######################################################"
echo " "
echo "# Default SAS_CCFPATH is ${sas_ccfpath} #"
echo "# The observations to filter are read from ${obsid_file} #"
echo "# The reprocessed files are read and written at ${base_path}/<obsid> #"
echo "# Only FLAG and PATTERN are filtered here: the solar activity is left for later #"
echo "# Only imaging mode is filtered, and an exposure taken with a closed filter is skipped #"
echo "# Every other exposure of every camera is filtered, and one already filtered is not done again #"
echo "# The filter of each exposure is read with ${python_bin} #"
echo "# The output of each observation is appended to ${base_path}/logs/<obsid>_flag-pattern.log #"
echo "# ${njobs} observations are processed at the same time #"
echo "# Change the script for customization. #"
echo " "
echo "######################################################"

obsids=$(cat "$obsid_file")

echo "Correcting EPIC observations for FLAG and PATTERN..."
for obsid in $obsids; do
	# One subshell per observation: it gets its own working directory and
	# its own SAS_ODF and SAS_CCF, which the parallel jobs would otherwise
	# overwrite for each other.
	(
	cd "$obsid" || exit # the observation was not reprocessed
	export SAS_CCF="$base_path/${obsid}/ccf.cif"
	target_sas=$(ls *SUM.SAS)
	export SAS_ODF="$base_path/${obsid}/$target_sas"
	for cam in "${cameras[@]}"; do
		for infile in *"${cam}"*ImagingEvts.ds; do # only imaging mode
			[ -e "$infile" ] || continue # the camera was not used in this observation
			# The name of the event file itself, revolution and observation
			# included, so that every file written here says which exposure
			# of which observation it comes from
			stem="${infile%.ds}"
			if [ -f "flag-pattern_${stem}.fits" ]; then
				echo "Skipping ${stem}, it is already filtered"
				continue
			fi
			# Closed is the optical blocking filter shut, CalClosed is that
			# plus the internal calibration source illuminating the CCDs.
			# Neither observes the sky, and the SOC does not process them
			# either: see the SAS watchout on PPS products for CalClosed
			# observations, www.cosmos.esa.int/web/xmm-newton/sas-watchout
			filter=$($python_bin -c "from astropy.io import fits; import sys; print(fits.getheader(sys.argv[1], 1)['FILTER'])" "$infile")
			if [[ "$filter" == "Closed" || "$filter" == "CalClosed" ]]; then
				echo "Skipping ${stem}, its filter is ${filter}"
				continue
			fi
			if [[ "$cam" == "EPN" ]]; then
				expression="#XMMEA_EP && ((PI>=${energy_band[0]} && PI<=$intermezzo && PATTERN==0) || (PI>$intermezzo && PI<=${energy_band[1]} && PATTERN<=4)) && (FLAG==0)"
			else
				expression="#XMMEA_EM && (PI>=${energy_band[0]} && PI<=${energy_band[1]} && PATTERN<=12)"
			fi
			evselect table="$infile" withfilteredset=Y filteredset="flag-pattern_${stem}.fits" destruct=Y keepfilteroutput=T expression="$expression"
			fkeyprint "$infile"[1] LIVETIME
			fkeyprint flag-pattern_${stem}.fits[1] LIVETIME
		done
	done
	) >> "logs/${obsid}_flag-pattern.log" 2>&1 &
	# Wait for a slot before launching the next observation
	while [[ $(jobs -rp | wc -l) -ge $njobs ]]; do wait -n; done
done
wait
echo "Done correcting EPIC observations for FLAG and PATTERN."
