#!/bin/bash

# Initialize HEASoft and CIAO
heainit
ciaoinit

base_path="/home/marina/Doctorado/2026/HMXB_project/Observations/Chandra_repro"
data_path="/home/marina/Doctorado/2026/HMXB_project/Observations/Chandra"
obsid_file="${data_path}/obsids_chandra.txt"
njobs=12 # observations reprocessed at the same time

cd "$base_path"
mkdir -p logs # one log per observation, since the jobs run at the same time

echo "Beginning reprocessing of Chandra data..."
for obsid in $(cat "$obsid_file"); do
	# One subshell per observation, so that the working directory and the
	# parameter files are its own: concurrent CIAO tools writing the same
	# .par file corrupt each other's parameters.
	(
	if [ -n "$(ls -A "$obsid" 2>/dev/null)" ]; then
		echo "Skipping ${obsid}, it is already reprocessed"
		exit
	fi
	mkdir -p "pfiles_${obsid}"
	export PFILES="${base_path}/pfiles_${obsid};${PFILES#*;}" # private first, system second
	chandra_repro indir="${data_path}/${obsid}" outdir="${obsid}"
	rm -rf "pfiles_${obsid}"
	) >> "logs/${obsid}_repro.log" 2>&1 &
	# Wait for a slot before launching the next observation
	while [[ $(jobs -rp | wc -l) -ge $njobs ]]; do wait -n; done
done
wait
echo "Done with reprocessing Chandra data."
