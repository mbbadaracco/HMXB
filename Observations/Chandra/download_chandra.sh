#!/bin/bash

# Initialize HEASoft and CIAO
heainit
ciaoinit

base_path="/home/marina/Doctorado/2026/HMXB_project/Observations/Chandra"
cd "$base_path"
mkdir -p logs # one log per observation

for obsid in $(cat "${base_path}/obsids_chandra.txt"); do
	if [ -n "$(ls -A "$obsid" 2>/dev/null)" ]; then
		echo "Skipping ${obsid}, already downloaded"
		continue
	fi
	echo "Downloading ${obsid}..."
	download_chandra_obsid "$obsid" >> "logs/${obsid}_download.log" 2>&1
done
echo "Done with Chandra data downloading."
