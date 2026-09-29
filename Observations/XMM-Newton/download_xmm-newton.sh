#!/bin/bash

# Initialize HEASoft and SAS
heainit
sasinit

base_path="/home/marina/Doctorado/2026/HMXB_project/Observations/XMM-Newton"
cd "$base_path"
mkdir -p logs # one log per observation

for obsid in $(cat "${base_path}/obsids_xmm.txt"); do
	if [ -n "$(ls -A "$obsid" 2>/dev/null)" ]; then
		echo "Skipping ${obsid}, already downloaded"
		continue
	fi
	echo "Downloading ${obsid}..."
	mkdir -p "$obsid"
	cd "$obsid"
	{
	curl -f -o "${obsid}.tar" "https://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?obsno=${obsid}&level=ODF" &&
	tar -xvf "${obsid}.tar" &&
	tar -xvf *"${obsid}".TAR
	} >> "${base_path}/logs/${obsid}_download.log" 2>&1
	status=$?
	cd "$base_path"
	# A transfer cut short leaves a truncated tar and half an ODF. Remove
	# the directory, so that the observation is downloaded again the next
	# time instead of being skipped as if it were complete.
	if [ $status -ne 0 ]; then
		echo "FAILED ${obsid}, its incomplete directory is removed"
		rm -rf "$obsid"
	fi
done
echo "Done with XMM-Newton data downloading."
