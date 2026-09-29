#!/bin/bash

# Download the pipeline spectral products of the XMM-Newton observations that
# hold a spectrum of one of our sources, so that the 5XMM-DR15 fits can be
# repeated on the catalogue's own data. SRSPEC is the source spectrum, BGSPEC
# the background and SRCARF the ancillary response; the pipeline distributes no
# RMF, the spectra name the canned one in their RESPFILE keyword.

base_path="/home/marina/Doctorado/2026/HMXB_project/work/xmm_pps"
obsid_file="/home/marina/Doctorado/2026/HMXB_project/work/Output/xmm_spectra_detections.tsv"
products=("SRSPEC" "BGSPEC" "SRCARF")

mkdir -p "$base_path"
cd "$base_path"
mkdir -p logs # one log per observation

obsids=$(tail -n +2 "$obsid_file" | cut -f2 | sort -u)
echo "$(echo "$obsids" | wc -l) observations to download"

for obsid in $obsids; do
	if [ -n "$(ls -A "$obsid" 2>/dev/null)" ]; then
		echo "Skipping ${obsid}, already downloaded"
		continue
	fi
	echo "Downloading ${obsid}..."
	for product in "${products[@]}"; do
		{
		curl -f -o "${obsid}_${product}.tar" "https://nxsa.esac.esa.int/nxsa-sl/servlet/data-action-aio?obsno=${obsid}&level=PPS&name=${product}&extension=FTZ" &&
		tar -xvf "${obsid}_${product}.tar" &&
		rm -f "${obsid}_${product}.tar"
		} >> "logs/${obsid}_pps.log" 2>&1
		if [ $? -ne 0 ]; then
			echo "FAILED ${obsid} ${product}"
			rm -f "${obsid}_${product}.tar"
		fi
	done
done
echo "Done with XMM-Newton pipeline spectra downloading."
