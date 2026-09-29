#!/bin/bash

# Broad-band images and source detection for the XMM-Newton observations of Observations/Images/.
# One EP8 image, 0.2-12 keV, per camera and exposure. Nothing already on disk is made again.

heainit
sasinit

base_path="/home/marina/Doctorado/2026/HMXB_project/Observations/Images"
repro_path="/home/marina/Doctorado/2026/HMXB_project/Observations/XMM-Newton_repro"
python_bin="/home/marina/Doctorado/2026/HMXB_project/work/env/bin/python"
sas_ccfpath="/home/marina/Software/SAS/CCF"
scratch_path="${base_path}/scratch"
export SAS_CCFPATH=${sas_ccfpath}

cameras=("EMOS1" "EMOS2" "EPN")
bin_size=80        # Ensures that 4arcsec=1pixel
njobs=12           # observations imaged at the same time
njobs_detect=4     # observations detected at the same time
min_free_gb=60     # refuse to start without this much room on the filesystem holding the scratch
band="EP8"; band_lo=200; band_hi=12000   # the full EPIC band of 5XMM-DR15; the events were already filtered to it


# 81 observations cover more than one system of the sample, and their products are the same products: the detection runs on the whole field, not on a region around the system. The first system to make one is linked to by the rest, so the bytes exist once. Two jobs starting the same observation at the same time can still both make it; work/dedup_image_products.py links those afterwards.
link_existing() {
	local name="$1" other
	other=$(ls "${base_path}"/*/"${PWD##*/}"/"$name" 2>/dev/null | grep -vx "${PWD}/${name}" | head -1)
	[ -n "$other" ] || return 1
	ln -f "$other" "$name"
}

cd "$base_path" || exit 1
mkdir -p logs "$scratch_path"

# SAS_TMPDIR defaults to /tmp, and /tmp here is a 32 GB tmpfs that a dozen concurrent jobs fill.
free_gb=$(df -BG --output=avail "$scratch_path" | tail -1 | tr -dc '0-9')
if [ "${free_gb:-0}" -lt "$min_free_gb" ]; then
	echo "only ${free_gb} GB free on the filesystem holding ${scratch_path}, ${min_free_gb} GB wanted: stopping"
	exit 1
fi

echo "######################################################"
echo " "
echo "# Observations read from ${base_path}/<system>/observations.txt, events from ${repro_path}/<obsid> #"
echo "# One ${band} image (${band_lo}-${band_hi} eV) per camera and exposure #"
echo "# Detection with edetect_chain on the ${band} images of the three cameras #"
echo "# Scratch in ${scratch_path}/<obsid>, never /tmp; ${free_gb} GB free there #"
echo "# Logs in ${base_path}/logs/<obsid>_xmm-{images,detect}.log #"
echo "# ${njobs} observations imaged at a time, ${njobs_detect} detected #"
echo "# Anything already present is not made again #"
echo " "
echo "######################################################"

# Column 2 of the table is the XMM-Newton observation, "-" on a Chandra row.
pairs=""
for table in */observations.txt; do
	system="${table%/observations.txt}"
	[ -d "${system}/XMM-Newton_images" ] || continue
	for obsid in $(awk '!/^#/ && $2 != "-" && $2 != "XMM-NewtonObsID" {print $2}' "$table"); do
		pairs="${pairs}${system}|${obsid}"$'\n'
	done
done

echo "Generating EPIC images..."
while IFS='|' read -r system obsid; do
	[ -n "$obsid" ] || continue
	# One subshell per pair, with its own SAS_ODF and SAS_CCF, which parallel jobs would otherwise overwrite for each other.
	(
	obsdir="${repro_path}/${obsid}"
	if [ ! -d "$obsdir" ]; then
		echo "Skipping ${obsid} for ${system}: not reprocessed"
		exit
	fi
	export SAS_CCF="${obsdir}/ccf.cif"
	target_sas=$(ls "${obsdir}"/*SUM.SAS 2>/dev/null | head -1)
	if [ -z "$target_sas" ]; then
		echo "Skipping ${obsid} for ${system}: no SUM.SAS"
		exit
	fi
	export SAS_ODF="$target_sas"
	export SAS_TMPDIR="${scratch_path}/${obsid}_img"
	export TMPDIR="$SAS_TMPDIR"
	mkdir -p "$SAS_TMPDIR"
	cd "${system}/XMM-Newton_images" || exit
	# One image per filtered event file, that is per camera and per exposure, as the reprocessing and the filtering keep them.
	for cam in "${cameras[@]}"; do
		for evt in "${obsdir}"/flag-pattern_*"${cam}"*ImagingEvts.fits; do
			[ -e "$evt" ] || continue
			stem=$(basename "$evt" .fits)
			stem="${stem#flag-pattern_}"
			
			img="${stem}_${band}_img.fits"
			if [ -f "$img" ]; then
				echo "Skipping ${img}, it already exists"
				continue
			fi
			if link_existing "$img"; then
				echo "Linked ${img} from the system that already imaged this observation"
				continue
			fi
			evselect table="$evt" xcolumn=X ycolumn=Y imagebinning=binSize ximagebinsize=$bin_size yimagebinsize=$bin_size withimageset=true imageset="$img" expression="(PI in [${band_lo}:${band_hi}])"
		done
	done
	rm -rf "$SAS_TMPDIR"
	) >> "logs/${obsid}_xmm-images.log" 2>&1 &
	while [[ $(jobs -rp | wc -l) -ge $njobs ]]; do wait -n; done
done < <(echo "$pairs")
wait
echo "Done generating EPIC images."

echo "Detecting sources in EPIC observations..."
while IFS='|' read -r system obsid; do
	[ -n "$obsid" ] || continue
	(
	obsdir="${repro_path}/${obsid}"
	export SAS_CCF="${obsdir}/ccf.cif"
	target_sas=$(ls "${obsdir}"/*SUM.SAS 2>/dev/null | head -1)
	[ -n "$target_sas" ] || exit
	export SAS_ODF="$target_sas"
	export SAS_TMPDIR="${scratch_path}/${obsid}_det"
	export TMPDIR="$SAS_TMPDIR"
	mkdir -p "$SAS_TMPDIR"
	att_file=$(ls "${obsdir}"/*AttHk.ds 2>/dev/null | head -1)
	cd "${system}/XMM-Newton_images" || exit
	
	if [ -f "${obsid}_emllist.fits" ]; then
		echo "Skipping detection for ${obsid}, its source list already exists"
		rm -rf "$SAS_TMPDIR"
		exit
	fi
	# Same observation, same detection: take the products of whichever system ran it.
	other=$(ls "${base_path}"/*/XMM-Newton_images/"${obsid}_emllist.fits" 2>/dev/null | grep -vx "${PWD}/${obsid}_emllist.fits" | head -1)
	if [ -n "$other" ]; then
		for f in "$(dirname "$other")"/*"${obsid}"*; do
			[ -f "${PWD}/$(basename "$f")" ] || ln -f "$f" "$(basename "$f")"
		done
		echo "Linked the detection products of ${obsid} from ${other%/*}"
		rm -rf "$SAS_TMPDIR"
		exit
	fi
	# One image per camera, the first exposure of each, and the event list it came from.
	detect_images=""
	detect_events=""
	for cam in "${cameras[@]}"; do
		for evt in "${obsdir}"/flag-pattern_*"${cam}"*ImagingEvts.fits; do
			[ -e "$evt" ] || continue
			stem=$(basename "$evt" .fits)
			stem="${stem#flag-pattern_}"
			img="${stem}_${band}_img.fits"
			if [ -f "$img" ] && [[ "$detect_images" != *"_${cam}_"* ]]; then
				detect_images="${detect_images} ${img}"
				detect_events="${detect_events} ${evt}"
			fi
		done
	done
	if [ -z "$detect_images" ]; then
		echo "No ${band} image for ${obsid}: detection skipped"
	elif [ -z "$att_file" ]; then
		echo "No attitude file for ${obsid}: detection skipped"
	else
		n=$(echo $detect_images | wc -w)
		pimin=$(printf "${band_lo} %.0s" $(seq $n))
		pimax=$(printf "${band_hi} %.0s" $(seq $n))
		ecf=$(printf "1.0 %.0s" $(seq $n))
		edetect_chain imagesets="${detect_images# }" eventsets="${detect_events# }" attitudeset="$att_file" pimin="$pimin" pimax="$pimax" ecf="$ecf" eboxl_list="${obsid}_eboxlist_l.fits" eboxm_list="${obsid}_eboxlist_m.fits" eml_list="${obsid}_emllist.fits" esp_nsplinenodes=16 esen_mlmin=10 eml_ecut=0.68 eml_fitextent=yes psfmodel=ellbeta
		# The source list as a table anyone can read without SAS
		"$python_bin" - "${obsid}_emllist.fits" "${obsid}_emllist.csv" <<-'PYEOF'
		import sys
		from astropy.io import fits
		from astropy.table import Table
		with fits.open(sys.argv[1]) as hdul:
		    data = Table(hdul[1].data)
		    names = [n for n in data.colnames if len(data[n].shape) <= 1]
		data[names].to_pandas().to_csv(sys.argv[2], index=False)
		PYEOF
	fi
	rm -rf "$SAS_TMPDIR"
	) >> "logs/${obsid}_xmm-detect.log" 2>&1 &
	while [[ $(jobs -rp | wc -l) -ge $njobs_detect ]]; do wait -n; done
done < <(echo "$pairs")
wait
rmdir "$scratch_path" 2>/dev/null
echo "Done detecting sources in EPIC observations."
