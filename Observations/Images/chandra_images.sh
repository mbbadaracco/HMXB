#!/bin/bash

# Broad-band images and source detection for the Chandra observations of Observations/Images/.
# ACIS gets csc-b, 0.5-7 keV, the CSC detection band. HRC has no energy column and gets csc-w, its 0.1-10 keV wide band, unfiltered.
# Nothing already on disk is made again.

heainit
ciaoinit

base_path="/home/marina/Doctorado/2026/HMXB_project/Observations/Images"
repro_path="/home/marina/Doctorado/2026/HMXB_project/Observations/Chandra_repro"
python_bin="/home/marina/Doctorado/2026/HMXB_project/work/env/bin/python"
scratch_path="${base_path}/scratch"

njobs=12         # observations imaged at the same time
njobs_detect=4   # observations detected at the same time, fewer because each holds its wavelet intermediates at once
min_free_gb=60   # refuse to start without this much room on the filesystem holding the scratch
psfenergy=1.5    # in keV, the energy the PSF map is built at
scales="1.0 2.0 4.0 8.0 16.0"   # wavdetect wavelet scales, in pixels
acis_band="csc-b"; acis_lo=500; acis_hi=7000
hrc_band="csc-w"
# Detection runs on a copy blocked by this. An unblocked ACIS field is 8192x8192 and wavdetect works in double precision, so each intermediate is 537 MB; blocking by 4 divides that by sixteen at the cost of 1.97 arcsec per pixel. Set 1 for full resolution and give it 2.7 GB of scratch per observation.
detect_block=4


# 81 observations cover more than one system of the sample, and their products are the same products: the detection runs on the whole field, not on a region around the system. The first system to make one is linked to by the rest, so the bytes exist once. Two jobs starting the same observation at the same time can still both make it; work/dedup_image_products.py links those afterwards.
link_existing() {
	local name="$1" other
	other=$(ls "${base_path}"/*/"${PWD##*/}"/"$name" 2>/dev/null | grep -vx "${PWD}/${name}" | head -1)
	[ -n "$other" ] || return 1
	ln -f "$other" "$name"
}

cd "$base_path" || exit 1
mkdir -p logs "$scratch_path"

# CIAO writes its working files in ASCDS_WORK_PATH, which ciaoinit sets to /tmp, and /tmp here is a 32 GB tmpfs that twelve concurrent jobs fill.
free_gb=$(df -BG --output=avail "$scratch_path" | tail -1 | tr -dc '0-9')
if [ "${free_gb:-0}" -lt "$min_free_gb" ]; then
	echo "only ${free_gb} GB free on the filesystem holding ${scratch_path}, ${min_free_gb} GB wanted: stopping"
	exit 1
fi

echo "######################################################"
echo " "
echo "# Observations read from ${base_path}/<system>/observations.txt, events from ${repro_path}/<obsid> #"
echo "# One image per observation: ${acis_band} for ACIS, ${hrc_band} for HRC #"
echo "# Detection with wavdetect on that image blocked by ${detect_block} #"
echo "# Scratch in ${scratch_path}/<obsid>, never /tmp; ${free_gb} GB free there #"
echo "# Logs in ${base_path}/logs/<obsid>_chandra-{images,detect}.log #"
echo "# ${njobs} observations imaged at a time, ${njobs_detect} detected #"
echo "# Anything already present is not made again #"
echo " "
echo "######################################################"

# Column 1 of the table is the Chandra observation, "-" on an XMM row.
pairs=""
for table in */observations.txt; do
	system="${table%/observations.txt}"
	[ -d "${system}/Chandra_images" ] || continue
	for obsid in $(awk '!/^#/ && $1 != "-" && $1 != "ChandraObsID" {print $1}' "$table"); do
		pairs="${pairs}${system}|${obsid}"$'\n'
	done
done

echo "Generating Chandra images..."
while IFS='|' read -r system obsid; do
	[ -n "$obsid" ] || continue
	# One subshell per pair, with its own PFILES: concurrent CIAO tools sharing a parameter file corrupt each other.
	(
	evt=$(ls "${repro_path}/${obsid}"/*repro_evt2.fits 2>/dev/null | head -1)
	if [ -z "$evt" ]; then
		echo "Skipping ${obsid} for ${system}: no reprocessed evt2 file"
		exit
	fi
	cd "${system}/Chandra_images" || exit
	mkdir -p "pfiles_${obsid}"
	export PFILES="${PWD}/pfiles_${obsid};${PFILES#*;}"
	export ASCDS_WORK_PATH="${scratch_path}/${obsid}_img"
	export TMPDIR="$ASCDS_WORK_PATH"
	mkdir -p "$ASCDS_WORK_PATH"
	
	instrume=$("$python_bin" -c "from astropy.io import fits; import sys; print(fits.getheader(sys.argv[1], 1).get('INSTRUME', '?'))" "$evt")
	if [[ "$instrume" == HRC* ]]; then img="${obsid}_${hrc_band}_img.fits"; else img="${obsid}_${acis_band}_img.fits"; fi
	if [ -f "$img" ]; then
		echo "Skipping ${img}, it already exists"
	elif link_existing "$img"; then
		echo "Linked ${img} from the system that already imaged this observation"
	elif [[ "$instrume" == HRC* ]]; then
		dmcopy "${evt}[events][bin x=::1,y=::1][IMAGE]" "$img" clobber=no
	else
		dmcopy "${evt}[events][energy=${acis_lo}:${acis_hi}][bin x=::1,y=::1][IMAGE]" "$img" clobber=no
	fi
	rm -rf "pfiles_${obsid}" "$ASCDS_WORK_PATH"
	) >> "logs/${obsid}_chandra-images.log" 2>&1 &
	while [[ $(jobs -rp | wc -l) -ge $njobs ]]; do wait -n; done
done < <(echo "$pairs")
wait
echo "Done generating Chandra images."

echo "Detecting sources in Chandra observations..."
while IFS='|' read -r system obsid; do
	[ -n "$obsid" ] || continue
	(
	cd "${system}/Chandra_images" || exit
	band="$acis_band"
	[ -f "${obsid}_${band}_img.fits" ] || band="$hrc_band"
	img="${obsid}_${band}_img.fits"
	if [ ! -f "$img" ]; then
		echo "No ${acis_band} or ${hrc_band} image: detection skipped for ${obsid}"
		exit
	fi
	
	if [ -f "${obsid}_source_list.fits" ]; then
		echo "Skipping detection for ${obsid}, its source list already exists"
		exit
	fi
	# Same observation, same detection: take the products of whichever system ran it.
	other=$(ls "${base_path}"/*/Chandra_images/"${obsid}_source_list.fits" 2>/dev/null | grep -vx "${PWD}/${obsid}_source_list.fits" | head -1)
	if [ -n "$other" ]; then
		for f in "$(dirname "$other")/${obsid}_"*; do
			[ -f "${PWD}/$(basename "$f")" ] || ln -f "$f" "$(basename "$f")"
		done
		echo "Linked the detection products of ${obsid} from ${other%/*}"
		exit
	fi
	# Every tool below runs with clobber=no, so the leftovers of an attempt that failed would make the retry fail on them instead of on the original cause. The source list is written last, so anything else here is from an unfinished attempt.
	rm -f "${obsid}_psfmap.fits" "${obsid}_src.reg" "${obsid}_source_cell.fits" "${obsid}_imagefile.fits" "${obsid}_bkg.fits"
	mkdir -p "pfiles_det_${obsid}"
	export PFILES="${PWD}/pfiles_det_${obsid};${PFILES#*;}"
	export ASCDS_WORK_PATH="${scratch_path}/${obsid}_det"
	export TMPDIR="$ASCDS_WORK_PATH"
	mkdir -p "$ASCDS_WORK_PATH"
	det="${obsid}_${band}_block${detect_block}_img.fits"
	[ -f "$det" ] || dmcopy "${img}[bin x=::${detect_block},y=::${detect_block}]" "$det" clobber=no
	mkpsfmap "$det" "${obsid}_psfmap.fits" energy=${psfenergy} ecf=0.68 clobber=no
	# Neither mkpsfmap nor wavdetect returns a failing status, so the only way to know is to look for what they should have written.
	if [ ! -f "${obsid}_psfmap.fits" ]; then
		echo "mkpsfmap wrote no PSF map for ${obsid}: detection skipped"
		rm -rf "pfiles_det_${obsid}" "$ASCDS_WORK_PATH"
		exit
	fi
	wavdetect infile="$det" outfile="${obsid}_source_list.fits" scellfile="${obsid}_source_cell.fits" imagefile="${obsid}_imagefile.fits" defnbkgfile="${obsid}_bkg.fits" scales="$scales" regfile="${obsid}_src.reg" expfile=none psffile="${obsid}_psfmap.fits" ellsigma="3.0" clobber=no
	if [ ! -f "${obsid}_source_list.fits" ]; then
		echo "wavdetect wrote no source list for ${obsid}: see above for the reason"
		rm -rf "pfiles_det_${obsid}" "$ASCDS_WORK_PATH"
		exit
	fi
	# The source list as a table anyone can read without CIAO
	"$python_bin" - "${obsid}_source_list.fits" "${obsid}_source_list.csv" <<-'PYEOF'
	import sys
	from astropy.io import fits
	from astropy.table import Table
	with fits.open(sys.argv[1]) as hdul:
	    data = Table(hdul[1].data)
	    names = [n for n in data.colnames if len(data[n].shape) <= 1]
	data[names].to_pandas().to_csv(sys.argv[2], index=False)
	PYEOF
	rm -rf "pfiles_det_${obsid}" "$ASCDS_WORK_PATH"
	) >> "logs/${obsid}_chandra-detect.log" 2>&1 &
	while [[ $(jobs -rp | wc -l) -ge $njobs_detect ]]; do wait -n; done
done < <(echo "$pairs")
wait
rmdir "$scratch_path" 2>/dev/null
echo "Done with detection of sources in Chandra observations."
