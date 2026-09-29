#!/bin/bash

# Build a ds9 mosaic of every observation of one HMXB, tiled in order of orbital phase from 0 to ~0.99, Chandra and XMM-Newton together, each tile centred on the source and zoomed so that the two missions are shown at a fixed ratio of sky scale.
#
#   ./ds9_phase_mosaic.sh <system_folder> [zoom_chandra] [--by-mjd] [--dry-run] [--exit]
#   ./ds9_phase_mosaic.sh 1FGL_J1018.6-5856
#   ./ds9_phase_mosaic.sh LS_I+61_303 4 --exit
#   ./ds9_phase_mosaic.sh SS_433 --by-mjd
#
# <system_folder> is a directory of Observations/Images/, the system identifier with underscores for spaces. zoom_chandra is 3 by default and the XMM-Newton and HRC zooms follow from it. --dry-run prints the XPA calls and launches nothing. --exit closes ds9 once the image is saved.
#
# --by-mjd tiles every observation in time order and labels each tile with its MJD, for the systems with no phase to order by; Output/mosaic_systems_no_phase.txt lists them. A system with no phased observation falls back to it on its own.
#
# One ds9 per system, driven through XPA rather than through a single enormous command line, so that every step is one call and a failure says which one.
#
# Pan and zoom, never crop: ds9's crop discards the rest of the image from the display, so a cropped frame cannot be panned or zoomed outside the box and its limits are recomputed over what is left, which flattens a sparse EPIC tile to one colour.
#
# The three plate scales differ, ACIS 0.492 arcsec/pixel, HRC 0.1318 and EPIC binned to 4, so zoom_chandra is the only knob and the other two zooms are derived from it: an HRC tile covers the same sky as an ACIS tile, an XMM-Newton tile field_ratio_xmm times as much. The 10 arcsec circle on the source is the scale of the tile.

base_path="/home/marina/Doctorado/2026/HMXB_project/Observations/Images"
positions="/home/marina/Doctorado/2026/HMXB_project/work/Output/hmxb_sample_positions.vot"
python_bin="/home/marina/Doctorado/2026/HMXB_project/work/env/bin/python"

# The one band each mission is imaged in: 0.5-7 keV for ACIS, 0.1-10 keV unfiltered for HRC, 0.2-12 keV for EPIC.
chandra_band="csc-b"
# HRC has no energy column and gets one unfiltered wide-band image. Three systems have HRC and nothing else: 4U 1954+319, IGR J00370+6122, IGR J18483-0311.
chandra_band_hrc="csc-w"
xmm_band="EP8"
# An observation with all three cameras is drawn as one RGB frame, 161 of the 192; the other 31 have no pn and fall back to a single grey frame. All three channels or none: a two-channel RGB would put a colour on the source meaning "this camera was missing".
xmm_rgb="yes"
xmm_rgb_channels=("EMOS1 red" "EMOS2 green" "EPN blue")
# Which EPIC camera to show when the RGB cannot be built, most sensitive first.
xmm_cameras=("EPN" "EMOS1" "EMOS2")

# zscale collapses on sparse Poisson counts and leaves the tile flat; a high percentile with a square-root stretch is what shows a faint source against an empty detector.
scale_mode="99.5"
# A high percentile on smoothed sparse data saturates all three channels at once and washes the tile out to white, so RGB channels use minmax instead.
scale_mode_rgb="minmax"
scale_fn="sqrt"
# An EPIC image binned to 4 arcsec has very few counts per pixel and needs a small Gaussian to read as a source. Chandra is left unsmoothed, its PSF being narrower than a display pixel.
smooth_xmm="yes"
smooth_xmm_radius=2
# HRC is shown blocked. Unblocked it is 32768x32768 and 2.1 GB, and smoothing one put ds9 at 25 GB and 98 per cent of a core with nothing drawn. chandra_images.sh already writes the blocked copy its detection uses, 8192x8192 and 129 MB, whose 0.527 arcsec pixel is within seven per cent of the ACIS one. The block factor is read from the file name.
hrc_use_blocked="yes"
# Blocking by 4 sums sixteen pixels, so an HRC frame ends with counts per pixel like an ACIS one, which is not smoothed either. A full-resolution HRC image is never smoothed whatever this says.
smooth_hrc="no"
smooth_hrc_radius=2
cmap="blue"
# HRC gets its own colormap: its field is background where an ACIS field is empty, a million counts against a few thousand.
cmap_hrc="heat"
marker_radius_arcsec=10            # radius of the circle on the HMXB position
marker_radius="${marker_radius_arcsec}\""
zoom_chandra=3          # zoom of a Chandra tile; raise it to make the field smaller and the source larger
# How much more sky an XMM-Newton tile covers than a Chandra one. The EPIC PSF is about six arcsec against half an arcsec on axis for ACIS, so at equal field the EPIC tile is a handful of blown-up squares.
field_ratio_xmm=3
# Not locked: "lock frame wcs" forces every frame onto the sky scale of the current one and throws away the per-mission zoom. The script prints the line to send when they should be tied.
lock_frames="no"
# Off, so each frame keeps its own contrast. Nothing about the stretch is shared: lock scalelimits stays no and each frame computes its own limits, scale scope being local.
lock_colorbar="no"
# Measured here: a ds9 that registers can take over two minutes to reach the name server, and one that fails to reach it never recovers, so the wait is long and a ds9 that has not appeared is replaced.
xpa_wait=180
xpa_attempts=2
# Grid type "analysis" draws lines, axes and numerics inside the image; "publication" puts them in a margin a tile has no room for and drops the numerics. The gaps adapt to the zoom either way.
grid_type="analysis"
grid_color="cyan"
grid_numerics_color="yellow"
# Type analysis draws no frame title, so the phase goes in a text region above the source. In screen pixels, not arcsec: a region sits on the sky, so a fixed angular offset leaves the tile as the zoom rises. 70 px also clears the 10 arcsec circle at the default zoom.
label_offset_px=70
# tile_cols=0 costs every possible number of rows and takes the cheapest: how far a tile is from square, plus empty_cell_cost per empty cell in the last row. Squareness alone would put five frames in one row of five.
# saveimage photographs the window, so a small window gives a small picture. win_w=win_h=0 sizes it to one monitor; the desktop here spans two, so the size comes from the monitor xrandr marks primary and not from the whole display.
win_w=0
win_h=0
monitor="primary"
win_margin=60           # pixels kept back for the window decorations and the panel
# -fullscreen works but takes over the screen; the window sized to the monitor gives the same room without that.
use_fullscreen="no"
tile_cols=0
empty_cell_cost=0.15    # how much one empty cell in the last row is worth against the shape of a tile
scalebar_arcsec=10                 # length of the labelled scale bar drawn in each tile
# Gap under the scale bar, in screen pixels because it separates two drawn things: a Chandra tile is 0.16 arcsec per pixel and an XMM-Newton one 0.49, so one angular gap cannot serve both.
bar_label_gap_px=30

system="$1"
[ -n "$2" ] && [ "$2" != "--dry-run" ] && [ "$2" != "--exit" ] && [ "$2" != "--by-mjd" ] && zoom_chandra="$2"
dry_run=no
close_when_done=no
order_by="phase"
for a in "$@"; do [ "$a" = "--dry-run" ] && dry_run=yes; [ "$a" = "--exit" ] && close_when_done=yes; [ "$a" = "--by-mjd" ] && order_by="mjd"; done

if [ -z "$system" ] || [ ! -f "${base_path}/${system}/observations.txt" ]; then
	echo "usage: $0 <system_folder> [zoom_chandra] [--by-mjd] [--dry-run] [--exit]"
	echo "no such system: ${system}"
	exit 1
fi

# Size the window to one monitor when asked to.
if [ "$win_w" -eq 0 ] || [ "$win_h" -eq 0 ]; then
	read -r mon_w mon_h mon_x mon_y < <(xrandr --listmonitors 2>/dev/null | awk -v want="$monitor" '
		NR > 1 {
			name = $NF
			geo = $3
			primary = ($2 ~ /\*/)
			split(geo, a, "+"); split(a[1], b, "x")
			gsub(/\/.*/, "", b[1]); gsub(/\/.*/, "", b[2])
			if ((want == "primary" && primary) || want == name) { print b[1], b[2], a[2], a[3]; exit }
		}')
	if [ -z "$mon_w" ]; then
		# No xrandr, or no monitor of that name: fall back to the whole display, then to a sane default.
		read -r mon_w mon_h < <(xdpyinfo 2>/dev/null | awk '/dimensions:/ {split($2, a, "x"); print a[1], a[2]; exit}')
		mon_x=0; mon_y=0
	fi
	win_w=${mon_w:-1400}
	win_h=${mon_h:-1000}
	win_w=$(( win_w - win_margin ))
	win_h=$(( win_h - win_margin ))
	win_x=${mon_x:-0}
	win_y=${mon_y:-0}
fi

cd "$base_path" || exit 1
table="${system}/observations.txt"
id=$(head -1 "$table" | sed 's/^# *//')

# The position is our adopted Gaia counterpart, not the X-ray one: it is the best position we have for the system itself.
read -r ra dec < <("$python_bin" -c "
import sys
from astropy.table import Table
t = Table.read('${positions}')
row = t[[str(x).strip() == '''${id}''' for x in t['ID']]]
print('%.6f %.6f' % (row['RA_ICRS'][0], row['DE_ICRS'][0]) if len(row) else '')
")
if [ -z "$ra" ]; then
	echo "no position for '${id}' in ${positions}"
	exit 1
fi

bar_offset_arcsec=20       # how far below and to the left of the source the scale bar sits, on Chandra
bar_offset_arcsec_xmm=20   # and on XMM-Newton. The same angular offset in both puts the bar right next to the source in an XMM tile, which covers field_ratio_xmm times more sky; scaling it by that ratio instead pushed it out towards the corner.
# What one arcsecond is worth on each detector.
acis_scale=0.492
epic_scale=4
# The HRC pixel, 0.1318 arcsec from CDELT1 of the images, is 3.73 times finer than ACIS, so zoom_hrc is derived to give an HRC tile the same field as a Chandra one.
hrc_scale=0.1318
acis_per_arcsec=$(awk "BEGIN {printf \"%.2f\", 1 / ${acis_scale}}")
epic_per_arcsec=$(awk "BEGIN {printf \"%.2f\", 1 / ${epic_scale}}")
hrc_per_arcsec=$(awk "BEGIN {printf \"%.2f\", 1 / ${hrc_scale}}")
# A tile of W screen pixels at zoom z shows W/z detector pixels, that is W*scale/z arcsec. Asking the XMM field to be field_ratio_xmm times the Chandra one gives zoom_xmm = zoom_chandra * epic_scale / (acis_scale * field_ratio_xmm).
zoom_xmm=$(awk "BEGIN {printf \"%.3f\", ${zoom_chandra} * ${epic_scale} / (${acis_scale} * ${field_ratio_xmm})}")
# The field a tile covers depends on the size of the tile, which depends on the window and the number of columns, so it is reported in arcsec per screen pixel instead: that number is exact whatever the layout.
arcsec_per_px_chandra=$(awk "BEGIN {printf \"%.3f\", ${acis_scale} / ${zoom_chandra}}")
arcsec_per_px_xmm=$(awk "BEGIN {printf \"%.3f\", ${epic_scale} / ${zoom_xmm}}")
zoom_hrc=$(awk "BEGIN {printf \"%.3f\", ${zoom_chandra} * ${hrc_scale} / ${acis_scale}}")
# The gap under the scale bar, the same number of screen pixels on either mission. An HRC tile needs no third value: zoom_hrc is derived to give it the arcsec per screen pixel of a Chandra tile.
bar_label_gap_c=$(awk "BEGIN {printf \"%.4f\", ${bar_label_gap_px} * ${arcsec_per_px_chandra}}")
bar_label_gap_x=$(awk "BEGIN {printf \"%.4f\", ${bar_label_gap_px} * ${arcsec_per_px_xmm}}")
# And the phase label, likewise the same number of screen pixels above the source on either mission.
label_offset_c=$(awk "BEGIN {printf \"%.4f\", ${label_offset_px} * ${arcsec_per_px_chandra}}")
label_offset_x=$(awk "BEGIN {printf \"%.4f\", ${label_offset_px} * ${arcsec_per_px_xmm}}")

# The scale bar is a plain line with a label rather than a ds9 ruler. One arcsecond of right ascension is 1/3600/cos(dec) degrees, which is why the declination enters. The regions themselves are commented out below.
read -r bar_c_x0 bar_c_x1 bar_c_y bar_c_lx bar_c_ly bar_x_x0 bar_x_x1 bar_x_y bar_x_lx bar_x_ly label_y_c label_y_x < <("$python_bin" -c "
import math
ra, dec = ${ra}, ${dec}
cosd = math.cos(math.radians(dec))
out = []
for field, gap in ((${bar_offset_arcsec}, ${bar_label_gap_c}), (${bar_offset_arcsec_xmm}, ${bar_label_gap_x})):
    off = field / 3600.0
    x0 = ra + off / cosd
    x1 = x0 - (${scalebar_arcsec} / 3600.0) / cosd
    y = dec - off
    out += ['%.8f' % x0, '%.8f' % x1, '%.8f' % y, '%.8f' % ((x0 + x1) / 2), '%.8f' % (y - gap / 3600.0)]
out.append('%.8f' % (dec + ${label_offset_c} / 3600.0))
out.append('%.8f' % (dec + ${label_offset_x} / 3600.0))
print(' '.join(out))
")

# In phase order the unphased rows are skipped; in --by-mjd every row is taken, ordered by the mid-point of the observation. The first field is what the tiles are ordered by and what is written on them.
n_phased=$(awk '!/^#/ && $1 != "ChandraObsID" && NF == 4' "$table" | wc -l)
n_unphased=$(awk '!/^#/ && $1 != "ChandraObsID" && NF == 3' "$table" | wc -l)
# A system with nothing to order by phase would otherwise refuse to draw observations it does have.
if [ "$order_by" = phase ] && [ "$n_phased" -eq 0 ]; then
	echo "${id} has no observation with an orbital phase: ordering by MJD instead"
	order_by="mjd"
fi
if [ "$order_by" = mjd ]; then
	rows=$(awk '!/^#/ && $1 != "ChandraObsID" && NF >= 3 {print $3, $1, $2, $3}' "$table" | sort -n)
	n_rows=$(echo "$rows" | grep -c . )
	n_skipped=0
else
	rows=$(awk '!/^#/ && $1 != "ChandraObsID" && NF == 4 {print $4, $1, $2, $3}' "$table" | sort -n)
	n_rows=$(echo "$rows" | grep -c . )
	n_skipped="$n_unphased"
fi

echo "######################################################"
echo " "
echo "# ${id} at RA=${ra} Dec=${dec}, from ${positions} #"
if [ "$order_by" = mjd ]; then
	echo "# ${n_rows} observations tiled in time order, labelled with the MJD; ${n_phased} of them carry a phase #"
else
	echo "# ${n_rows} observations with a phase, tiled in phase order; ${n_skipped} without one, skipped #"
fi
echo "# Chandra shown in ${chandra_band}, or ${chandra_band_hrc} for an HRC observation, XMM-Newton in ${xmm_band} #"
echo "# An XMM-Newton observation with all three cameras is one RGB frame: EMOS1 red, EMOS2 green, EPN blue #"
echo "# Without all three it falls back to one grey frame, camera preference ${xmm_cameras[*]} #"
echo "# Tiles centred on the source, Chandra at zoom ${zoom_chandra}, XMM-Newton at zoom ${zoom_xmm}, HRC at ${zoom_hrc} unblocked #"
echo "# An HRC observation is shown blocked when the blocked copy exists: 129 MB instead of 2.1 GB, and its zoom follows the block factor #"
echo "# That is ${arcsec_per_px_chandra} arcsec per screen pixel on Chandra and ${arcsec_per_px_xmm} on XMM-Newton, a 1:${field_ratio_xmm} scale #"
echo "# Nothing is cropped: zoom out and pan past the edges of any frame freely #"
#echo "# Each tile carries a ${scalebar_arcsec} arcsec scale bar at the lower left, labelled just below it #"
echo "# The phase is written ${label_offset_px} screen pixels above the source, so it is inside the tile at any zoom #"
echo "# The ${marker_radius_arcsec} arcsec circle on the source is the scale of the tile; the scale bar is commented out #"
echo "# 1 arcsec is ${acis_per_arcsec} ACIS pixels, ${hrc_per_arcsec} HRC pixels, ${epic_per_arcsec} EPIC image pixels #"
echo "# Scale ${scale_fn}/${scale_mode}, colormap ${cmap} and ${cmap_hrc} for HRC, a ${marker_radius} circle on the source #"
echo "# The mosaic is saved to ${base_path}/${system}/${system}_phase_mosaic.png #"
echo "# The ds9 window is ${win_w}x${win_h} at +${win_x:-0}+${win_y:-0}, and that is what the saved png shows #"
echo "# Change the script for customization. #"
echo " "
echo "######################################################"

title="mosaic_${system}"
ds9() { if [ "$dry_run" = yes ]; then echo "xpaset -p ${title} $*"; else xpaset -p "$title" "$@"; fi; }

if [ "$dry_run" = no ]; then
	# XPA addresses a program by name through a name server, xpans, which listens on port 14285. ds9 opens its own socket and announces itself to that server; with no server running the announcement goes nowhere, ds9 works perfectly well and is simply unreachable by name, which is the "did not register" failure. ds9 does not start the server here, so start it if it is not up. It is a daemon: one is enough for every ds9 on the machine, and starting a second is harmless.
	# A process called xpans existing is not the same as a name server answering, and pgrep cannot tell the two apart: a wedged server holds the port, replies to nothing, and makes the guard skip the restart that would fix it. Ask the server a question instead.
	if ! timeout 5 xpaget xpans >/dev/null 2>&1; then
		echo "the XPA name server was not answering: restarting xpans"
		pkill -x xpans 2>/dev/null
		xpans >/dev/null 2>&1 &
		sleep 2
	fi
	registered=no
	for attempt in $(seq "$xpa_attempts"); do
		if [ "$use_fullscreen" = yes ]; then
			command ds9 -title "$title" -fullscreen &
		else
			command ds9 -title "$title" -geometry "${win_w}x${win_h}+${win_x:-0}+${win_y:-0}" &
		fi
		ds9_pid=$!
		for i in $(seq "$xpa_wait"); do [ "$(xpaaccess -n "$title" 2>/dev/null)" = "1" ] && break; sleep 1; done
		if [ "$(xpaaccess -n "$title" 2>/dev/null)" = "1" ]; then
			registered=yes
			break
		fi
		echo "ds9 did not register with XPA after ${xpa_wait}s, attempt ${attempt} of ${xpa_attempts}: replacing it"
		kill $ds9_pid 2>/dev/null
		sleep 2
	done
	if [ "$registered" != yes ]; then
		echo "no ds9 registered with XPA after ${xpa_attempts} attempts of ${xpa_wait}s."
		echo "  ds9 itself is probably fine: check with pgrep -x ds9 and look for its window. What fails is the connection it makes to the name server, not the program."
		echo "  Check the server with xpaget xpans, which lists what is registered and answers even when nothing is; if it hangs or errors, pkill -x xpans and start a fresh one in the background."
		echo "  Check whether ds9 ever reached it: ss -tnp | grep 14285 shows the connection while a registered ds9 is up."
		echo "  Raise xpa_wait or xpa_attempts in this script if the machine is loaded."
		exit 1
	fi
fi

# Global look: no panner or magnifier, they waste room in a mosaic, and a coordinate grid so the axes carry RA and Dec.
ds9 view panner no
ds9 view magnifier no
ds9 view buttons yes
ds9 view colorbar yes
ds9 frame delete all

frame=0
missing=0
while read -r key cid xid mjd; do
	# The tile is labelled with whatever the mosaic is ordered by, the phase or the MJD, so that the number on the picture and the sequence it sits in are the same quantity.
	tile_text="$key"
	# Reset per row, not per branch: a Chandra row that inherited the previous observation's channel list would be drawn as that observation's RGB instead of as itself.
	img=""
	rgb_files=""
	cmap_this="$cmap"
	if [ "$cid" != "-" ]; then
		img="${system}/Chandra_images/${cid}_${chandra_band}_img.fits"
		label="Chandra ${cid}"
		zoom="$zoom_chandra"
		smooth="no"
		smooth_radius="$smooth_hrc_radius"
		bar_x0="$bar_c_x0"; bar_x1="$bar_c_x1"; bar_y="$bar_c_y"; bar_lx="$bar_c_lx"; bar_ly="$bar_c_ly"; label_y="$label_y_c"
		if [ ! -f "$img" ]; then
			# An HRC observation: its own band and its own zoom, the HRC sky pixel not being the ACIS one. The blocked copy is preferred over the full-resolution image, and the zoom follows the block factor read from its name so that the tile covers the same field either way.
			hrc_blocked=$(ls "${system}/Chandra_images/${cid}_${chandra_band_hrc}"_block*_img.fits 2>/dev/null | head -1)
			if [ "$hrc_use_blocked" = yes ] && [ -n "$hrc_blocked" ]; then
				img="$hrc_blocked"
				hrc_block=$(basename "$img" | sed 's/.*_block\([0-9]*\)_img\.fits/\1/')
				label="Chandra/HRC ${cid} blocked ${hrc_block}"
				zoom=$(awk "BEGIN {printf \"%.3f\", ${zoom_chandra} * ${hrc_scale} * ${hrc_block} / ${acis_scale}}")
				smooth="$smooth_hrc"
				cmap_this="$cmap_hrc"
			elif [ -f "${system}/Chandra_images/${cid}_${chandra_band_hrc}_img.fits" ]; then
				img="${system}/Chandra_images/${cid}_${chandra_band_hrc}_img.fits"
				label="Chandra/HRC ${cid}"
				zoom="$zoom_hrc"
				smooth="no"
				cmap_this="$cmap_hrc"
			fi
		fi
	else
		if [ "$xmm_rgb" = yes ]; then
			# All three channels, or none: a two-channel RGB would put a colour on the source that means "this camera was missing".
			for pair in "${xmm_rgb_channels[@]}"; do
				set -- $pair
				cand=$(ls "${system}/XMM-Newton_images/"*"_$1_"*"_${xmm_band}_img.fits" 2>/dev/null | grep -- "_${xid}_" | head -1)
				[ -z "$cand" ] && { rgb_files=""; break; }
				rgb_files="${rgb_files}${cand}|$2 "
			done
		fi
		if [ -n "$rgb_files" ]; then
			img="${rgb_files%%|*}"
			label="${xid}"
		else
			for cam in "${xmm_cameras[@]}"; do
				cand=$(ls "${system}/XMM-Newton_images/"*"_${cam}_"*"_${xmm_band}_img.fits" 2>/dev/null | grep -- "_${xid}_" | head -1)
				[ -n "$cand" ] && { img="$cand"; label="${xid}"; break; }
			done
		fi
		zoom="$zoom_xmm"
		smooth="$smooth_xmm"
		smooth_radius="$smooth_xmm_radius"
		bar_x0="$bar_x_x0"; bar_x1="$bar_x_x1"; bar_y="$bar_x_y"; bar_lx="$bar_x_lx"; bar_ly="$bar_x_ly"; label_y="$label_y_x"
	fi
	if [ -z "$img" ] || [ ! -f "$img" ]; then
		echo "missing image for ${label} at ${order_by} ${key}: skipped"
		missing=$((missing + 1))
		continue
	fi
	frame=$((frame + 1))
	if [ -n "$rgb_files" ]; then
		# An RGB frame holds one image per channel. The stretch is set per channel, after the channel is selected, so each camera is scaled on its own data; a colormap does not apply to an RGB frame.
		ds9 frame new rgb
		for entry in $rgb_files; do
			f="${entry%%|*}"
			ch="${entry##*|}"
			ds9 rgb channel "$ch"
			ds9 fits "$(readlink -f "$f")"
			ds9 scale "$scale_fn"
			ds9 scale scope local
			ds9 scale mode "$scale_mode_rgb"
		done
		ds9 rgb lock scale no
		ds9 rgb lock colorbar no
	else
		ds9 frame new
		ds9 fits "$(readlink -f "$img")"
		ds9 scale "$scale_fn"
		ds9 scale scope local
		ds9 scale mode "$scale_mode"
		ds9 cmap "$cmap_this"
	fi
	ds9 pan to "$ra" "$dec" wcs fk5
	ds9 zoom to "$zoom"
	# Smoothing is inherited by a new frame, so it has to be turned off explicitly or one smoothed XMM frame smooths every frame after it.
	if [ "$smooth" = yes ]; then
		ds9 smooth function gaussian
		ds9 smooth radius "$smooth_radius"
		ds9 smooth yes
	else
		ds9 smooth no
	fi
	ds9 grid yes
	ds9 grid type "$grid_type"
	ds9 grid axes type interior
	ds9 grid numerics type interior
	ds9 grid grid color "$grid_color"
	ds9 grid axes color "$grid_color"
	ds9 grid tickmarks color "$grid_color"
	ds9 grid numerics color "$grid_numerics_color"
	tile_label="fk5; text(${ra},${label_y}) # color=white textangle=0 text={${tile_text}} font=\"helvetica 9 bold roman\""
	marker="fk5; circle(${ra},${dec},${marker_radius}) # color=red width=2"
	#bar="fk5; line(${bar_x0},${bar_y},${bar_x1},${bar_y}) # line=0 0 color=red width=3"
	#bar_label="fk5; text(${bar_lx},${bar_ly}) # color=red textangle=0 text={${scalebar_arcsec}\"} font=\"helvetica 10 bold roman\""
	if [ "$dry_run" = yes ]; then
		echo "echo '${tile_label}' | xpaset ${title} regions"
		echo "echo '${marker}' | xpaset ${title} regions"
		#echo "echo '${bar}' | xpaset ${title} regions"
		#echo "echo '${bar_label}' | xpaset ${title} regions"
	else
		echo "$tile_label" | xpaset "$title" regions
		echo "$marker" | xpaset "$title" regions
		#echo "$bar" | xpaset "$title" regions
		#echo "$bar_label" | xpaset "$title" regions
	fi
	echo "frame ${frame}: ${order_by} ${key}  ${label}  MJD ${mjd}  $(basename "$img")"
done < <(echo "$rows")

if [ "$frame" -eq 0 ]; then
	echo "no images found for ${id}: nothing to tile. Have the image scripts finished?"
	echo "  ${n_phased} of its observations carry a phase and ${n_unphased} do not; --by-mjd tiles both."
	[ "$dry_run" = no ] && kill $ds9_pid 2>/dev/null
	exit 1
fi

# Tile the frames, filling the window as squarely as the count allows.
if [ "$tile_cols" -gt 0 ]; then
	cols="$tile_cols"
else
	cols=$(awk -v n="$frame" -v w="$win_w" -v h="$win_h" -v pen="$empty_cell_cost" 'BEGIN {
		best = -1
		for (r = 1; r <= n; r++) {
			c = int((n + r - 1) / r)                      # the columns r rows need
			if (int((n + c - 1) / c) != r) continue       # those columns would not fill r rows: not a layout
			a = (w / c) / (h / r)                         # the shape of one tile, 1 being square
			d = (a >= 1) ? log(a) : -log(a)               # how far from square, symmetric in wide and tall
			cost = d + pen * (c * r - n)
			if (best < 0 || cost < best) { best = cost; bestc = c }
		}
		print bestc
	}')
fi
rows_n=$(( (frame + cols - 1) / cols ))
ds9 tile mode grid
ds9 tile grid mode manual
ds9 tile grid layout $cols $rows_n
ds9 tile yes
if [ "$lock_frames" = yes ]; then
	ds9 lock frame wcs
else
	ds9 lock frame none
fi
ds9 lock scalelimits no
if [ "$lock_colorbar" = yes ]; then
	ds9 lock colorbar yes
else
	ds9 lock colorbar no
fi
# ds9 8.6 starts in mode "none", where the mouse is inert: no click to select a tile, and the wheel only ever zooms the frame that was already current.
ds9 mode pointer
ds9 frame first

out="${system}/${system}_phase_mosaic.png"
# saveimage photographs the window, so the window has to be on screen and unobscured; there is no headless mode in ds9.
ds9 saveimage png "$(readlink -f "${system}")/$(basename "$out")"
echo "${frame} frames tiled ${cols} across in ${rows_n} rows, ${missing} observations without an image -> ${out}"
echo "To tie the frames together on the sky: xpaset -p ${title} lock frame wcs"
echo "To let them move independently again:  xpaset -p ${title} lock frame none"
echo "To change the stretch on all of them:  xpaset -p ${title} lock scalelimits yes"
[ "$close_when_done" = yes ] && ds9 exit
