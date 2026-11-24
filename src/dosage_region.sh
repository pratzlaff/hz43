#! /bin/bash

# Plots source region CHIP coordinate range on cumulative HRC dosage map.

set -eo pipefail

[ $# -eq 1 ] || {
    echo "Usage: $0 ObsID" >&2
    exit 1
}
obsid="$1"

SCRIPTDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

. ~/.bash_aliases
shopt -s expand_aliases

ciao

evt2=$(ls /data/legs/rpete/flight/hz43/data/$obsid/tg_reprocess/hrcf*$obsid*_evt2.fits)
cumu=/data/legs/rpete/flight/hz43/src/hrci0c_2026-10.fits.gz

region=$(dmlist "$evt2[region]" data,raw | head -1 | perl -ne '/(Circle[^|]+)/; print $1')
dmstat "$evt2[sky=$region][col chip]" 1>/dev/null
read xmin ymin <<<$(pget dmstat out_min | sed 's/,/ /')
read xmax ymax <<<$(pget dmstat out_max | sed 's/,/ /')

# not working: rectangle($xmin,$ymin,$xmax,$ymax)
width=$(echo "$xmax - $xmin" | bc)
height=$(echo "$ymax - $ymin" | bc)
xcen=$(echo "($xmax+$xmin)/2" | bc)
ycen=$(echo "($ymax+$ymin)/2" | bc)

dosage_reg=dosage_$obsid.reg
dosage_img=dosage_area_$obsid.png

rm -f "$dosage_reg" "$dosage_img"

cat <<EOP > $dosage_reg
# Region file format: CIAO version 1.0

box($xcen,$ycen,$width,$height)
EOP

ds9 $cumu -bin factor 16 -zoom .0625 -pan to $xcen $ycen -scale 99.5 -region "$dosage_reg" #-saveimage "$dosage_img" -exit
