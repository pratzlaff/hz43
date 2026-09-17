#! /bin/bash

dir=/data/legs/rpete/flight/hz43

set -eo pipefail

. ~/.bash_aliases
shopt -s expand_aliases

ciao

obsids=$(\grep '^[0-9]' "$dir"/obsids | \grep -i hrc-i | perl -anle 'print $F[0]')

for o in $obsids
do

    evt2=$(\ls "$dir"/data/$o/tg_reprocess/hrcf*${o}N*_evt2.fits)

    read srcx srcy r <<<$(dmlist "$evt2"'[region]' data,raw | \grep circle | perl -anle 'print "$F[2] $F[3] $F[4]"')

    dmstat "$evt2[sky=circle($srcx, $srcy, $r)][col chipx]" verbose=0
    xmin=$(pget dmstat out_min)
    xmax=$(pget dmstat out_max)

    dmstat "$evt2[sky=circle($srcx, $srcy, $r)][col chipy]" verbose=0
    ymin=$(pget dmstat out_min)
    ymax=$(pget dmstat out_max)

    printf "$o: CHIPX=[%.0f, %.0f], CHIPY=[%.0f, %.0f]\n" $xmin $xmax $ymin $ymax

done
