#! /bin/bash

set -e

date=$1
mkdir -p $date
for fn in GFA/reduced/v0005_guide/2020${date}/*/guide-*_catalog-00002.fits; do
    echo $fn
    E=$(echo $fn | awk 'BEGIN{FS="/"} {print $5}')
    echo "Exposure $E"
    outdir=${date}/
    echo "Processing $fn"
    echo python3 ../bin/desi_fit_guide_star_coordinates -i $fn -o ${outdir}fm-$E.json --cat ${outdir}fm$E.fits;
    python3 ../bin/desi_fit_guide_star_coordinates -i $fn -o ${outdir}fm-$E.json --cat ${outdir}fm$E.fits;
done
