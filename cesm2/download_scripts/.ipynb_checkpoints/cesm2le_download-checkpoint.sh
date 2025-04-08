#!/usr/bin/env bash

# location of download scripts
script_dir="/home/tykukla/aglime-swap-cdr/cesm2/download_scripts/"
outdir="s3://carbonplan-carbon-removal/aglime-swap/cesm2/data/preprocessed/"
# outdir="/home/tykukla/aglime-swap-cdr/cesm2/data/preprocessed/"

for f in "$script_dir"/cesm2le*.py;
do
    echo "now collecting ${f}"
    python3 ${f}
done

# --- start venv
# (may need to activate this in terminal `conda activate cdr-scepter1p0_env`)
conda run -n cdr-scepter1p0_env

# --- move all files to the data directory
# mv b.e21*.nc $outdir
# s5cmd --dry-run mv "${script_dir}/b.e21*.nc" "${outdir}"
s5cmd mv "${script_dir}/b.e21*.nc" "${outdir}"
