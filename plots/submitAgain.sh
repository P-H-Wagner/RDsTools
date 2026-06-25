#!/bin/bash

# Check argument
if [ $# -ne 1 ]; then
    echo "Usage: $0 <datetime>"
    echo "Example: $0 27_05_2026_13_36_12"
    exit 1
fi

datetime=$1
base="./cmsplots_binned/${datetime}"

if [ ! -d "$base" ]; then
    echo "Directory not found: $base"
    exit 1
fi

for fullpath in "$base"/submitter_sf*.sh; do

    # Handle case where no files match
    [ -e "$fullpath" ] || {
        echo "No submitter scripts found in $base"
        exit 1
    }

    echo "Submitting: $fullpath"
    sbatch -p short --cpus-per-task=8 --mem-per-cpu=4000M --time=60 "$fullpath"

done
