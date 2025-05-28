#!/usr/bin/env bash

# Run this script with the argument "github_action" to save SIF files in a tar archive
# Run without arguments to simply build all singularity .def files locally

build_mode=$1
wd=$(pwd)

# Create output directory for .sif images
mkdir -p sif_images

# List to store built images
sif_list=( )

# Loop through all .def files in singularity/ subdirectories
for def_file in singularity/*/*.def; do
    # Extract image name from directory
    imgname=$(basename "$(dirname "$def_file")")
    sif_output="sif_images/${imgname}.sif"
    sif_list+=("$sif_output")

    echo ██████████████████▓▒░   Building ${imgname}   ░▒▓██████████████████
    
    # Build with fakeroot if available
    if command -v singularity &> /dev/null; then
        singularity build --fakeroot "$sif_output" "$def_file"
    else
        echo "❌ Singularity not found! Aborting."
        exit 1
    fi
done

# If used in GitHub Actions, create an archive
if [[ ${build_mode} == 'github_action' ]]; then
    echo "Saving SIF images to cache archive..."
    tar -cf singularity-images.tar "${sif_list[@]}"
    echo Archive size: $(stat --printf="%s" singularity-images.tar)
fi