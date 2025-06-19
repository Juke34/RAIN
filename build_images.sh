#!/usr/bin/env sh

# Run this script with the argument "github_action" in order to save images in an archive for caching
# Pass no arguments to run the script in "normal" build mode suitable for a local machine

# Read first argument as the "build mode"
# The build is used, for instance, for special build commands
# for Github actions.
build_mode=$1

# get architecture
arch=$(uname -m)
# set architecture to docker buildx
docker_arch_option=""

# save original working directory
wd=$(pwd)

# list of image names
image_list=( )

for dir in docker/*
do 
    cd ${dir}
    imgname=$(echo $dir | rev | cut -d/ -f1 | rev)
    image_list+=(${imgname})

    echo ██████████████████▓▒░   Building ${imgname}   ░▒▓██████████████████
    
    # Reditools2 does not compile on arm64, force using amd64 compilation
    if [[ "$arch" == arm* || "$arch" == "aarch64" ]]; then
        echo "Reditools2 does not compile on arm64, force using amd64 compilation"
        docker_arch_option=" --platform linux/amd64"
    fi

    docker build ${docker_arch_option} -t ${imgname} .
    
    # back to the original working directory
    cd $wd
done

if [[ ${build_mode} == 'github_action' ]]; then
    echo "Saving docker images to cache..."
    docker save ${image_list[@]} -o docker-images.tar
    echo Archive size: $(stat --printf="%s" docker-images.tar)
fi
