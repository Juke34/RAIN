#!/usr/bin/env sh

# get architecture
arch=$(uname -m)
# set architecture to docker buildx
docker_arch_option=""

# save original working directory
wd=$(pwd)
for dir in docker/*
do 
    cd ${dir}
    imgname=$(echo $dir | rev | cut -d/ -f1 | rev)
    
    # Reditools2 does not compile on arm64, force using amd64 compilation
    if [[ $dir =~ "reditools2" ]];then
        if [[ "$arch" == arm* || "$arch" == "aarch64" ]]; then
            echo "Reditools2 does not compile on arm64, force using amd64 compilation"
            docker_arch_option=" --platform linux/amd64"
        fi
    fi
    
    docker build ${docker_arch_option} -t ${imgname} .
    
    # back to the original working directory
    cd $wd
done
