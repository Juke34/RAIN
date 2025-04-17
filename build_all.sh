#!/usr/bin/env sh

build_mode=$1

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

    echo ██████████████████▓▒░   Building ${imgname}   ░▒▓██████████████████
    
    # Reditools2 does not compile on arm64, force using amd64 compilation
    if [[ $dir =~ "reditools2" ]];then
        if [[ "$arch" == arm* || "$arch" == "aarch64" ]]; then
            echo "Reditools2 does not compile on arm64, force using amd64 compilation"
            docker_arch_option=" --platform linux/amd64"
        fi
    fi

    # Check "github_action" mode to enable cache
    if [[ ${build_mode} == 'github_action' ]]; then
        docker buildx build --cache-from type=local,src=/tmp/.buildx-cache --cache-to type=local,dest=/tmp/.buildx-cache ${docker_arch_option} -t ${imgname} .
    else
        docker build ${docker_arch_option} -t ${imgname} .
    fi
    
    # back to the original working directory
    cd $wd
done
