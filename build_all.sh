#!/usr/bin/env sh

wd=$(pwd)
for dir in docker/*
do 
    cd ${dir}
    imgname=$(echo $dir | rev | cut -d/ -f1 | rev)
    docker build -t ${imgname} .
    cd $wd
done
