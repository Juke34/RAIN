#!/usr/bin/env sh

for dir in docker/*
do 
cd ${dir}
docker build -t ${dir} .
cd ..
done
