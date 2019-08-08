#!/bin/bash

if [[ -z $1 ]]
then
    echo "usage: $0 <jobID>"
    exit
fi

jid=$1

mkdir -p out
cd out
rucio download user.${USER}.${jid}.*.results.root download user.${USER}.${jid}.'*'.results.txt user.${USER}.'*'.${jid}.*.log.tgz

tag=`cd user.${USER}; ls *.tgz |cut -d . -f3 |uniq`
mkdir -p $tag/root
mkdir -p $tag/txt
mkdir -p $tag/log
mv user.${USER}/*.root $tag/root
mv user.${USER}/*.txt $tag/txt
mv user.${USER}/*.tgz $tag/log

rmdir user.${USER}

cd $tag/log

for file in `ls *.tgz`
do
    tar -xzf $file
    rm $file
done
