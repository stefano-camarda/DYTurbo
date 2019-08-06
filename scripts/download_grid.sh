#!/bin/bash

if [[ -z $1 ]]
then
    "Usage: $0 <jobID>"
    exit
fi

jid=$1

mkdir -p out
cd out
rucio download user.*.${jid}.*.results.root download user.*.${jid}.*.results.txt user.*.*.${jid}.*.log.tgz
