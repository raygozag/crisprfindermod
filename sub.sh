#!/bin/bash

SID=$RANDOM
ID=$1
CV=$(pwd)
cd /tmp/ramdisk
mkdir ds${SID}
cd ds${SID}
~/Downloads/CRISPRFinder-v3.2.pl ${CV}/sub${ID}.fasta ${CV}/results/
#mv /tmp/ramdisk/${SID}${SID}/* ${CV}/results/

rm -rf /tmp/ramdisk/ds${SID}
