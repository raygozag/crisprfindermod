#!/bin/bash

set -e
CNT=0
for i in $(seq 0 $1); do
 ./sub.sh ${i}&
 echo ${i}
 CNT=$((CNT+1))
 j=$((CNT%$2))
 if [ ${j} = 0 ]
 then
     wait
 fi
done