#!/bin/bash

snpcheck_file=$1
snp_file=$2

while read line; do 
    pat=$(echo ${line} | sed 's/^/\^/' | sed 's/$/\t/')
    grep -E -m 1 ${pat} ${snpcheck_file}
done < ${snp_file}
