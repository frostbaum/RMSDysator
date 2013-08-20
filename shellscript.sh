#!/bin/bash
rmsd=""
while [[ ! $rmsd =~ ^[0-9]*\.?[0-9]+$ ]]; do
    read -p "Enter RMSD LIMIT [Angstrom] " rmsd
done
echo $rmsd > tmp.txt

cnt=0
for file in *.xyz; do
    let cnt++
    filenum=$( [[ $file =~ ([0-9]+).*\.xyz ]] && echo ${BASH_REMATCH[1]} )
    echo $filenum $file >> tmp.txt
done
sed -i "1s/^/$cnt\n/" tmp.txt
