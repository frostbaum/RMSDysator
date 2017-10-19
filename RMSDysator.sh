#!/bin/bash
# Path to program
relpath=RMSDysator_prog
origin=$( dirname "$(readlink -f "$0")" )

abspath=$origin/$relpath

# User input for matching parameters
rmsd=""
while [[ ! $rmsd =~ ^[0-9]*\.?[0-9]+$ ]]; do
    read -p "Enter RMSD LIMIT [Angstrom] " rmsd
done

econ=""
while [[ ! $econ =~ ^[0-9]*\.?[0-9]+$ ]]; do
    read -p "Enter energy conversion constant " econ
done

elim=""
while [[ ! $elim =~ ^[0-9]*\.?[0-9]+$ ]]; do
    read -p "Enter energy cutoff [converted units] " elim
done

anceps="1.0"
while [[ ! $anceps =~ ^[0-9]*\.?[0-9]+$ ]]; do
    read -p "Enter anchor tolerance [Angstrom] " anceps
done

#Make tmp.txt
echo $rmsd $elim $econ > tmp.txt
echo $anceps >> tmp.txt

cnt=0
for file in *.xyz; do
    filenum=( $( echo $file | grep -o -E '[0-9]+') )
    tmpnum[$cnt]=${filenum[@]: -1}
    tmpfile[$cnt]=$file
    echo ${tmpnum[$cnt]} $file >> tmp.txt
    let cnt++
done
nfiles=$cnt
sed -i "1s/^/$nfiles\n/" tmp.txt

#Execute program
$abspath

# User input for path and name of the representative structures
fnamedefault=RMSD_out
read -p "Output file prefix (*_######.xyz) [$fnamedefault]: " fname
fname=${fname:-$fnamedefault}

odirdefault=RMSD_grouping_results
read -p "Output directory [$odirdefault]: " odir
odir=${odir:-$odirdefault}

echo
echo ${odir}/${fname}_######.xyz

# Make appropriate directory if it doesn't exist
[[ -d $odir ]] || mkdir $odir

# read out.txt
mapfile -t outfile < out.txt

for ((i=0;i<${#outfile[@]};i++)); do
  read -ra arr <<<${outfile[$i]}
  outfile[$i]=${arr[4]}
done

# Order and copy structures
for ((k=0;k<${#outfile[@]};k++)); do
  fnum=$( printf  "%06d" $(( k + 1)) )
  for ((i=0;i<$nfiles;i++)); do
    if [ ${tmpnum[$i]} -eq ${outfile[$k]} ]; then
      cp ${tmpfile[$i]} "${odir}/${fname}_${fnum}.xyz"
      break
    fi
  done
done

