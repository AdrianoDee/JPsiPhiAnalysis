#!/usr/bin/env bash

if [[ $# -leq 1 ]] ; then
    echo 'No filename (1) and/or dir(2) given!'
    exit 0
fi


f=$1
d=$2

echo "Building python library with MC file list for ${f} in ${d}"

cd $d

for ff in crab*;
do
  echo $ff
  find $PWD/$ff/*/*/*.root -mmin +30 > $ff.py;
done

cat *.py > $f.py
sed -i -e "s/^/\"file:/g" $f.py
sed -i -e "s/$/\",/g" $f.py
sed -i -e "1s/^/${f} = [ /g" $f.py
sed '$s/,/]/' $f.py

cp $f.py /lustre/home/adrianodif/jpsiphi/2018/CMSSW_10_2_1/src/jpsiphi/jpsiphi/test/mclists/

git add /lustre/home/adrianodif/jpsiphi/2018/CMSSW_10_2_1/src/jpsiphi/jpsiphi/test/mclists/$f.py

git commit -m "Adding ${f} MC list"
git push

cd -
