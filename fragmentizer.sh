#! /bin/bash

# Bash script that launches the protein pipeline for fragment lengths
#between 3 and 20

i=3
UPPER=20
while [ $i -le $UPPER ]
do
  python -W ignore protein.py $i
  ((i++))
done
