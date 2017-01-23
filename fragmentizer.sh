#! /bin/bash

i=3
while [ $i -le 20 ]
do
  python -W ignore protein.py $i
  ((i++))
done



