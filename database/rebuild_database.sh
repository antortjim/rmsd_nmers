#! /bin/bash

# Are you currently renovating? > $1
cd top

echo "$(ls -1 | wc -l) files already downloaded :)"
python ../rebuild_database.py

mv obsolete/*/* .

for f in ./*
do
  echo $f
  grep -v '^DBREF' $f > ./temporary
  mv ./temporary $f
done
