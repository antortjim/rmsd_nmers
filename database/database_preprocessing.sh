cd top100H

rm 1etmH

for f in ./*
do
  echo $f
  grep -v '^DBREF' $f > temporary
  mv temporary $f
done

