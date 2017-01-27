i=1
while [ $i -le 2 ]
do
  egrep '[0-9]{1,2} [0-9]{1,2} [0-9]\.[0-9]* ubox' sequence${i}_dp.eps | \
        awk '{print $1, $2, $3}' > \
        out_${i}/bp_probabilities
  
  SEQUENCE=$(grep -A 1 '/sequence { (\\' sequence${i}_dp.eps | \
         grep -v '/sequence { (\\')
  
  echo ${SEQUENCE::-1} > out_${i}/sequence.txt

  Rscript parser.R $i
  ((i++))
done
