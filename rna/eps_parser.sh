i=1
while [ $i -le 2 ]
do
  #search for the pattern: 1 to 2 digits, space, 1 to 2 digits space 1 digit.several digits (decimal number) space ubox
  egrep '[0-9]{1,2} [0-9]{1,2} [0-9]\.[0-9]* ubox' sequence${i}_dp.eps | \
        #keep the first three columns, that is the three numbers
        awk '{print $1, $2, $3}' > \
        # output that to file
        out_${i}/bp_probabilities
  
#  SEQUENCE=$(grep -A 1 '/sequence { (\\' sequence${i}_dp.eps | \
#         grep -v '/sequence { (\\')
#  
#  echo ${SEQUENCE::-1} > out_${i}/sequence.txt
#
  ((i++))
done

Rscript plotter.R
