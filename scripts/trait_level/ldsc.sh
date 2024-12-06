#!/bin/bash

while getopts "l:s:r:o:" opt; do
  case $opt in
    l) ldsc="$OPTARG"
    ;;
    s) summ="$OPTARG"
    ;;
    r) ref="$OPTARG"
    ;;
    o) output="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

#source activate ldsc


##format
awk -F '\t' '{print $3"\t"$4"\t"$5"\t"$11"\t"$12}' ${summ}.gemma > ${summ}.ldsc
sed -i '1c SNP\tA1\tA2\tZ\tN' ${summ}.ldsc


## summary data for ldsc
gzip ${summ}.ldsc

## heritability
python ${ldsc} --rg ${summ1}.ldsc.gz,${summ2}.ldsc.gz --ref-ld-chr ${ref} --w-ld-chr ${ref} --out ${output}
rm ${summ1}.ldsc.gz

