#! /bin/bash

plink_in=$1
out_file_name=$2


head -n 1 $plink_in > headers.temp

tail -n +2 $plink_in > no_headers.temp

cut -f1 -d' '  no_headers.temp > samples.temp

cut -f7- -d' '  no_headers.temp > data.temp

cut -f7- -d' '  headers.temp > new_headers.temp

paste samples.temp data.temp > sample_data.temp

cat new_headers.temp sample_data.temp > $2.raw

sed -i 's/ /\t/g' $2.raw 

rm *temp




