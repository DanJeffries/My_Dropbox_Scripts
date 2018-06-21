#!bin/bash

## Run from dir contaaining original VCF

## Change headers as plink doesn't like too many "_"s 
sed 's/_/-/g' Diagnostic_SNPs.vcf > Diagnostic_SNPs_altered.vcf


## Convert to plink.raw
~/RAD_programs/PLINK2/plink --vcf Diagnostic_SNPs_altered.vcf --recodeA --allow-extra-chr



echo "All done" 
