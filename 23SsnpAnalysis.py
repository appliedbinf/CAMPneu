#!/usr/bin/env python3
import subprocess
import sys

vcf = sys.argv[1]
reference = sys.argv[2]
start = sys.argv[3]
end = sys.argv[4]
regions_file = sys.argv[5]

if reference in vcf:
    result = subprocess.Popen(f"bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' {vcf} | head -1 | cut -f 1", shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = result.communicate()
    chrom = stdout.decode().strip()
    filename_all = vcf.replace(".fna.sorted.vcf.gz","all23S.subset.txt")
    filename_reg = vcf.replace(".fna.sorted.vcf.gz","identified.snps.txt")
    subprocess.run(f"bcftools view -r {chrom}:{start}-{end} --no-header {vcf} > {filename_all}", shell=True)
    subprocess.run(f"bcftools view -R {regions_file} --no-header {vcf} > {filename_reg}", shell=True)
    
#
#SNPS to look for are
#snp1 = C217T      NC_000912.1	120272	120273
#snp2 = T1112G     NC_000912.1	121167	121168
#snp3 = A2063G/T/C NC_000912.1	122118	122119
#snp4 = A2064G     NC_000912.1	122119	122120
#snp5 = A2067G     NC_000912.1  122122  122123
#snp6 = A2431G     NC_000912.1	122486	122487
#snp7 = C2611G     NC_000912.1	122666	122667
#snp8 = C2617A     NC_000912.1  122672  122673
#snp9 = 
