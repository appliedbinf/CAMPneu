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
    filename_all = vcf.replace(".fna.sorted.vcf.gz","all23S.subset.vcf")
    filename_reg = vcf.replace(".fna.sorted.vcf.gz","identified.snps.vcf")
    subprocess.run(f"bcftools view -r {chrom}:{start}-{end} --no-header {vcf} > {filename_all}", shell=True)
    subprocess.run(f"bcftools view -R {regions_file} --no-header {vcf} > {filename_reg}", shell=True)
    
#
#SNPS to look for are
#snp1 = C217T
#snp2 = T1112G
#snp3 = A2063G
#snp4 = A2064G
#snp5 = A2431G  
#snp6 = C2611G



