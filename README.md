# CAMPneu
Comprehensive Analysis of Mycoplasma Pneumoniae

This scripts utilizes a nextflow pipeline that runs multiple existing bioinformatics tools to classify and characterize raw paired reads. The script can classify raw samples based on the reference genomes provided. 

Conda installation of all packages:
"""
conda install -n CAMPneu -c bioconda -c conda-forge art bcftools bedtools fastani freebayes minimap2 ncbi-amrfinderplus nextflow raxml samtools snippy spades unicycler 
"""


The script does the following:
1. Assemble raw reads based on de novo assembly
2. ANI check to classify reads to reference types\
3. Map the raw reads to the best reference
4. Find snps in the isolates
5. Run snippy for core alignment
6. RaXML for phylogeny
