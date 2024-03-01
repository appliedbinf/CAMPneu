# CAMPneu
***C**omprehensive **A**nalysis of **M**ycoplasma **Pneu**moniae*

This scripts utilizes a nextflow pipeline that runs multiple existing bioinformatics tools to classify and characterize raw paired reads. The script can classify raw samples based on the reference genomes provided. 

### Command to run the nextflow script:
```
nextflow run CAMPneu.nf
```

### Conda installation of all packages:
```
conda install -n CAMPneu -c bioconda -c conda-forge art bcftools bedtools fastani freebayes minimap2 ncbi-amrfinderplus nextflow raxml samtools snippy spades unicycler 
```

### NextFlow script step-by-step workflow:
1.	De Novo Assembly raw reads: Raw reads are assembles using SPADES denovo assembler
2.	Estimation of best reference using FastANI: The assemblies generated in the steps above are then aligned to the type1 and type 2 references to infer the best reference/most similar type to the sample/isolate assembly
3.	The raw reads are then aligned to the best reference using minimap2 to generate the sam alignment file. 
4.	Sam files are then converted to bamfiles and sorted using samtools.
5.	Freebayes is then used to detect genetic variants with respect to just one reference. 
6.	The snps in the vcf files generated in the freebayes process are then filtered to extract the snps from the 23S ribosomal RNA (the snps mentioned in the paper are in the 23s rRNA region)
7.	Phylogenetic analysis are done through snippy where the core alignment is then visualized using RAXML.

### Required data:
1.	Raw reads
2.	Type 1 and Type 2 reference files

### Additional scripts:
**23SsnpAnalysis.py**: This scripts takes the output VCF files generated in the freebayes process and filters it to include snps in the 23S ribosomal RNA in the reference genomes. The scripts requires the co-ordinates of the ribosomal RNA region to subset the VCF file. The script also looks for a set of snps that are provided in a bed file. 
