# CAMPneu
***C**omprehensive **A**nalysis of **M**ycoplasma **Pneu**moniae*

CAMPneu is a Nextflow bioinformatic pipeline that is reproducible, scalable, and suitable for a wide range of computation environments. 
While extensible, early drafts of CAMPneu are targeted for Illumina paired-end sequence data with the objectives of 
(1) determining if the specimen belongs to the M. pneumoniae species
(2) classification of the subtype (type1 or type2) of M. pneumoniae
(3) identification of known SNPs conferring macrolide-resistance and any other AMR-related genes present within the sample.

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

### Required inputs: 
1. Illumina paired-end sequences
2. Mpneumoniae Type 1(GCF_000027345.1) and Type 2 (GCF_001272835.1) reference files
3. 23SsnpAnalysis.py: Python script for VCF manipulation and analysis 

### Additional scripts:
**23SsnpAnalysis.py**: This scripts takes the output VCF files generated in the freebayes process and filters it to include snps in the 23S ribosomal RNA in the reference genomes. The scripts requires the co-ordinates of the ribosomal RNA region to subset the VCF file. The script also looks for a set of snps that are provided in a bed file. 
