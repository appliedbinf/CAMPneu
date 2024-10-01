# CAMPneu
***C**omprehensive **A**nalysis of **M**ycoplasma **Pneu**moniae*

CAMPneu is a Nextflow bioinformatic pipeline that is reproducible, scalable, and suitable for a wide range of computation environments. 
While extensible, early drafts of CAMPneu are targeted for Illumina paired-end sequence data with the objectives of 
(1) determining if the specimen belongs to the M. pneumoniae species
(2) classification of the subtype (type1 or type2) of M. pneumoniae
(3) identification of known SNPs conferring macrolide-resistance present within the sample.

### How to Run

1. Installation of the conda package
```
conda install -n campneu -c bioconda -c conda-forge -c appliedbinf campneu 
conda activate campneu 
```

2. Help message
```
CAMPneu.nf --help
```

3. Run command
```
CAMPneu.nf --input <fastq_reads_dir> --output <output_dir>
```
         
Required arguments:     
  --input     Path to the Paired Fastq Reads directory  
  --output    Directory where process outputs are saved          
Optional arguments:  
  --snpFile   Path to the custom SNP bed file
  --help      Print this message and exit
```


### NextFlow script step-by-step workflow:	
1. **Kraken2 Taxonomic Classification:** Classifies input sequences based on a pre-built database.
2. **Quality Control with Fastp:** Profiles and filters reads to ensure high-quality data.
3. **Coverage Assessment with Samtools:** Calculates mean depth to evaluate sequencing coverage.
4. **De Novo Assembly with Unicycler:** Reconstructs microbial genomes without a reference.
5. **ANI Calculation:** Determines the best match by comparing the assembled genomes to reference genomes.
6. **Alignment with Minimap2:** Aligns reads to the best-matched reference genome.
7. **Variant Calling with FreeBayes:** Identifies SNPs and genetic variations against a type 1 reference.
8. **Macrolide-Resistant SNP Identification:** Detects SNPs associated with macrolide resistance

### Cut Off Thresholds ###
The pipeline sets specific thresholds for input paired reads/samples. Any reads or samples that do not meet these thresholds are marked as failed.
1. Kraken2 Percentage of Reads assigned to *M. Pneumoniea* > 90
2. Average Q score > 30
3. Coverage > 10x
4. ANI to reference > 95
5. SNP call quality > 100; Depth > 10

### Required inputs: 
1. Illumina paired-end sequences
2. 23SsnpAnalysis.py: Python script for VCF manipulation and analysis 

### Additional scripts:
**23SsnpAnalysis.py**: This scripts takes the output VCF files generated in the freebayes process and filters it to include snps in the 23S ribosomal RNA in the reference genomes. The scripts requires the co-ordinates of the ribosomal RNA region to subset the VCF file. The script also looks for a set of snps that are provided in a bed file. 

### Outputs:
The scripts generates output directories for each process which have the files generated in the process

#### Process Outputs: 
1. Kraken: kraken reports and kraken summaries for all the paired end reads 
2. fastp: fastp reports and quality filtered paired end reads
3. Coverage_check: samtools coverage report and coverage filtered paired end reads
4. assembly: assembled fasta of the QC filtered samples and empty fasta of the failed samples
5. fastANI: fastANI report
6. bestReference: fastANI report with only the subtyped reference for the sample

#### Summary
1. Sample_reports: Reports for each sample summarizing QC and type information
2. Summary: Report for the entire run w=summarizing which samples have Passed or failed the QC and the SNPs identified for macrolide resistance

