# CAMPneu
***C**omprehensive **A**nalysis of **M**ycoplasma **Pneu**moniae*

CAMPneu is a Nextflow bioinformatic pipeline that is reproducible, scalable, and suitable for a wide range of computation environments. 
While extensible, early drafts of CAMPneu are targeted for Illumina paired-end sequence data with the objectives of 
(1) determining if the specimen belongs to the M. pneumoniae species
(2) classification of the subtype (type1 or type2) of M. pneumoniae
(3) identification of known SNPs conferring macrolide-resistance and any other AMR-related genes present within the sample.

### Command to run the nextflow script:

CAMPneu has been set up to explicitly run on rosalind/scicomp resources for which prerequites are required to be done. Prior to running the script, the nextflow and conda/mamba modules need to be loaded.

```
module load nextflow
module load miniconda3/20230728 
mamba install -n campneu -c bioconda -c conda-forge -c appliedbinf campneu 
mamba activate campneu 
```

Once NextFlow and Conda are activated, the "CAMPneu' conda package will be created and activated from within the nextflow processes. Then the script can be run.

```
nextflow run CAMPneu.nf --help

Setup Run: Download minikraken database (https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240112.tar.gz) and 
           M.pneumoniae M129 (type 1) and FH (type 2) reference files
Command: nextflow run CAMPneu.nf --download_db

Default Run: Run all the processes in pipeline using the downloaded database and references
Command: nextflow run CAMPneu.nf --input <fastq_reads_dir> --output <output_dir>

Custom Run: Run all the processes in pipeline using the downloaded database and references with a user specified SNP bed file
Command: nextflow run CAMPneu.nf --input <fastq_reads_dir> --output <output_dir> --snpFile <snp_bed_file>
         
Required arguments:     
  --input     Path to the Paired Fastq Reads directory  
  --output    Directory where process outputs are saved          
Optional arguments:  
  --snpFile   Path to the custom SNP bed file
  --help      Print this message and exit
```

### Conda installation of all packages:
```
conda install -c bioconda -c conda-forge appliedbinf::campneu  
```

### NextFlow script step-by-step workflow:	
1. Taxonomic classification with Kraken2 to test for contaminants or erroneous data
2.	Quality score check is done using fastp and reads below the quality score of 30 are dropped
3.	Samtools Coverage is used to fail samples with coverages less than 10x
4.	Raw reads are assembled using Unicycler denovo assembler. 
5.	Estimation of best reference using FastANI: The assemblies generated in the steps above are then aligned to the type1 and type 2 references to infer the best reference/most similar type to the sample/isolate assembly
6.	The raw reads are then aligned to the best reference using minimap2 to generate the sam alignment file. 
7.	Sam files are then converted to bamfiles and sorted using samtools.
8.	Freebayes is then used to detect genetic variants with respect to just one reference. 
9.	The snps in the vcf files generated in the freebayes process are then filtered to extract the snps from the 23S ribosomal RNA (the snps mentioned in the paper are in the 23s rRNA region)

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

