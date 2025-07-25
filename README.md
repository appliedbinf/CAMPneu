# CAMPneu
***C**omprehensive **A**nalysis of **M**ycoplasma **Pneu**moniae*

CAMPneu is a Nextflow bioinformatic pipeline that is reproducible, scalable, and suitable for a wide range of computation environments. 
While extensible, CAMPneu is currently designed for Illumina paired-end sequence data with the objectives of: 
1. Determining if the specimen belongs to the M. pneumoniae species
2. Classification of the subtype (type1 or type2) of M. pneumoniae
3. Identification of known SNPs conferring macrolide-resistance present within the sample
4. *Under evaluation: Identification of SNPs identified from near-neighbors that confer resistance to tetracycline and fluoroquinolone antibiotics*

**System Requirements:**

CAMPneu requires systems to have the following installed/available:
1. Nextflow (to be used with the singularity profile)
2. Conda OR Singularity

CAMPneu is designed to work with both Conda and Singularity containers, offering flexibility and reproducibility in computational environments.

**CONDA:**

Conda excels at managing dependencies and creating isolated environments. Conda is also easy to use across different operating systems and is ideal for setting up reproducible environments on local machines.

1. Installation using Conda
```
conda install -n campneu -c bioconda -c conda-forge -c appliedbinf campneu 
conda activate campneu
```

2. Run command
```
CAMPneu.nf --input <fastq_reads_dir> --output <output_dir> -profile conda
```

3. Help message
```
CAMPneu.nf --help
```

**SINGULARITY:**

Singularity ensures consistency and portability across systems and is tailored for high-performance computing (HPC) environments offering enhanced efficiency.

The conda installed version of CAMPneu can also be run using singularity but if the user does not have access to conda, they can clone the git repository (nextflow is required for this approach).

1. Git installation
```
git clone https://github.com/appliedbinf/CAMPneu.git
```

2. Run the program from the project repository:
```
nextflow run CAMPneu.nf --input <fastq_reads_dir> --output <output_dir> -profile singularity
```

3. Help message
```
nextflow run CAMPneu.nf --help
```

**Script Input Requirements**

Required arguments:   
```
  --input     Path to the Paired Fastq Reads directory  
  --output    Directory where process outputs are saved          
```
Optional arguments:
``` 
  --help      Print this message and exit
```

### NextFlow script step-by-step workflow:	
1. **Kraken2 Taxonomic Classification:** Classifies input sequences based on a pre-built database.
2. **Quality Control with Fastp:** Profiles and filters reads to ensure high-quality data.
3. **Coverage Assessment with Samtools:** Calculates mean depth to evaluate sequencing coverage.
4. **De Novo Assembly with Unicycler:** Reconstructs microbial genomes without a reference. This is only used to screen for AMR genes with AMRFinder.
5. **ANI Calculation:** Determines the best match by comparing the assembled genomes to reference genomes.
6. **Alignment with Minimap2:** Aligns reads to the type1 reference genome.
7. **Variant Calling with FreeBayes:** Identifies SNPs and genetic variations against a type 1 reference.
8. **Macrolide-Resistant SNP Identification:** Detects SNPs in the 23S rRNA gene associated with macrolide resistance
9. **Tetracycline-Resistant SNP Identification:** Detects SNPs in the 16S rRNA gene associated with tetracycline resistance.
10. **Quinolone-resistant SNP Identification:** Uses snpEff to identify specific mutations known to confer fluoroquinolone-resistance in closely related species.

### Cut Off Thresholds ###
The pipeline sets specific thresholds for input paired reads/samples. Any reads or samples that do not meet these thresholds are marked as failed.
1. Kraken2 Percentage of Reads assigned to *M. Pneumoniea* > 90
2. Average Q score > 30
3. Coverage > 10x
4. ANI to reference > 95
5. SNP call quality > 100; Depth > 10

### AMR Detection
All genome locations are based on:  
*M. pneumoniae* M129; NCBI Accession: NC_000912.1

#### Macrolide SNPS
CAMPNeu will report a potentially resistance strain if a variant is detected to ANY base in the below positions.
| Gene Target | Gene Location | Genome Location | Variant | Note|
| ----------- | ------------- | ------------ | ------- | --- |
| 23S | 217 | 120273 | C > T | M129 has a G here |
| 23S | 1112 | 121168 | T > G | M129 has a C here | 
| 23S | 2063 | 122119 | A > G/T/C | |
| 23S | 2064 | 122120 | A > G | |
| 23S | 2431 | 122487 | A > G | |
| 23S | 2611 | 122667 | C > G | M129 has a T here |
| 23S | 2617 | 122673 | C > G |  |

#### Tetracycline SNPs:
CAMPNeu will check for the exact base change and will not report a variant if it doesn't confirm to the expectation listed in the table below.
| Gene Target | Gene Location | Genome Location | Variant | Note |
| ----------- | ------------- | ------------ | ------- | ----- |
| 16S | 968 | 119280 | T > C | M129 has a G at this location | 
| 16S | 1193 | 119505 | G > A | M129 has a T at this location |

#### Fluoroquinolone SNPs:
CAMPNeu uses snpEff to check for the specific amino acid changes listed below.
| Gene Target | Gene Location (nt) | Gene Location (aa) | Genome Location | Variant (nt) | Variant (aa)
| ----- | ----- | -------- | ------- | ----- | ------ |
| GyrA | 295 | 99 | 5115 | G > A | Asp > any |
| GyrB | 1327 | 443 | 4195 | G > A | Asp > any |
| GyrB | 1391 | 464 | 4259 | G > A | Arg > Lys |
| GyrB | 1448 | 483 | 4316 | A > G | Glu > Gly |
| ParC | 241 | 81 | 158614 | G > T | Gly > Cys |
| ParC | 248 | 83 | 158621 | C > T | Aps > any |
| ParC | 259 | 87 | 158632 | G > A | Asp > any |
| ParE | 1345 | 449 | 157811 | C > T | Pro > Ser |

### Required inputs: 
1. Illumina paired-end sequences

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
