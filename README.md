# scRNA_QC
Quality control of single cell RNA sequencing data


A pipeline for mapping and quality assessment single cell RNA-seq data

### See usage and parameters by typing:
python pipeline_master.py -h 

### What you need to run the program:
Directory of single cell sequencing data (parameter -i)
Mapping tool (parameter -m) - at least one of: bowtie, bowtie2, STAR, BWA, gsnap, salmon
Reference genome - either an index that fits the mapping tool (parameter -g) or as fasta and gtf (use -f and -gtf)
Quantification tool - at least one of: Tophat (run with bowtie or bowtie2), htseq-count or cufflinks 
Ulitity tools: samtools, gtf_splicesites, iit_store, picard-tools, bam2fastq, iit_store
R (with access to the internet for library installation and access to biomarRt)

<!--- list is to be continued - provide links to 3rd party tools -->



# Running the program
From the provided data path (input dir relative to provided root - see config file) the pipeline expects a directory named "raw" containing all the files in the experiment.
Files should be named with a hash and a number, that identifies the cell. e.g. for cell number 45. 
SingleCell/#45_1.fq
The two latter "_1" and ".fq" are identifiers for paired read number and fastq-file extention, respectively; and can be set in the config file.

For the most simple running of the pipeline create a directory (eg. /data1/scRNAdata/MouseKO/) for running of the pipeline. Store fastq files in a directory called "raw" (/data1/scRNAdata/MouseKO/raw). Store a fasta file and a gtf files of the reference gennome (e.g. /data1/reference_genomes/mm9.fasta, /data1/reference_genomes/mm9.gtf). Set appropriate paths in the config file (see config_example.txt)


## The config file explained
Please refer to the config_example.txt in order to understand the structure. 

[DIRECTORIES]
--- | ---
ROOT_DIR={path}                    | *prefix for input directory. Could be a single cell project directory. Following fastq files must be in {ROOT_DIR}/{INPUT DIRECTORY}/raw/\*
TEMP_DIR={path}                    | *all temporary files will be stored here - could be a fast non-backed up drive *
REF_DIR={path}                     | *path to reference genomes should contain directories named according to /{species_name}/{mapper_name} (-g and -m parameters). this can be build uding the. The reference genomes can be build by the pipeline, just requires setting -f parameter*
TOOLS_MAPPING =                    | *path to mapping tools - will be added to path and must thus contain executables in child directories*
TOOLS_QUANTIFICATION = {path}      | *prefix for quantification tools.*
MAPPING_ROOT = mapped              | *prefix for mapping tools*
QUANTIFICATION_ROOT = counts       | *directory name for count data*
PREPROCESSED_ROOT = preprocessed   | *directory name for preprocessed files*
SORTED_BAM_ROOT = sorted_bam       | *directory name for sorted bam files.*
SAM_ROOT = sam                     | *directory name for samfiles*


[EXTENSIONS]
EXT_BAM = .bam
EXT_FASTQ= .fq
EXT_FASTQ_GZ = .fq.gz
EXT_SAM = .sam
EXT_SORTED = .sorted
EXT_COUNTS = .counts
EXT_METRICS = .metrics
EXT_MERGED = .merged
EXT_MARKED_DUPL= .marked.duplicates
EXT_FASTA = .fa
EXT_LOG = .log
EXT_GTF= .gtf
EXT_GSNAP_SPLICE = .splice_sites
EXT_GSNAP_IIT = .iit
EXT_DICT = .dict
EXT_SUMMARY = .stat
FORWARD_STRAND = _1                | *Identifier used if analysing paired end reads.*
REVERSE_STRAND = _2                 

[SOFTWARE]
GSNAP=/bin/gsnap                  | *give path relative to MAPPING_ROOT*
GSNAP_BUILD=/bin/gmap_build       | *-*
BOWTIE1=/bowtie                   | *-*
BOWTIE1_BUILD=/bowtie-build       | *-*
BOWTIE2=/bowtie2                  | *-*
BOWTIE2_BUILD=/bowtie2-build      | *-*
BWA=/bwa                          | *-*
BWA_BUILD=/bwa                    | *-*
STAR = /source/STAR               | *-*
STAR_BUILD=/source/STAR           | *-*
SALMON = /bin/salmon              | *-*
SALMON_BUILD = /bin/salmon        | *-*
HTSEQTOOL=/scripts/htseq-count                                     | *give path relative to QUANTIFICATION_ROOT*
TOPHAT=/tophat2                                                    | *-*
CUFFLINKS=/cufflinks                                               | *-*
SAMTOOLS=/homes/ti1/tools/samtools-0.1.19/samtools                 | *-*
GTF_SPLICE_TOOL = /bin/gtf_splicesites                             | *-*
GTF_IIT_TOOL = /bin/iit_store                                      | *-*
PICARD_TOOL = /nfs/research2/teichmann/tools/picard-tools-1.113    | *give full path*
BAM_2_FASTQ=/homes/ti1/tools/bam2fastq-1.1.0/bam2fastq             | *-*





#Running the celluline pipeline on the pre-installed Amazon instance (AMI). 

In order to use the amazon instance you need to setup an AWS account [AMAZON AWS](https://aws.amazon.com). This requires a credit card number and a phone number for verification. There is a free trial for all new registrants, to try out the service, and institutions with teaching responsibilities might be eligible for further free credits. The AMI can also be exported as a VMware image for use on local clusters. [See Amazons guide](http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ExportingEC2Instances.html)

The AWS instance handles jobs via [GNU parallel](http://www.gnu.org/software/parallel/), and this job handling (pipeline parameter --cluster aws) is also well suited for dedicated servers, instances or farms which are not dependent on LSF queueing systems using bsub.

After creating a profile go to the console and choose ES2 and "Launch Instance" and then Community AMIs, where you search for celloline.

##Running the pipeline with test data.
We advice that you always try the pipeline with a small toy set for testing that everything works. This could be a subset of the cells and only the first few thousand reads. This can ususlly be run on the cheap t2.large instance, and will make sure that everything runs before upgrading the instance.

To get started with a test run quickly, start by making an index-file for the mapper of your choice:
```
pipeline.py -m gmap -f /home/ubuntu/references/Mouse/chr9.fa -gtf /home/ubuntu/references/Mouse/chr9.gtf -g chr9_fa -c /home/ubuntu/projects/example/config.txt -q htseq-count -i example --cluster aws --ram 10 -cpu 2
```
In the example we used only a small file with chromosome9. The full mouse genome is also preloaded in the same directory, but will take longer to run. Here we are providing fasta of the genome (-f) and gtf-file for the same genome (-gtf). Therefor the pipline will run only for creating a genome index - this only has to be done once for every mapper and genome. Hereafter the index-files can be retrieved by the pipeline.

Thereafter run the pipeline 
```
pipeline.py -m gmap -g chr9_fa -c /home/ubuntu/projects/example/config.txt -q htseq-count -i example --cluster aws --ram 10 -cpu 2
```

This is the same command without the genome fasta and gtf-file.

The pipeline will use the mapper *gmap* with the geome called *chr9_fa* with the configuration file */home/ubuntu/projects/example/config.tx*, quantification tool *htseq-count* using the *example* project and the GNU parallel job-handler (*aws*) requiring 10 mb of RAM free before starting new jobs on the server and allows for using 2 cpus for each job. 

## The AMI is structured as follows:
The configuration-file that you provide on runing the pipeline determine where the pipline will look for and also store files. It will be neccesary to attach more storage to the instance, when using full size files. For this purpose the setting of these directories in the top if th config is handy. Alternatively symbolic links can be used to achieve the same. In order to learn how to attach storage to the instance [EBS](http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ebs-attaching-volume.html)

###projects###
In the projects-directory is found and example project, and this is where new projects can be created. Each project directory should, when starting the pipeline, contain a directory called */raw* where uncompressed fastq-files are stored. These files should end with a *’#’* followed by the number of the cell, and an *”_”* followed by 1/2 indicating forward/reverse strand for paired end reads. Thus, in a paired read experiment cell number 84 would have two files that could be called
`singleCell#84_1.fastq and singleCell#84_2.fastq`

This can be comprised of a collection of symbolic links. Collecting and renaming of files depends on the naiming convention used by your sequencing facility, but some inspiration on some lines of code that could brung your about this task with numbering the cells quickly can be found in the (xxx), which must be altered to suit your filenames. We suggest not to rename the original files, in order to avoid loosing information or identity of the cells.
The project-directory also contains your **configuration file**. Please copy this file from the example. Please inspect this file, and set all the directories in the top of the file. 

###tools
Here are symbolic links the all the tools installed in the pipeline, to avoid a croveded bin.

###references
Here, the reference genomes are stored. On the AMI fasta and gtf-files for the mouse genome is stored along with smaller files only containing chomosome 9,for testing. Here the index-files for the mappers will also be stored.

###temp
This is the directory where temporary files are stored. This presently includes count-matices and statistics. This directory should be able to contain mapping and quantification output from all cells, and is better placed on a fast storage device with larger capacity (e.g. Amazon EBS), when running the pipeline on full size files. Many files in the temporary directory maybe deleted after run, for storage space considerations.

# Using your own data
We suggest to run a few files though the pipeline, create index files etc. using the smallest possible instance type ()
When using your own data in the pipeline you will need more storage and will have to attach an [EBS](http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/ebs-attaching-volume.html). We suggest to use the smallest possible instance (usuly) for making index and copying files and then later upgrade the instance to a larger one.



