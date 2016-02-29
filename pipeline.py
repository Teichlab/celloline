#!/usr/local/bin/ipython

__author__ = 'Tomislav Ilicic'
__copyright__ = 'Copyright (C) 2015 Ambrose Carr'
__license__ = 'MIT'
__version__ = 0.0

import os
import sys
import glob
import re
import ConfigParser
import subprocess
import shutil
import logging

import os, fnmatch
import numpy
import pandas as pd
from operator import itemgetter
from itertools import groupby
from os.path import basename
#LIST ALL DIRECTORIES IN SELECTED DIRECTORY. NOT RECURSIVE
def list_dirs(dir):
    walk = os.walk(dir)
    
    if (len(list(os.walk(dir))) == 0):
        return []
    else : 
        return next(os.walk(dir))[1]

def print_error (message):
    print "[ERROR]: " + message;
    sys.exit(1);


#CHECK IF OBJECT HAS BEEN SPECIFIED
def check_object_specified (type, object, objects_available):

    if (len(objects_available) == 0):
        print_error("No " + type + " are available. Please create one.");
    else:
        #IF OBJECT IS NOT AVAILABLE
        if (object not in objects_available):
            counter = 0;
            objects_available_num = [];
            for key in objects_available:
                counter =  counter + 1
                objects_available_num.append(str(counter) + ". " + str(key))
            #AND SUGGEST TO USER
            object = str(object)
            print_error(type + " \"" + object + "\" not available or not specified. Create or choose from:\n" + "\n".join(objects_available_num));

print
#TOP LEVEL FUNCTION EXECUTING THE PIPELINE
def run(args):

    logging.basicConfig(level=logging.INFO)

    if (not len(sys.argv) > 1):
        return False
    logging.info("Staring pipeline")
    code = 0    
    #USER WANTS TO CREATE REFERENCE GENOME INDEX FILE
    if (args.fasta != None and  args.gtf != None):
        logging.info("Genome index creation has started")
        code = run_ref_genome_index(args.mapping, args.fasta, args.gtf,args.cluster)
        if(code == 1):
                logging.error("Genome index creation unsuccessful")
                logging.error("Pipeline broke.")
                return(1)
        logging.info("Genome index creation complete")
    #ALL OTHER THINGS
    else:

        #INPUT NAME/DIR TO BE SET AT THIS POINT
        if (args.input != None):
            logging.info("Collecting input files")
            files = collect_input(args.input)
        else:
            print_error("Input folder (-i) not specified")

        user_out=ConfigSectionMap("DIRECTORIES")['temp_dir']+args.input+"/"+args.output
        #PIPELINE MAIN OUTPUT DIRECTORIES
        preprocessing_root       = user_out+"/" + ConfigSectionMap("DIRECTORIES")['preprocessed_root']
        mapping_root        = user_out+"/" + ConfigSectionMap("DIRECTORIES")['mapping_root']
        counts_root          = user_out+"/" + ConfigSectionMap("DIRECTORIES")['quantification_root']
        stats_root  = user_out+"/" + ConfigSectionMap("DIRECTORIES")['stats_root']
        #MAPPING CHOOSED
        if (args.mapping != None):
            logging.info("Mapping started")
            code = run_mapping(args.mapping, args.mapping_args, mapping_root, stats_root, args.input, files, args.genome, args.cluster,args.overwrite)
            if(code == 1):
                logging.error("Mapping unsuccessful")
                logging.error("Pipeline broke.")
                return(1)
            else:
                logging.info("Mapping sucessful")
        
        #QUANTIFICATION CHOOSED
        if (args.quantification != None):
            logging.info("Quantification started")
            code = run_quantification(args.quantification, args.quant_args, counts_root, mapping_root, args.input, files[1], args.genome, args.cluster, args.overwrite)
            if (code == 1): 
                logging.error("Quantification unsuccessful")
                logging.error("Pipeline broke.")
                return(1)
            else:
                logging.info("Quantification sucessful") 
        if (args.range  != None):
            range = args.range.split(",")
            range = map(int, range)
            files_process_index = set(range).intersection(files_process_index)
    logging.info("Pipeline is done.")
    return code

def write_stats(used_genome, sample_name,files_process_index, stats_root, mapping_root, cluster, overwrite):

    GTF = ConfigSectionMap("DIRECTORIES")['ref_dir'] + used_genome + "/" +  used_genome + ConfigSectionMap("EXTENSIONS")['ext_gtf']
    log_file = stats_root + "/stats.%I.log"
    sorted_sam_files = mapping_root + "/sam/" + sample_name + "_{#}" + ConfigSectionMap("EXTENSIONS")['ext_sam']
    _write_mapping_stats(GTF, files_process_index, sorted_sam_files, stats_root, sample_name, log_file, cluster, overwrite)

def _write_mapping_stats(GTF_file, files_process_index, sorted_sam_files,  stats_root, sample_name, log_file, cluster, overwrite):


    stats_output_merged = stats_root + "/" + sample_name + ".stats"
    logging.info("Generating mapping statistics.")
    expected_output_stats = list_expected_output(files_process_index, ".stats", stats_root, sample_name)
    if(overwrite == False and pre_checks(expected_output_stats) == True and os.path.exists(stats_output_merged)):
        return (0)
    stats_output = stats_root + "/" + sample_name + "_{#}.stats"

    commands = [] 
    commands.append("python /homes/ti1/code/git_code/st-method-comparison/generate_mapping_stats.py") 
    commands.append("-i")
    commands.append(sorted_sam_files)
    commands.append("-o")
    commands.append(stats_output)
    commands.append("-n")
    commands.append(sample_name + "_{#}.counts")
    commands.append("-g")
    commands.append(GTF_file)

    code = execute(files_process_index, commands, log_file, cluster)

    if (code !=0):
        return code
    commands = []
    commands.append("echo")
    
    s=""
    for i in range(0,10): 
        s += "S_"+str(i) + ","
    s += "S_10+"

    commands.append("cell,total,mapped,unmapped,unique,multi,intergenic,intragenic,exonic,intronic,ambigious,exonicM,alignments,multi-intergenic,multi-intragenic,multi-exonic,multi-intronic,multi-ambigious,perfect,partly_perfect,mapped_no_correct,"+s+",I,D,INDEL")
    commands.append(">")
    commands.append(stats_output_merged)
    proc = subprocess.Popen(" ".join(commands), stdout=subprocess.PIPE, shell=True)
    output = proc.stdout.read()
    
    commands = []
    commands.append("cat")
    commands.append(stats_root + "/" + sample_name + "_*.stats")
    commands.append(">>")
    commands.append(stats_output_merged)

    proc = subprocess.Popen(" ".join(commands), stdout=subprocess.PIPE, shell=True)
    output = proc.stdout.read()

    return (code)

    
#TOP LEVEL FUNCTION TO RUN QUANTIFICATION
def run_quantification(quantification, quant_args, counts_root, mapping_root, sample_name, files_process_index, genome, cluster, overwrite):
    available_quants = list_dirs(ConfigSectionMap("DIRECTORIES")['tools_quantification']);
    available_genomes = list_dirs(ConfigSectionMap("DIRECTORIES")['ref_dir']);
    #THRO ERROR IF QUANTIFICATION TOOL IS NOT AVAILABLE
    check_object_specified("Quantification tool (-q)", args.quantification, available_quants)

    #THROW ERROR IF GENOME NOT AVAILABLE
    if (args.genome == None):
        print_error("Reference genome (-g) option not specified")
    check_object_specified("Reference genome (-g)", genome, available_genomes);

    q_args = ""
    if (quant_args != None):
        q_args = quant_args[0].split(" ")
        del quant_args[0]
        q_args = q_args + quant_args


    return(_run_quantification(quantification, counts_root, mapping_root, sample_name, files_process_index, cluster, overwrite, genome, q_args))

#RUNS ALL QUANTIFICATION RELATED STUFF
def _run_quantification(quantifier, counts_dir, mapping_root, sample_name, files_process_index, cluster, overwrite, genome, quant_args):

    logging.info("Quantification tool:" + quantifier)
    logging.info("Passed parameters:" + quant_args)
    #Count dir example: /counts/htseq
    #Check if counts already exist
    quant_output_dir = counts_dir + "/" + quantifier.replace(".", "_") 
    expected_output_counts = list_expected_output(files_process_index, ConfigSectionMap("EXTENSIONS")['ext_counts'], quant_output_dir, sample_name)
    output_file = quant_output_dir + "/" + sample_name + "_" + quantifier.replace(".", "_")  + ConfigSectionMap("EXTENSIONS")['ext_counts']
    #If doesn't exist, create dir and run
    if(pre_checks(expected_output_counts) == False or overwrite == True or os.path.exists(output_file)==False):
        make_dir(quant_output_dir)
        
        #Check if required sorted bam files exist
        sorted_bam_dir = mapping_root + "/" + ConfigSectionMap("DIRECTORIES")['sorted_bam_root']
        sorted_bam_files = sorted_bam_dir + "/" + sample_name + "_{#}" + ConfigSectionMap("EXTENSIONS")['ext_sorted']
        #HACK
        #sorted_bam_files = mapping_root + "/sam/" + sample_name + "_{#}" + ConfigSectionMap("EXTENSIONS")['ext_sam']
        expected_input_sorted_bam = list_expected_output(files_process_index, ConfigSectionMap("EXTENSIONS")['ext_sorted'], sorted_bam_dir, sample_name)
        log_file = quant_output_dir + "/" + sample_name + "_%I" + ConfigSectionMap("EXTENSIONS")['ext_log'] 
        
        output_counts = quant_output_dir + "/" + sample_name + "_{#}" + ConfigSectionMap("EXTENSIONS")['ext_counts']
        #If sorted bam files do not exist, return 1
        if(not pre_checks(expected_input_sorted_bam, verbose=True)):
            print("Sorted bam files needed for quantification, but do not exist under dir:" + sorted_bam_dir)
            return (1)
        
        #Run selected quantification tool
        if("HTSeq" in quantifier):
            code = execute_htseq(sample_name,quantifier, quant_args, genome, sorted_bam_files, output_counts, files_process_index, log_file, cluster)
            merge(expected_output_counts, output_file,[1], -1)
        elif("cufflinks" in quantifier):
            code = 0
            code = execute_cufflink(sample_name, quantifier, quant_args, genome, sorted_bam_files, quant_output_dir, output_counts, files_process_index, log_file, cluster)
            merge(expected_output_counts, output_file,[9], 1)
        #Return code
        if(code != 0): 
            logging.error("Please check log files for more information: "+ quant_output_dir)
        if(code != 0 or not post_checks(expected_output_counts)):
            return 1    
        else:
            output_file = quant_output_dir + "/" + sample_name + "_" + quantifier.replace(".", "_")  + ConfigSectionMap("EXTENSIONS")['ext_counts']
    
    #Return code
    return 0

def merge(files, output_file, columns_list, skip_row_i):

    #COLUMNS TO EXTRACT FROM FILE
    columns = columns_list
    #SKIP ROW
    if (skip_row_i != -1):
        data = pd.read_table(files[0], skiprows=skip_row_i)
    else:
        data = pd.read_table(files[0], header=None)
    d = data.iloc[:, columns]
    row_names = data.iloc[:,0]
    merged = pd.DataFrame(d)
    file_names = []
    file_name =  os.path.basename(files[0])
    file_names.append(file_name)
    for i in  range(1, len(files)):
        f = files[i]
        file_name =  os.path.basename(f)
        file_names.append(file_name)


        if (skip_row_i != -1):
            data = pd.read_table(f, skiprows=skip_row_i)
        else:
            data = pd.read_table(f, header=None)

        d = data.iloc[:, columns]
        merged = pd.concat([merged,d], axis=1)

    merged.columns = file_names
    merged.index = row_names
    merged.index.name = "Genes"
    merged.to_csv(output_file, sep="\t")
    return (merged)
    
        
#EXECUTES CUFFLINK: USING SORTED BAM FILES BY NAME
def execute_cufflink(sample_name, quantifier, quant_args, genome_name, sorted_bam_files, quant_output_dir, output_counts, files_process_index, log_file, cluster):
    CUFFLINKS = ConfigSectionMap("DIRECTORIES")['tools_quantification'] + "/"+ quantifier + "/" + ConfigSectionMap("SOFTWARE")['cufflinks']
    gtf_file = ConfigSectionMap("DIRECTORIES")['ref_dir'] + "/" + genome_name + "/" + genome_name + ConfigSectionMap("EXTENSIONS")['ext_gtf']
    commands = []
    commands.append(CUFFLINKS)
    commands.append("-p")
    commands.append(args.cpu)
    commands.append("--output-dir")
    commands.append( quant_output_dir + "/{#}")
    commands.append("-G")
    commands.append(gtf_file)
    commands.append(sorted_bam_files)
    code = execute(files_process_index, commands, log_file, cluster)
     #Tophat produces sam files where you cannot change the name. So one needs to manually change name and location
    if(code == 0):
       move_files(files_process_index, quant_output_dir + "/{#}/genes.fpkm_tracking", quant_output_dir+"/"+sample_name, ConfigSectionMap("EXTENSIONS")['ext_counts'])
    return 0


#EXECUTES HTSEQI: USING SORTED BAM FILES BY NAME
def execute_htseq(sample_name, quantifier, quant_args, genome_name, sorted_bam_files, output_counts, files_process_index, log_file, cluster):
    HTSEQ = ConfigSectionMap("DIRECTORIES")['tools_quantification'] + "/"+ quantifier + "/" + ConfigSectionMap("SOFTWARE")['htseqtool']
    gtf_file = ConfigSectionMap("DIRECTORIES")['ref_dir'] + "/" + genome_name + "/" + genome_name + ConfigSectionMap("EXTENSIONS")['ext_gtf']
    commands = []
    commands.append(HTSEQ)
    commands.append("-f bam")
    commands.append("-r pos")
    commands.append("-a 0")
    commands.append("-s no")
    commands.append(" ".join(convert_to_default_map_param(quant_args)))
    commands.append(sorted_bam_files)
    commands.append(gtf_file) 
    commands.append(">")
    commands.append(output_counts)
    return(execute(files_process_index, commands, log_file, cluster))


#TOP LEVEL FUNCTION MAPPING ALGORITHM
def run_mapping(mapper,mapping_args, mapping_root, stats_root, sample_name, files, genome, cluster, overwrite):
    files_process_index  = files[1]
    files_index_holder = files[2]
    paired = files[3]
    read_length = files[4]

    available_mapper = list_dirs(ConfigSectionMap("DIRECTORIES")['tools_mapping']);
    available_genomes = list_dirs(ConfigSectionMap("DIRECTORIES")['ref_dir']);
    mapper_for_ref = mapper

    #THE PROBLEM WITH TOPHAT IS THAT IT RELIES ON BOWTIE
    #HENCE IT CAN'T BE TREATED INDIVIDUALLY.
    #THEREFORE IF TOPHAT WAS SELECTED, IT LOOKS UP WHICH
    #BOWTIE VERSION IS THE NEWEST AND GETS ALL AVAILABLE
    #REFERENCE GENOME FOR BOWTIE
    #if "tophat" not in mapper:
    if("tophat" in mapper):
        mapper_for_ref= get_latest_mapper("bowtie2")
        logging.info("Tophat will use Bowtie version: "+ mapper_for_ref)

    check_object_specified("Mapping tool (-m)", mapper_for_ref, available_mapper)

    #ONLY LIST GENOMES THAT ARE AVAIALBLE FOR A SPECIFIC MAPPING TOOL
    mapper_available_genomes = []
    for key in available_genomes:
        g = list_dirs(ConfigSectionMap("DIRECTORIES")['ref_dir'] + "/" + key);
        if (mapper_for_ref in g):
            mapper_available_genomes.append(key)

    #THROW ERROR IF GENOME NOT AVAILABLE
    if (args.genome == None):
        print_error("Reference genome (-g) option not specified")
    check_object_specified("Reference genome (-g)", genome, mapper_available_genomes);
    map_args = ""
    if (mapping_args != None):
        map_args = mapping_args[0].split(" ")
        del mapping_args[0]
        map_args = map_args + mapping_args

    files_index_holder = sorted(files_index_holder)
    #EXECUTE MAPPING
    code = (_run_mapping_all_steps(mapper,map_args, genome, files_index_holder, files_process_index, paired, read_length, mapping_root, sample_name, cluster, overwrite))


    if (code == 0):
        make_dir(stats_root)
        write_stats(genome, sample_name, files_process_index, stats_root, mapping_root, args.cluster, overwrite)
    
    return(code)
#TOP LEVEL FUNCTION TO CREATE REF GENOME INDEX FILES
def run_ref_genome_index(mapper, fasta_file, gtf, cluster):
    #CHECK IF GTF AND GENOME FASTA FILE HAVE BEEN SET

    available_mapper = list_dirs(ConfigSectionMap("DIRECTORIES")['tools_mapping']);
    if (not os.path.isfile(fasta_file)):
        print_error("FASTA sequencing file\" " + fasta_file + "\" does not exist");

    if (gtf == None):
        print_error("GTF file (-gtf) not specified.");

    if (not os.path.isfile(gtf)):
        print_error("GTF file \'" + gtf + "\' does not exist.");

    if (mapper == None):
        print_error("Mapping tool (-m) not specified.");
    check_object_specified("Mapping tool (-m)", mapper, available_mapper);
    code = _run_ref_genome_index_all_steps(mapper, fasta_file, gtf, cluster)
    return(code)
    
#RUNS THE MAPPING OPTION
# + SORTING AND CREATING BAM FILES
def _run_mapping_all_steps(mapper, mapper_args, used_genome, files_index_holder, files_process_index,paired, read_length, mapping_root, sample_name, cluster, overwrite):

    #MAPPINGOUTPUT DIRECTORIES
    mapping_dir         = mapping_root+"/"+ConfigSectionMap("DIRECTORIES")['sam_root'];
    sorted_dir          = mapping_root+"/"+ConfigSectionMap("DIRECTORIES")['sorted_bam_root'];

    logging.info("Mapper: " + mapper)
    logging.info("Passed parameters: " + " ".join(mapper_args))
    code = 0
    #1. MAPPING
    #CHECK IF FILES ALREADY EXIST
    expected_output_sam = list_expected_output(files_process_index, ConfigSectionMap("EXTENSIONS")['ext_sam'], mapping_dir, sample_name)
    if(pre_checks(expected_output_sam) == False or overwrite == True):
        make_dir(mapping_dir)
        used_genome_dir =  ConfigSectionMap("DIRECTORIES")['ref_dir'] + used_genome + "/" + mapper
        code = exec_mapping(mapper, mapper_args, read_length, used_genome, files_index_holder, files_process_index,paired,mapping_dir, sample_name, cluster)
        if(code != 0):
            logging.error("Please check log files for more information: "+ mapping_dir)
            return 1
        if(not post_checks(expected_output_sam)):
            return 1
    
    logging.info("Sorting (as bam)...")
    #2. CONVERT TO BAM AND SORT IT
    #CHECK IF FILES ALREADY EXIST
    expected_output_sorted = list_expected_output(files_process_index, ConfigSectionMap("EXTENSIONS")['ext_sorted'], sorted_dir, sample_name)
    if(pre_checks(expected_output_sorted) == False or overwrite == True):
        make_dir(sorted_dir)
        code = sort_sam_to_bam(files_process_index, mapping_dir, sorted_dir, sample_name, cluster)
        if(code != 0):
            logging.error("Please check log files for more information: "+ mapping_dir)
        if(not post_checks(expected_output_sorted)):
            return 1
    logging.info("Output-directory:" + mapping_root) 
    return 0

#SORT A SAM FILE AND CONVERT IT TO A BAM FILE    
def sort_sam_to_bam(files_process_index, mapper_dir, sorted_dir, sample_name, cluster):

    output_bam_file = sorted_dir + "/" + sample_name + "_{#}" + ConfigSectionMap("EXTENSIONS")['ext_sorted']
    input_sam_file = mapper_dir  + "/" + sample_name + "_{#}" + ConfigSectionMap("EXTENSIONS")['ext_sam']
    commands = []
    commands.append("java -jar")
    commands.append(ConfigSectionMap("SOFTWARE")['picard_tool'] + "/SortSam.jar")
    commands.append("INPUT="+input_sam_file);
    commands.append("OUTPUT="+ output_bam_file);
    #commands.append("SORT_ORDER=queryname");
    commands.append("SORT_ORDER=coordinate");
    commands.append("CREATE_INDEX=true");
    commands.append("VALIDATION_STRINGENCY=LENIENT");

    log_file = sorted_dir + "/sorting_%I"+ ConfigSectionMap("EXTENSIONS")['ext_log']
    code = execute(files_process_index, commands, log_file, cluster, True)
    return(code)


#CREATE A DIRECTORY RECURSIVELY
def make_dir(dir):
    if not os.path.exists(dir):
        os.makedirs(dir)

#RUNS CREATING INDEX FILE FOR A GENOME
def _run_ref_genome_index_all_steps(mapper, fasta_file, gtf, cluster):
        
    code = 0
    genome_name = basename(fasta_file)
    genome_name = genome_name.replace(".", "_")

    genome_root_dir = ConfigSectionMap("DIRECTORIES")['ref_dir'] + "/" + genome_name
    genome_mapper_dir = genome_root_dir + "/" + mapper
    exists = False

    #SINCE TOPHAT IS SPECIAL, DONT CREATE A DIR
    if ("tophat" in mapper):
        #CHECK TOPHAT TRANSCRITPOME INDEX IF IT EXISTS (under /bowtie2 index files)
        bowtie2_mapper = get_latest_mapper("bowtie2")
        bowtie2_transcritome_dir = (genome_root_dir + "/"+ bowtie2_mapper + "/transcriptome_data")
        genome_mapper_dir = bowtie2_transcritome_dir
    
    #CHECK IF GENOME MAPPING DIR ALREADY EXISTS
    if(os.path.isdir(genome_mapper_dir)):
        logging.warning("Genome index exists. Delete manually: "+ genome_mapper_dir)
        exists = True
    else:
        if("tophat" not in mapper):
            make_dir(genome_mapper_dir)
        else:
            make_dir(genome_root_dir)

    gtf_dest = genome_root_dir + "/" + genome_name + ConfigSectionMap("EXTENSIONS")['ext_gtf']
    if(not os.path.exists(gtf_dest)):
        #Copy GTF file to ref genome root dir
        shutil.copy2(gtf, gtf_dest)

    fasta_file_dest = genome_root_dir + "/" + genome_name + ConfigSectionMap("EXTENSIONS")['ext_fasta']
    if(not os.path.exists(fasta_file_dest)):
        #Copy fasta file to ref genome root dir
        shutil.copy2(fasta_file, fasta_file_dest)

    #TO DO
    #CREATE SEQUENCE DIRECTORY FILE
    seq_directory_file = fasta_file_dest + ConfigSectionMap("EXTENSIONS")['ext_dict']
    if(not os.path.exists(seq_directory_file)):
        logging.info("Creating sequence directory")
        code = create_seq_dir(fasta_file_dest, cluster)
        if(code != 0 or not os.path.exists(seq_directory_file)):
            logging.error("File not created:" + seq_directory_file)
        logging.info("Creating sequence directory was successful") 

    log_file = genome_root_dir + "/" + mapper + "_index.log"
    if (exists):
        return 0

    logging.info("Creating index files...(will take a while)")
    if ("bowtie-1" in mapper):
        code = create_bowtie1_index(mapper, fasta_file, gtf, genome_mapper_dir, genome_name,cluster, log_file)
    elif ("bowtie2" in mapper):
        code = create_bowtie2_index(mapper, fasta_file, gtf, genome_mapper_dir, genome_name,cluster, log_file)
    elif ("bwa" in mapper):
        code = create_bwa_index(mapper, fasta_file, gtf, genome_mapper_dir, genome_name,cluster, log_file)
    elif ("gmap" in mapper):
        code = create_gmap_index(mapper, fasta_file, gtf, genome_mapper_dir, genome_name,cluster, log_file)
    elif ("STAR" in mapper):
        code = create_star_index(mapper, fasta_file, gtf, genome_mapper_dir, genome_name,cluster, log_file)
    elif ("tophat" in mapper):
        code = create_tophat_index(mapper, fasta_file_dest, gtf_dest, genome_mapper_dir, genome_name, cluster, log_file)

    if (code == 0):
        logging.info("Creating index files was successful:" + genome_mapper_dir)
    else:
        logging.error("Creating index files failed. Check logs:" + log_file)

    return(code)

#CREATES SEQUENCE DIRECTORY FILE NEEDED FOR QUANTIFICATION
def create_seq_dir(fasta_file_dest, cluster):
    
    log_file = fasta_file_dest + ConfigSectionMap("EXTENSIONS")['ext_dict'] + ConfigSectionMap("EXTENSIONS")['ext_log']
    picard = ConfigSectionMap("SOFTWARE")['picard_tool']
    commands = []
    commands.append("java -jar")
    commands.append(picard+"/CreateSequenceDictionary.jar")
    commands.append("REFERENCE=" + fasta_file_dest)
    commands.append("OUTPUT=" + fasta_file_dest + ConfigSectionMap("EXTENSIONS")['ext_dict'])
    return(execute([], commands, log_file, cluster))
    

#LIST EXPECTED OUTPUT
def list_expected_output(files_process, ext, out_dir, sample_name):
    expected = []

    for i in files_process:
        #dir/sampleName_index_.ext
        expected.append(out_dir + "/" + sample_name + "_" + str(i) + ext)
    return expected    

#EXECUTE A GIVEN MAPPING ALGORITHM ON A GIVEN CLUSTER SYSTEM
def exec_mapping(mapper, mapper_args, read_length, used_genome, files_index_holder, files_process_index,paired,mapping_dir, sample_name, cluster):
    code = 0 
    if ("bowtie-1" in mapper):
        code = exec_bowtie1_mapping(mapper, mapper_args, read_length, used_genome, files_index_holder, files_process_index,paired,mapping_dir, sample_name, cluster)
    elif ("bowtie2" in mapper):
        code = exec_bowtie2_mapping(mapper, mapper_args, read_length, used_genome, files_index_holder, files_process_index,paired,mapping_dir, sample_name, cluster)
    elif ("bwa" in mapper):
        code = exec_bwa_mapping(mapper, mapper_args, read_length, used_genome, files_index_holder, files_process_index,paired,mapping_dir, sample_name, cluster)
    elif ("gmap" in mapper):
        code = exec_gmap_mapping(mapper, mapper_args, read_length, used_genome, files_index_holder, files_process_index,paired,mapping_dir, sample_name, cluster)
    elif ("STAR" in mapper):
        code = exec_star_mapping(mapper, mapper_args, read_length, used_genome, files_index_holder, files_process_index,paired,mapping_dir, sample_name, cluster)
    elif ("tophat" in mapper):
        code = exec_tophat_mapping(mapper, mapper_args, read_length, used_genome, files_index_holder, files_process_index,paired,mapping_dir, sample_name, cluster)

    return(code)

def exec_star_mapping(mapper, mapper_args, read_length, used_genome, files_index_holder, files_process_index,paired,mapping_dir, sample_name, cluster):
    star = ConfigSectionMap("DIRECTORIES")['tools_mapping'] + "/"+ mapper + ConfigSectionMap("SOFTWARE")['star']
    genome_dir =  ConfigSectionMap("DIRECTORIES")['ref_dir'] + used_genome + "/" + mapper
    code = 0
    commands = []
    GTF = ConfigSectionMap("DIRECTORIES")['ref_dir'] + used_genome + "/" +  used_genome + ConfigSectionMap("EXTENSIONS")['ext_gtf']
    log_file1 = mapping_dir + "/" + sample_name + "_%I.mapping.log"
    output_sam_file = mapping_dir + "/" + sample_name + "_{#}" + ConfigSectionMap("EXTENSIONS")['ext_sam']
    commands.append(star)
    commands.append("--outFileNamePrefix")
    commands.append(mapping_dir + "/" + sample_name + "_{#}.")
    commands.append("--genomeDir")
    commands.append(genome_dir)
    commands.append("--readFilesIn")
    commands.append(" ".join(files_index_holder))
    commands.append("--genomeLoad")
    commands.append("NoSharedMemory")
    commands.append("--outSAMunmapped")
    commands.append("Within")
    commands.append("--outStd")
    commands.append("SAM")
    #TEMP FOLDER WHERE EVERYTHING IS STORED
    #NEEDS TO BE DELETED OTHERWISE ERROR FROM STAR
    TEMP = mapping_dir  + "/" + sample_name + "_{#}" + "/"
    for i in files_process_index:
        f = TEMP.replace("{#}", str(i))
        if(os.path.isdir(f)):
            shutil.rmtree(f)
    #commands.append("--outTmpDir")
    #commands.append(TEMP)

    options = convert_to_default_map_param(mapper_args)
    splice = "-splice" in options

    if (splice):
       commands.append("--sjdbGTFfile")
       commands.append(GTF) 
       options.remove("-splice")

    commands.append(" ".join(options))
    commands.append(">")
    commands.append(output_sam_file)

    return(execute(files_process_index, commands, log_file1, cluster, True))

def exec_gmap_mapping(mapper, mapper_args, read_length, used_genome, files_index_holder, files_process_index,paired,mapping_dir, sample_name, cluster): 
    gmap = ConfigSectionMap("DIRECTORIES")['tools_mapping'] + "/"+ mapper + ConfigSectionMap("SOFTWARE")['gsnap']
    genome_dir =  ConfigSectionMap("DIRECTORIES")['ref_dir'] + used_genome + "/" + mapper
    code = 0
    commands = [] 
   
    fasta_size = os.stat(ConfigSectionMap("DIRECTORIES")['ref_dir'] + used_genome + "/" + used_genome + ConfigSectionMap("EXTENSIONS")['ext_fasta']).st_size
    fasta_size_GB = fasta_size / 1024.0 /1024.0 / 1024.0
    if (fasta_size_GB >= 5.0) :
        gmap = ConfigSectionMap("DIRECTORIES")['tools_mapping'] + "/"+ mapper + ConfigSectionMap("SOFTWARE")['gsnapl']
    output_sam_file = mapping_dir + "/" + sample_name + "_{#}" + ConfigSectionMap("EXTENSIONS")['ext_sam'] 
    commands.append(gmap)
    commands.append("-A sam -B 5")
    commands.append("-t")
    commands.append(args.cpu)
   
    options = convert_to_default_map_param(mapper_args)
    splice = "-splice" in options
    splice_f = genome_dir + "/" + used_genome + ConfigSectionMap("EXTENSIONS")['ext_gsnap_splice']
    print (splice_f)
    #Only add splice file if not empty. Otherwise gmap results in a weired error
    if(splice == True and os.stat(splice_f).st_size > 0):
        commands.append("-s")
        commands.append(genome_dir + "/" + used_genome)
        commands.append("-N")
        commands.append(1)
        options.remove("-splice") 

    commands.append(" ".join(options))
    commands.append("-d") 
    commands.append(used_genome)
    commands.append("-D")
    commands.append(genome_dir)

    #COMMAND
    log_file1 = mapping_dir + "/" + sample_name + "_%I.mapping.log"
    #Peform two alignments and sampe command and convert to .sam file
    if(paired):
        commands.append(files_index_holder[0])
        commands.append(files_index_holder[1])
    else:
        commands.append(files_index_holder[0]) 

    commands.append(">")
    commands.append(output_sam_file)
    return(execute(files_process_index, commands, log_file1, cluster, True))

def exec_bowtie1_mapping(mapper, mapper_args, read_length, used_genome, files_index_holder, files_process_index,paired,mapping_dir, sample_name, cluster): 
    BOWTIE1 = ConfigSectionMap("DIRECTORIES")['tools_mapping'] + "/"+ mapper + ConfigSectionMap("SOFTWARE")['bowtie1']
    genome_dir =  ConfigSectionMap("DIRECTORIES")['ref_dir'] + used_genome + "/" + mapper
    code = 0
    commands = [] 
   
    output_sam_file = mapping_dir + "/" + sample_name + "_{#}" + ConfigSectionMap("EXTENSIONS")['ext_sam'] 
    commands.append("export BOWTIE_INDEXES=" +genome_dir+ " && ")
    commands.append(BOWTIE1)
    commands.append("-S")
    commands.append(genome_dir + "/" + used_genome)
 
    #COMMAND
    log_file1 = mapping_dir + "/" + sample_name + "_%I.mapping.log"
    #Peform two alignments and sampe command and convert to .sam file
    if(paired):
        commands.append("-1") 
        commands.append(files_index_holder[0])
        commands.append("-2")
        commands.append(files_index_holder[1])
    else:
        commands.append(files_index_holder[0]) 
    commands.append(" ".join(convert_to_default_map_param(mapper_args)))
    commands.append(">")
    commands.append(output_sam_file)
    return(execute(files_process_index, commands, log_file1, cluster, True))

def exec_bowtie2_mapping(mapper, mapper_args, read_length, used_genome, files_index_holder, files_process_index,paired,mapping_dir, sample_name, cluster): 
    BOWTIE1 = ConfigSectionMap("DIRECTORIES")['tools_mapping'] + "/"+ mapper + ConfigSectionMap("SOFTWARE")['bowtie2']
    genome_dir =  ConfigSectionMap("DIRECTORIES")['ref_dir'] + used_genome + "/" + mapper
    code = 0
    commands = [] 
   
    output_sam_file = mapping_dir + "/" + sample_name + "_{#}" + ConfigSectionMap("EXTENSIONS")['ext_sam'] 
    commands.append("export BOWTIE2_INDEXES=" +genome_dir+ " && ")
    commands.append(BOWTIE1)
    commands.append("-x")
    commands.append(genome_dir + "/" + used_genome)
 
    #COMMAND
    log_file1 = mapping_dir + "/" + sample_name + "_%I.mapping.log"
    #Peform two alignments and sampe command and convert to .sam file
    if(paired):
        commands.append("-1") 
        commands.append(files_index_holder[0])
        commands.append("-2")
        commands.append(files_index_holder[1])
    else:
        commands.append(files_index_holder[0]) 
    commands.append(" ".join(convert_to_default_map_param(mapper_args)))
    commands.append(">")
    commands.append(output_sam_file)
    return(execute(files_process_index, commands, log_file1, cluster, True))

def pre_checks(expected_output, verbose=False):
    
    not_exists = []
    for  exp in expected_output:
        if(not os.path.exists(exp)):
            not_exists.append(exp)
    if (len(not_exists) == 0):
        if(not verbose):
            logging.info("Output already exists:" + ",".join(expected_output))
        return True
    return False

def post_checks(expected_output):
       
    not_exists = []
    for  exp in expected_output:
        if(not os.path.exists(exp)):
            not_exists.append(exp)

    if (len(not_exists) != 0):
        logging.error("Following files don't exists:" + ",".join(not_exists))
        return False
    return True

#CREATE BOWTIE2 INDEX FILES
def create_gmap_index(mapper, fasta_file, gtf, genome_dir, genome_name, cluster, log_file):

    build = ConfigSectionMap("DIRECTORIES")['tools_mapping'] + "/"+ mapper + ConfigSectionMap("SOFTWARE")['gsnap_build']
    gtf_iit_tool = ConfigSectionMap("DIRECTORIES")['tools_mapping'] + "/"+ mapper + ConfigSectionMap("SOFTWARE")['gtf_iit_tool']
    gtf_splice_tool = ConfigSectionMap("DIRECTORIES")['tools_mapping'] + "/"+ mapper + ConfigSectionMap("SOFTWARE")['gtf_splice_tool']

    gtf_out = genome_dir + "/" +genome_name +  ConfigSectionMap("EXTENSIONS")['ext_gsnap_splice']
    
    logging.info("Creating splice sites file...")
    commands = []
    commands.append("cat")
    commands.append(gtf)
    commands.append(" | " + gtf_splice_tool + " -")
    commands.append(">")
    commands.append(gtf_out)
    code = execute([], commands, log_file, cluster)
    if(code == 0):
        logging.info("Successful")
    else:
        logging.error("Unsuccessful")
        return code    
    
    splice_file = genome_dir + "/" + genome_name + ConfigSectionMap("EXTENSIONS")['ext_gsnap_iit']; 
    logging.info("Creating splice sites file II ...")
    commands = []
    commands.append("cat")
    commands.append(gtf_out)
    commands.append(" | " + gtf_iit_tool)
    commands.append("-o " + splice_file)
    code = execute([], commands, log_file, cluster)
    if(code == 0):
        logging.info("Successful")
    else:
        logging.error("Unsuccessful")
        return code 

    commands = []
    commands.append(build)
    commands.append("-d")
    commands.append(genome_name)
    commands.append(fasta_file)
    commands.append("-D")
    commands.append(genome_dir)
    code = execute([], commands, log_file, cluster, True)

    return(code)



#CREATE BOWTIE2 INDEX FILES
def create_bowtie2_index(mapper, fasta_file, gtf, genome_dir, genome_name, cluster, log_file):
    bowtie_build = ConfigSectionMap("DIRECTORIES")['tools_mapping'] + "/"+ mapper + ConfigSectionMap("SOFTWARE")['bowtie2_build']

    commands = []
    commands.append("export BOWTIE_INDEXES=" + genome_dir + " && (")
    commands.append(bowtie_build)
    commands.append(fasta_file)
    commands.append(genome_name)
    commands.append(";mv " + genome_name + "*.bt2 "+ genome_dir + ")")
    code = execute([], commands, log_file, cluster, True)

    return(code)

#CREATE BOWTIE1 INDEX FILES
def create_bowtie1_index(mapper, fasta_file, gtf, genome_dir, genome_name, cluster, log_file):
    bowtie_build = ConfigSectionMap("DIRECTORIES")['tools_mapping'] + "/"+ mapper + ConfigSectionMap("SOFTWARE")['bowtie1_build']

    commands = [] 
    commands.append("export BOWTIE_INDEXES=" + genome_dir + " && (")
    commands.append(bowtie_build)
    commands.append(fasta_file)
    commands.append(genome_name)
    commands.append(";mv " + genome_name + "*.ebwt "+ genome_dir + ")") 
    code = execute([], commands, log_file, cluster, True)

    return(code)

#CREATE BOWTIE1 INDEX FILES
def create_star_index(mapper, fasta_file, gtf, genome_dir, genome_name, cluster, log_file):
    build = ConfigSectionMap("DIRECTORIES")['tools_mapping'] + "/"+ mapper + ConfigSectionMap("SOFTWARE")['star_build']

    commands = []
    commands.append(build)
    commands.append("--runMode genomeGenerate")
    commands.append("--genomeDir")
    commands.append(genome_dir)
    commands.append("--genomeFastaFiles")
    commands.append(fasta_file)
    commands.append("--limitGenomeGenerateRAM 186718572928")
    commands.append("--genomeChrBinNbits 18")
    code = execute([], commands, log_file, cluster, True)

    return(code)


#CREATE BWA INDEX FILES
def create_bwa_index(mapper, fasta_file, gtf, genome_dir, genome_name, cluster, log_file):
    build = ConfigSectionMap("DIRECTORIES")['tools_mapping'] + "/"+ mapper + ConfigSectionMap("SOFTWARE")['bwa_build']
    commands = [] 
    commands.append("(")
    commands.append(build)
    commands.append("index")
    commands.append("-a bwtsw")
    commands.append("-p")
    commands.append(genome_name)
    commands.append(fasta_file)
    move = genome_name + "*.amb " +genome_name + "*.ann " +genome_name + "*.pac " +genome_name + "*.bwt " +genome_name + "*.sa "
    commands.append(";mv " + move + " "+ genome_dir + ")") 
    code = execute([], commands, log_file, cluster, True)

    return(code)


#CREATE TOPHAT INDEX FILES
def create_tophat_index(mapper, fasta_file, gtf, genome_dir, genome_name, cluster, log_file):
    

    bowtie2_mapper = get_latest_mapper("bowtie2")
    TOPHAT = ConfigSectionMap("DIRECTORIES")['tools_mapping'] + "/"+ mapper + ConfigSectionMap("SOFTWARE")['tophat']

    #IMPORTANT TO USE THE / AT THE END!!
    used_genome_dir =  ConfigSectionMap("DIRECTORIES")['ref_dir'] + genome_name + "/" + bowtie2_mapper + "/" 

    if(not os.path.isdir(used_genome_dir)):
        logging.error("Genome index: "+ genome_name + " not available for Tophat. Please select a bowtie2 mapper to generate genome index first, before creating transcriptome.")
        return(1)

    #TOPHAT NEEDS FASTA FILE NEXT TO OTHER BOWTIE2 INDEX FILES
    shutil.copy2(fasta_file, used_genome_dir)
    code = 0 
    commands = []
    commands.append("export BOWTIE2_INDEXES=" + used_genome_dir + " && ")
    commands.append("export PATH=\$PATH:"+ ConfigSectionMap("DIRECTORIES")['tools_mapping'] + "/"+ bowtie2_mapper + " && ")
    commands.append(TOPHAT) 
    commands.append("-G " + gtf) 
    commands.append("--transcriptome-index " + "transcriptome_data/known")
    commands.append(genome_name)
    
    #Important to change current working dir. As all files will be created there
    current_dir = os.path.dirname(os.path.realpath(__file__))
    os.chdir(used_genome_dir)
    code = execute([], commands, log_file, cluster, True)
    os.chdir(current_dir)

    return(code)

#MAP USING TOPHAT2 
#MAPS FIRST TO GENOME AND THEN TO TRANSCRIPTOME
def exec_tophat_mapping(mapper, mapper_args, read_length, used_genome, files_index_holder, files_process_index,paired,mapping_dir, sample_name, cluster):
    
    bowtie2_mapper = get_latest_mapper("bowtie2")
    log_file = mapping_dir + "/" + sample_name + "_%I.mapping.log"
    TOPHAT = ConfigSectionMap("DIRECTORIES")['tools_mapping'] + "/"+ mapper + ConfigSectionMap("SOFTWARE")['tophat']
   
    #IMPORTANT TO USE THE / AT THE END!!
    used_genome_dir =  ConfigSectionMap("DIRECTORIES")['ref_dir'] + used_genome + "/" + bowtie2_mapper + "/"
    code = 0
    commands = []
    commands.append("export BOWTIE2_INDEXES=" + used_genome_dir + " && ")
    commands.append("export PATH=\$PATH:"+ ConfigSectionMap("DIRECTORIES")['tools_mapping'] + "/"+ bowtie2_mapper + " && ")
    commands.append(TOPHAT) 
    commands.append("-o " + mapping_dir + "/{#}")
    #commands.append("--no-convert-bam")
    commands.append("-p")
    commands.append(args.cpu)

    options = convert_to_doublets_map_param(mapper_args)
    splice = "--splice" in options
    
    if (splice):
        options.remove("--splice")
        commands.append("--transcriptome-index " + used_genome_dir+"/transcriptome_data/known")
    commands.append(" ".join(options))

    commands.append(used_genome)
    commands.append(" ".join(files_index_holder))
    code = execute(files_process_index, commands, log_file, cluster, True)
    
    log_file = mapping_dir + "/" + sample_name + "_%I.merge_bam.log"
    #Tophat produces sam files where you cannot change the name. So one needs to manually change name and location
    if(code == 0):
        code = bam_to_sam(files_process_index, [mapping_dir + "/{#}/accepted_hits.bam", mapping_dir + "/{#}/unmapped.bam"], mapping_dir, sample_name, ConfigSectionMap("EXTENSIONS")['ext_sam'], cluster, log_file)
         
    return(code)

#CONVERT BAM TO FILES
def bam_to_sam(files_process_index, files, dest_dir, prefix, ext, cluster, log_file):
    input_bam_file = files
    output_sam_file = dest_dir  + "/" + prefix + "_{#}" +  ext
    commands = []
    commands.append("java -jar")
    commands.append(ConfigSectionMap("SOFTWARE")['picard_tool'] + "/MergeSamFiles.jar")
    for f in files:
        commands.append("INPUT="+f);
    commands.append("OUTPUT="+ output_sam_file);
    commands.append("SORT_ORDER=unsorted");
    commands.append("VALIDATION_STRINGENCY=LENIENT");
    code = execute(files_process_index, commands, log_file, cluster)
    return(code)
        

#CONCATINATE FILES
def concat_files(files_process_index, files, dest, ext):

    all_files = []
    files_process_index_t = list(files_process_index)
    for i in files_process_index_t:
        temp_f = []
        for f in range(0, len(files)):
            temp_f.append(files[f].replace("{#}", str(i)))
        all_files.append(temp_f)

    #CONCATINATE
    for i in range(0, len(files_process_index_t)):
        output_file = dest+"_"+str(files_process_index_t[i])+ext
        temp_f = all_files[i] 
        with open(output_file, 'w') as outfile:
            for fname in temp_f:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)
        outfile.close()

def move_files(files_process_index, source_dir, dest, ext):
   for i in files_process_index:
        s = source_dir.replace("{#}", str(i))
        d = dest+"_"+str(i)+ext
        shutil.move(s, d) 
   
#ONVERT PARAM TO DOUBLE WAY: e.g. --i X
def convert_to_doublets_map_param(mapper_args):
    m_args = [] 
    for arg in mapper_args:            
        m_args.append("--"+arg.replace("=", " "))
    return(m_args) 

#CONVERT PARAM TO DEFAULT WAY: e.g. -i X
def convert_to_default_map_param(mapper_args):
    m_args = []
    for arg in mapper_args:       
        m_args.append("-"+arg.replace("=", " "))
    return(m_args) 

#FOR A GIVEN ALGORITHM, GET THE LATEST VERSION AVAILABLE LOCALLY
def get_latest_mapper(mapper):
    available_mapper = list_dirs(ConfigSectionMap("DIRECTORIES")['tools_mapping']);
    available_mapper = sorted(available_mapper) 
    
    for i in range(0, len(available_mapper)):
        latest = available_mapper[i]
        if(mapper in latest):
            return(latest)

    return("")
#EXECUTE BWA SHORT OR LONG MAPPING
#DECIDES AUTOMATLICALLY IF SHORT OR LONG
def exec_bwa_mapping(mapper, mapper_args, read_length, used_genome, files_index_holder, files_process_index,paired,mapping_dir, sample_name, cluster):
    BWA = ConfigSectionMap("DIRECTORIES")['tools_mapping'] + "/"+ mapper + ConfigSectionMap("SOFTWARE")['bwa']
    used_genome_dir =  ConfigSectionMap("DIRECTORIES")['ref_dir'] + used_genome + "/" + mapper
    code = 0
    commands = []
   
    output_sam_file = mapping_dir + "/" + sample_name + "_{#}" + ConfigSectionMap("EXTENSIONS")['ext_sam'] 
    #IF READ SHORT, USE ALN COMMAND
    if (read_length  < 30):
        commands.append(BWA)
        commands.append("aln")
        commands.append("-t")
        commands.append(args.cpu)
        commands.append(" ".join(convert_to_default_map_param(mapper_args)))
        commands.append(used_genome_dir + "/" + used_genome)
            
        #COMMAND
        input_file_1 = files_index_holder[0]
        output_sai_file1 = mapping_dir + "/" + sample_name + "_{#}.sai"
        commands.append(input_file_1)
        commands.append(">")
        commands.append(output_sai_file1)
        log_file1 = mapping_dir + "/" + sample_name + "_%I.mapping.log"
        #Peform two alignments and sampe command and convert to .sam file
        if(paired):

            #Peform two alignments
            output_sai_file1 = mapping_dir + "/" + sample_name + "_{#}_1.sai"
            output_sai_file2 = mapping_dir + "/" + sample_name + "_{#}_2.sai"
    
            commands[5]=output_sai_file1
            log_file1 = mapping_dir + "/" + sample_name + "_%I_1.sai.log"
            commands[3] = files_index_holder[0]
            code1 = execute(files_process_index, commands, log_file1, cluster)     

            output_sai_file2 = mapping_dir + "/" + sample_name + "_{#}_2.sai"
            commands[5]=output_sai_file2
            log_file2 = mapping_dir + "/" + sample_name + "_%I_2.sai.log"
            commands[3] = files_index_holder[1]
            code2 = execute(files_process_index, commands, log_file2, cluster, True)
        
            #sample command to conver to .sam file
            commands = []
            commands.append(BWA)
            commands.append("sampe")
            commands.append(used_genome_dir + "/" + used_genome)
            commands.append(output_sai_file1)
            commands.append(output_sai_file2) 
            commands.append(files_index_holder[0])
            commands.append(files_index_holder[1])
            commands.append(">")
            commands.append(output_sam_file)
            log_file = mapping_dir + "/" + sample_name + "_%I.sampe.log"
            code3 = execute(files_process_index, commands, log_file, cluster, True)

            code = code1 + code2 +code3
        #Run one alignment (single-end) and convert to .sam file
        else:
            code1 = execute(files_process_index, commands, log_file1, cluster, True)
            commands = []
            commands.append(BWA)
            commands.append("samse")
            commands.append(used_genome_dir + "/" + used_genome)
            commands.append(output_sai_file1)
            commands.append(files_index_holder[0])
            commands.append(">")
            commands.append(output_sam_file)
            log_file = mapping_dir + "/" + sample_name + "_%I.sampe.log"
            code3 = execute(files_process_index, commands, log_file, cluster, True)
    
    #LONG READS
    else:
        commands = []
        commands.append(BWA)
        commands.append("mem")
        commands.append("-t")
        commands.append(args.cpu)
        commands.append(" ".join(convert_to_default_map_param(mapper_args)))
        commands.append(used_genome_dir + "/" + used_genome)
        commands.append(files_index_holder[0])
        log_file = mapping_dir + "/" + sample_name + "_%I.mapping.log"
        if(paired):
            commands.append(files_index_holder[1])
        commands.append(">")
        commands.append(output_sam_file)
        code = execute(files_process_index, commands, log_file, cluster, True)

    return code

#RECEIVE READ LENGTH OF FILE (2nd line)
def get_read_length(file):
    read_length = 0
    
    in_file = open(file, "r")
    line = in_file.readline()
    line = in_file.readline()
    #-1 because of new line
    read_length = len(line) - 1
    in_file.close()
    return(read_length) 

def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(os.path.join(root, name))
    return result


#FUNCTION THAT EXECUTES COMMAND ONTO CLUSTER
def execute(files_process_index, commands, log_file, cluster, speed=False):

    code = -1
    
    if (len(files_process_index) ==0 ):
        files_process_index.append("1")
    files_process_index =  sorted(files_process_index)
    files_process_index = list(files_process_index)
    
    ranges = []
    num_files = len(files_process_index)
    if (len(files_process_index) > 1):
        for k, g in groupby(enumerate(files_process_index), lambda (i,x):i-x):
            group = map(itemgetter(1), g)
            ranges.append(str(group[0]) + "-" + str(group[-1]))
    
        files_process_index = ranges
    files_process_index = [str(i) for i in files_process_index]
    
    #WRITE COMMAND TO LOG FILE
    
    logging.info("Executed Command: " + " ".join(commands)) # python will convert \n to os.linesep
    #CHOSE EITHER OF THE TWO CLUSTERING SYSYEMS
    if (cluster == "ebi"):
        code = run_on_EBI(num_files, files_process_index, commands, log_file, speed)
    elif (cluster == "aws"):
        code = run_on_AWS(files_process_index, commands, log_file)
    return code

#CONVERT PARAMETERS TO ONE COMMAND THAT CAN BE RUN ONTO THE EBI CLUSTER
def run_on_EBI(num_files, files_process_index, commands, log_file, speed=False):
        import time
        code = 0

        #BUILD COMMAND
        cluster_command = []
        edit_i = []
        for i in range(0,len(commands)):
            if('{#}' in commands[i]):
                edit_i.append(i)

        #REPLACE PLACEHOLDER WITH LSB JOBINDEX
        for i in edit_i:        
            commands[i] = commands[i].replace("{#}", "\\${LSB_JOBINDEX}")
#        select_statement = "-R \"select[panfs_nobackup_research]\""
        cluster_command.append("bsub")
    
        ram = "50000"
        cpu = "2"
        if (speed == True):
            if (args.ram != None):
                ram = args.ram
            if (args.cpu != None):
                cpu = args.cpu

        select_statement = "-R \"select[panfs_nobackup_research] rusage[mem=" + ram + "]\""
        cluster_command.append("-M " + ram)
        cluster_command.append("-n " + cpu)

        cluster_command.append(select_statement)
        cluster_command.append("-J \"["+ ",".join(files_process_index) + "]\"")
        cluster_command.append("-oo " + log_file)
        cluster_command.append("\"")
        cluster_command.append(" ".join(commands))
        cluster_command.append("\"")
        #EXECUTE PROCESS
        proc = subprocess.Popen(" ".join(cluster_command), stdout=subprocess.PIPE, shell=True)
        output = proc.stdout.read()
        #Feth jobid        
        job_id_sub=output.split()[1]
        job_id_sub=job_id_sub[1:len(job_id_sub)-1]
        logging.info("Submitted job ID: " + job_id_sub)

         #CHECK STATUS
        not_done = True
        startProgress("Process has started");
        while (not_done):
            cluster_command = []
            cluster_command.append("bjobs "+ job_id_sub)
            proc = subprocess.Popen(" ".join(cluster_command), stdout=subprocess.PIPE, shell=True)
            output = proc.stdout.read()

            status = output.split("\n")
            #IF NOTHING POPS UP, EVERYTHING IS DONE
            if (len(status) <=1):
                not_done = False
                endProgress()
            #PARSE STATUS OF JOBS 
            else:
                status = status[1:len(status)-1]
                running = []
                completed = []
                exit = []
                uknown = []
                pending = [] 
                #ADD EACH TYPE OF STATUS TO A SEPERATE ARRAY
                for i in range(0,len(status)):
                    stat =  status[i]
                    des = stat.split()
                    task_id = des[len(des)-4]
                    task_id = task_id[1:len(task_id)-1]
                    current_status = des[2]
                    if (current_status == "RUN"):
                        running.append(task_id)
                    elif (current_status == "EXIT"):
                        exit.append(task_id)
                    elif (current_status == "PEND"):
                        pending.append(task_id)
                    elif (current_status == "DONE"):
                        completed.append(task_id)
                    else: 
                        uknown.append(task_id)
                
                #CHECK IF COMPLETED
                #IF ERRORS RE-RUN ONE MORE TIME. THEN SHUT DOWN 
                #TO DO
                if (len(completed) == (num_files) or (len(running) == 0 and len(pending) == 0 and (len(exit) + len(uknown) == (num_files)))): 
                    not_done = False
                    code = 0
                    endProgress()
                    logging.info("Jobs Completed:" + str(len(completed)) + "/" + str((num_files)))
                    if(len(exit) > 0):
                        logging.error("Jobs Exited:" + str(len(exit)) + "/" + str((num_files)))
                        code = 1
                #PRINT STATUS
                else :
                    not_done = True
                    p = float(len(completed)+len(exit))/float((num_files))
                    progress(p*100)
                time.sleep(1)
    
        return code

#RUN PIPELINE ON AWS CLUSTERING SYSTEM
#REQUIRES THAT A CLUSTER IS ALREADY CREATED AND FILES ARE ALREADY THERE
#PROBLEM: RESOURCES ARE NOT ALLOCATED PROPERLY AND ALL PROCESSESS GET EXECUTED SIMULATENOUSLY
def run_on_AWS(files_process_index, commands, log_file):
    import tempfile
    import time
    cluster_command = []
    #GENERATE COMMAND
    with tempfile.NamedTemporaryFile(delete=False) as temp:
        temp.write('#!/bin/bash\n')
        temp.write('#$ -S /bin/bash\n')
        temp.write('#$ -cwd\n')
        temp.write('#$ -j y\n')
        temp.write('#$ -N ' + os.path.basename(log_file) + "\n")
        temp.write('#$ -t ' + ",".join(files_process_index)+"\n")
        resources = []
        if (args.ram != None):
                resources.append('s_vmem=' + args.ram)
        if (args.cpu != None):
               resources.append('s_core=' + args.cpu)
        if (len(resources) > 0):
                temp.write('#$ -l '+ ",".join(resources))                    
                temp.write("\n")
        #print log_file + '\${JOB_ID}.log'
        temp.write('echo \"job initiated at $(date)\"\n')
        temp.write(" ".join(commands) + " \n")
        temp.write('echo \"job ended at $(date)\"\n')
        temp.flush()

        #RUN COMMAND        
        cluster_command.append("qsub")
        cluster_command.append(temp.name)
        proc = subprocess.Popen(cluster_command, stdout=subprocess.PIPE)
        output = proc.stdout.read()
        print output
        job_id_sub=output.split(" ")[2]
        job_id_sub=job_id_sub.split(".")[0]
        print "Submitted qsub ID: " + job_id_sub
    
    #CHECK STATUS
    not_done = True
    while (not_done):
        cluster_command = []
        cluster_command.append("qstat")
        proc = subprocess.Popen(cluster_command, stdout=subprocess.PIPE)
        output = proc.stdout.read()

        status = output.split("\n")
        if (len(status) <=2):
            not_done = False
        else:
            status = status[2:len(status)-1]
            running = []
            for i in range(0,len(status)):
                stat =  status[i]
                des = stat.split()
                job_id =  des[0]
                task_id = des[len(des)-1]
                if (job_id == job_id_sub):
                    running.append(task_id)
            if (len(running) == 0):
                not_done = False
            else :
                not_done = True
                print "Running task-id: \n" + ",".join(running)
                time.sleep(180)
                
    return 0

############################################
############################################
#HELPING FUNCTIONS##########################
############################################
def startProgress(title):
    global progress_x
    sys.stdout.write(title + ": [" + "-"*40 + "]" + chr(8)*41)
    sys.stdout.flush()
    progress_x = 0

def progress(x):
    global progress_x
    x = int(x * 40 // 100)
    sys.stdout.write("#" * (x - progress_x))
    sys.stdout.flush()
    progress_x = x

def endProgress():
    sys.stdout.write("#" * (40 - progress_x) + "]\n")
    sys.stdout.flush()

def ConfigSectionMap(section):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                DebugPrint("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1

#
def collect_input(input_user):

    #THROW ERROR IF INPUT DIR NOT EXISTEND
    input_dir = ConfigSectionMap("DIRECTORIES")['root_dir'] + "/" + input_user + "/raw"
    paired = False
    #GREP ALL FILES ENDING IN EITHER: BAM, FASTQ, FASTQ.GZ
    if (not os.path.exists(input_dir)):
        print_error("Input folder \'" + input_dir + "\' does not exist. Please create and fill it with raw data.")
    supported_extensions = ("*"+ConfigSectionMap("EXTENSIONS")["ext_bam"],   "*"+ ConfigSectionMap("EXTENSIONS")["ext_fastq"], "*"+ ConfigSectionMap("EXTENSIONS")["ext_fastq_gz"])
    files = []
    for ext in supported_extensions:
        i = input_dir + "/" + ext
        files.extend(glob.glob(i))
    if (len(files) == 0):
        print_error("No input files found in: " + input_dir)
    type = []

    files_index_holder = set()
    files_process_index = set()
    #FOR EACH GREPPED FILE GET: CELL_NUM, PAIRED
    for file in files:
        file_name = "."+str(file.split(os.extsep, 1)[0])
        file_name =  os.path.basename(file_name)
        split=file_name.split("#");
        no_file_name=split[len(split)-1]
        descr=no_file_name.split("_")

        if len(descr) == 2:
            type.append(descr[1])

    types=numpy.unique(type)
    logging.info("Found:\n" + "\n".join(files))
    for file in files:
        
        #DISTUINGISH CELL VS LANE VS PAIR
        cell_num = -1
        lane_num = -1 
        pair = -1
        
        #RECEIVE LAST PART OF PATH -> FILE NAME
        origin_file = file
        file = file.split("/")
        file = file[len(file)-1]
        file_name = ('.').join(file.split('.')[:-1])
        file_name =  os.path.basename(file_name)
        ext = "."+str(file.split(os.extsep, 1)[1])

        #GET CELL NUMBER (number after #)
        split=file_name.split("#");
        no_file_name=split[len(split)-1]
        descr=no_file_name.split("_")
        #IF ONE UNDERSCORE IT IS CELL NUMBER 
        if len(descr) == 1:
            cell_num = descr[0]
        #IF TWO UNDERSCORE -> CELL and Forward or Reverse
        elif len(descr) == 2:
            cell_num = descr[0]
            if ("1" in types and "2" in types and len(types) == 2):
                pair = descr[1]
                paired = True
            else:
                lane_num = descr[1]

        regexp = re.compile(r'\d+')

        index_place_holder = "#" + str(cell_num)
        if regexp.search(cell_num) is not None:
            files_process_index.add(int(cell_num))
            files_index_holder.add(origin_file.replace(index_place_holder, "#{#}"))
        read_length = get_read_length(files[0])
    
    logging.info("Paired:" + str(paired))
    logging.info("Read-length:" + str(read_length))
    return [files, files_process_index, files_index_holder, paired,read_length]

if __name__ == "__main__":
    # process command line args
    # TODO make a switch for various datatypes using **kwargs
    import argparse
    p = argparse.ArgumentParser()
    group1 = p.add_argument_group('Data processing', 'Modules to process your RNA sequencing data')
    group1.add_argument('-i','--input', help='Input directory')
    group1.add_argument('-m','--mapping', help='Mapping algorithm')
    group1.add_argument('-margs','--mapping_args', help='Mapping algorithm additional arguments', nargs="+")
    group1.add_argument('-q','--quantification', help='Quantification algorithm')
    group1.add_argument('-q_args','--quant_args', help='Quantification algorithm additional arguments', nargs="+")
    group1.add_argument('-o','--output', help='fastq file containing R1 pairs', default="output")
    group1.add_argument('--overwrite', dest="overwrite",help='Overwrite files', action="store_true")
    group1.set_defaults(overwrite=False)
    group1.add_argument('-g','--genome', help='Reference genome')
    group1.add_argument('-r', '--range', help='Index range to run')
    group2 = p.add_argument_group('Ref index', 'Reference genome index files')
    group2.add_argument('-f','--fasta', help='Fasta files')
    group2.add_argument('-gtf', help='Gene annotation file')
    
    p.add_argument('-c', '--config', help='Config file')
    p.add_argument('-ram','--ram')
    p.add_argument('-cpu','--cpu')
    p.add_argument('-l', '--cluster', help='Cluster type') 
    args = p.parse_args()
    
    Config = ConfigParser.ConfigParser()
    if (args.config == None):
        print_error("Config file not specified");
    Config.read(args.config)
    code = run(args)

