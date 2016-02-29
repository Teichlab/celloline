#!/usr/bin/python
import sys 
import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt 
from os import listdir
from os.path import isfile, join
import seaborn as sb
import math
import scipy
import getopt
from collections import Counter, defaultdict, OrderedDict
from Bio import SeqIO
import glob
import fnmatch
import os

from intervaltree import Interval, IntervalTree
#LOAD GTF FILE
def load_GTF(gtf_file):

    gtf_index = defaultdict()
    with open(gtf_file) as f:
        for line in f:
             if (not line.startswith("#")):
                 entry = line.split("\t")
                 entry_addition = entry[8]
                 entry_addition = entry_addition.split(";")
                 entry_addition = entry_addition[0].split(" ")
                 gene_id = entry_addition[1]
               
                 type = entry[2] 
                 #TYPE(Gene, exon etc.), START, END, STRAND, gene_ID
                 info = [type, entry[3], entry[4], entry[6], gene_id]
        
                 #Build GTF INDEX
                 if (type != "" and entry[3]!= entry[4]):
                    index = IntervalTree()
                    if (entry[0] in gtf_index):
                         index = gtf_index[entry[0]]
                    index.addi(int(info[1]),int(info[2]),info) 
                    gtf_index[entry[0]] = index

    return (gtf_index)
 
def generate_mapping_stats(input_file, output_file, gtf_file, sample_name, count):

    print("Loading GTF file...") 
    gtf_dict = load_GTF(gtf_file) 
    print("Loaded.")

    #OUTPUT TABLE CONTAING STATS
    output_table = OrderedDict()
        
    #Dict indicating to which genes a specific read maps to
    #It is a temporary dict
    reads_mapped_to = defaultdict(str)
    exonic_mappings_temp = defaultdict(str)
    #Dict indicating which read is multi-mapped
    #It is a temporary dict
    multi_maps = defaultdict(str)
    sam_files = [input_file]
    exonic_multi_table = defaultdict(str)

    #MAPPABILITY
    output_table["total"] = 0 
    output_table["mapped"] = 0
    output_table["unmapped"] = 0
    output_table["unique"] = 0
    output_table["multi"] = 0
    

    #Transcript-annotated frequ
    output_table["intergenic"] = 0
    output_table["intragenic"] = 0
    output_table["exonic"] = 0
    output_table["intronic"] = 0
    output_table["ambigious"] = 0
    output_table["exonicM"] = 0

    output_table["alignments"] = 0
    output_table["multi-intergenic"] = 0
    output_table["multi-intragenic"] = 0

    output_table["multi-exonic"]=0
    output_table["multi-intronic"]=0
    output_table["multi-ambigious"]=0

    #ERROR
    output_table["perfect"] = 0
    output_table["partly_perfect"] = 0
    output_table["mapped_no_correct"] = 0
    for i in range(0,10):
        output_table["S_"+str(i)] = 0
    output_table["S_10+"] = 0
    output_table["I"] = 0
    output_table["D"] = 0
    output_table["INDEL"] = 0

    reads = Counter()
    multi_reads = defaultdict(str)

    exonic_multi_output = open(output_file+".exonic",'w')    
    #SAM PARSE SAM FILE
    for sam in sam_files:
        print("Parsing sam file...")
        with open(sam) as f:

            for line in f:
                split=line.split("\t")
                if (not line.startswith("@PG") and not line.startswith("@HD") and not line.startswith("@SQ") and len(split) >= 10):
                    read_name=split[0]
                    flagCode = int(split[1])
                    chrom = split[2]
                    pos = split[3] 
                    errors=split[5]
                    read=split[9]
                
                    errors_a = list(errors)
                    number = ""
                    num = 0
                    error_table = defaultdict(int)
                    name_and_flag = read_name
                    
                    #CHECK IF READ MAPPED OR UNMAPPED 
                    #IT US UNMAPPED
                    if(flagCode & 0x0004 != 0):
                        output_table["unmapped"] += 1
                        output_table["total"] += 1
                        error_table["*"] += 1
                    #IT IS MAPPED
                    else:
                        if (flagCode & 0x0001 !=0):                             #This is paired end sequencing
                            if (flagCode & 0x0040 != 0):                        #1st read
                                name_and_flag += ";first"
                            if (flagCode & 0x0080 != 0):                        #2nd read
                                 name_and_flag  += ";second"
                    
                        # CHECK TO WHICH GENE(S) IT MAPPED TO
                        genes_info, num_genes, num_exons = get_gene(gtf_dict, [chrom, pos])                           
                        #GENE COUNTS: only NON-overlapping genes and if read within exons
                        if (count and num_genes == 1 and num_exons > 0):
                          info = genes_info[0]                            
                          gene_id = info[4]
                          mapped_to = []
                          if (name_and_flag in reads_mapped_to):
                              mapped_to = reads_mapped_to[name_and_flag]
                          mapped_to.append(gene_id)
                          reads_mapped_to[name_and_flag] = mapped_to

                        output_table["alignments"] += 1.0
                        #STATS
                        if(name_and_flag not in reads):
                            reads[name_and_flag] += 1
                            output_table["unique"] += 1
                            output_table["total"] += 1
                            output_table["mapped"] += 1
                            

                            if(num_genes == 0):
                                output_table["intergenic"] += 1
                            elif (num_genes == 1):
                                output_table["intragenic"] += 1
                                if (num_exons==0):
                                    output_table["intronic"] += 1
                                else:
                                    output_table["exonic"] += 1 
                                    d = []
                                    if (name_and_flag in exonic_mappings_temp):
                                        d = exonic_mappings_temp[name_and_flag]
                                    d.append([genes_info[0], chrom, pos])
                                    exonic_mappings_temp[name_and_flag] = d
                            elif (num_genes > 1):
                                output_table["ambigious"] += 1
    
                        #READ IS MULTI-MAPPED
                        else:
                            if(reads[name_and_flag] == 1):
                                output_table["unique"] -= 1
                                output_table["multi"] += 1
                            reads[name_and_flag] += 1
                            d = []
                            #GET KNOWLEDGE IF FIRST MAPPING EXONIC OR INTRONIC
                            if (name_and_flag in exonic_mappings_temp):
                                d = exonic_mappings_temp[name_and_flag]
                            #output_table["alignments"] += 1.0 
                            if(num_genes == 0):
                                output_table["multi-intergenic"] +=(1)
                            elif (num_genes == 1):
                                output_table["multi-intragenic"] += (1)
                                if (num_exons==0):
                                    output_table["multi-intronic"] += (1)
                                else:
                                    output_table["multi-exonic"] += (1) 
                                    d.append([genes_info[0], chrom, pos])
                            elif (num_genes > 1):
                                output_table["multi-ambigious"] += (1)
                            #IF AT LEAST ONE EXONIC ALIGNMENT
                            if (len(d) > 0):
                                exonic_multi_table[name_and_flag] = d 
                        #PARSE MAPPING ERRORS         
                        for i in errors_a:
                            if (re.match("[0-9]",i)):
                                number +=(i)
                            elif(re.match("[A-Z]",i)):
                                num = int(number)
                                error_table[i] += num 
                                number = ""
                        #print output_table
                        #TABLE OF HOW MANY READS MAP PERFECT, PARTLY PERFECT, SUBSTITUINTS ETC
                        if("M" in  error_table and len(error_table)==1):
                            output_table["perfect"] += 1
                        elif("M" in error_table and len(error_table) > 1): 
                            output_table["partly_perfect"] += 1
                        elif("M" not in error_table and "*" not in error_table):
                            output_table["mapped_no_correct"] += 1
    
                        if("S" in error_table):
                            if(int(error_table["S"]) < 10):
                                output_table["S_"+str(error_table["S"])] += 1
                            else:
                                output_table["S_10+"] += 1
                        elif("S" not in error_table):
                            output_table["S_0"] += 1

                        if("I" in error_table):
                            output_table["I"] += 1

                        if("D" in error_table):
                            output_table["D"] += 1
    
                        if("I" in error_table or "D" in error_table):
                             output_table["INDEL"] += 1

  
    exonic_multi_output.close() 
    #WHEIGHT COUNTS  
    if (count):
        counts, counts_unique, counts_multi = weight_counts(reads_mapped_to)
        write_counts(output_file+".counts.unique", counts_unique, sample_name)
        write_counts(output_file+".counts", counts, sample_name)
        write_counts(output_file+".counts.multi", counts_multi, sample_name)

    o = ""
    exonicM = len(exonic_multi_table.keys())
    output_table["exonicM"] = exonicM
    write_stats(output_file, output_table, sample_name)
     
def write_stats(output_file, stats_table, o):
    #OUTPUT STATS
    f = open(output_file,'w')
    for k,v in stats_table.items():
        if (str(k) in ["unique", "multi", "intragenic", "intergenic", "exonic", "intronic", "ambigious", "exonicM"]):
            v = (v+0.0) / (stats_table["mapped"]+0.0)
            v = '%.2f' % (100.0*(v))
        if (str(k) in ["multi-intragenic", "multi-intergenic", "multi-exonic", "multi-intronic", "multi-ambigious"]):
            
            v = (v+0.0)
            if (stats_table["alignments"] !=0):
                v = v / (stats_table["alignments"]+0.0)  
            v = '%.2f' % (100.0*(v))
        val = v 
        print str(k) + ":" + str(val)
        o += "," + str(val)
    o += "\n" 
    f.write(o)
    f.close()


def write_counts(output_file, counts,o):
    f = open(output_file,'w')
    o="Gene,"+o
    for k,v in counts.items():
        o = str(k).replace("\"", "")+ "," + str(v) + "\n"
        f.write(o)
    f.close()


def weight_counts(reads_mapped_to):

    #Considers also multi-mapped reads
    #If multi-mapped, counts are shared across genes
    #Only considers genes that do not have any overlapping genes
    #Considers paired-end as single ended reads
    counts = defaultdict(float)
    counts_unique = defaultdict(float)
    counts_multi = defaultdict(float) 
    #OUTPUT GENE COUNTS
    for read in reads_mapped_to.keys():
        genes = reads_mapped_to[read]
        if (len(genes) == 1):
            counts_unique[genes[0]] += 1.0
            counts[genes[0]] += 1.0
        else:
            for i, g in enumerate(genes):
                counts_multi[g] += 1.0
                counts[g] += 1.0/len(genes)
    return counts, counts_unique, counts_multi


def count_genes(genes_info, reads_mapped_to, name_and_flag):
#GENE COUNTS: only NON-overlapping genes

    if (len(genes_info)== 1):
        info = genes_info[0]                            
        gene_id = info[4]
        mapped_to = []
        if (name_and_flag in reads_mapped_to):
            mapped_to = reads_mapped_to[name_and_flag]
        mapped_to.append(gene_id)
        reads_mapped_to[name_and_flag] = mapped_to

def get_gene(gtf_dict, pos_pair):
    
    num_genes = 0
    num_exons = 0
    
    if (pos_pair[0] not in gtf_dict):
        print ("Ignored pos: "+pos_pair[0])
        return([[], num_genes, num_exons])

    entries = gtf_dict[pos_pair[0]]
    pos = pos_pair[1]

    found = []
    found = entries.search(int(pos_pair[1]))

    list = []
    for e in found:
        info = e[2]
        if(info[0]== "gene"):
            list.append(info)
            num_genes+= 1
        elif (info[0] == "exon"):
            num_exons+=1
    return([list, num_genes, num_exons])

def main(argv):    
    input_file = ""
    output_dir = ""
    gtf_file = ""
    sample_name = ""
    count = False
    try:    
        opts, args = getopt.getopt(argv, "i:o:n:g:c:")
    except getopt.GetoptError:    
        sys.exit(2)    
        usage()
    for opt, arg in opts:    
        if opt in ("-i"):    
            input_file = arg    
        elif opt in ("-o"): 
            output_dir = arg 
        elif opt in ("-n"): 
            sample_name = arg 
        elif opt in ("-g"):
            gtf_file = arg 
        elif opt in ("-c"):
            count = True
    if input_file != "" and  output_dir != "" and gtf_file != "" and sample_name != "":  
        generate_mapping_stats(input_file, output_dir, gtf_file, sample_name, count)    
    else: 
        print ("Argument missing")
main(sys.argv[1:])
