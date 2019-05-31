#!/usr/bin/python
# -*- coding: utf-8 -*-

#This script generates reference sequences of each ref-taxa/species
#It extract gene sequences from the all-gene-sequence file downloaded from MetaRef
#according to a specific taxa/species gene list (can be a core genes list etc.)
#The output is a fna file contains all ref-taxa/species sequences
import sys
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from collections import defaultdict

if len( sys.argv ) != 4:
  print "Arguments invalid!!\n" \
       +"usage: python species_refseqs_generator.py [all_genes_seqs_file] [taxa_genes_list_file] [output_file]\n" \
       +"where all_genes_seqs_file is the file contains sequences of all different genes\n" \
       +"where taxa_genes_list_file is the sorted file whose each line has information of" \
       +"\"[taxa name]\t[gene index]\""# \
       #+"where is_iteration is a bool to choose iteration way"
  exit()
          
#which contains all metagenes sequences
#refgenes_file = "/home/huang/文档/my_research/data/metaRef/DNA-sequences-of-all-centroids/metaref.org/static/centroids_v.1.0.fna"
refgenes_file = sys.argv[1]
#which contains non-redundant sorted refgene ID
#species_gene_file = "/home/huang/文档/my_research/proceeding/practice/sim-Arthrobacter/test-taxa-AllGeneIDs-sorted"
species_gene_file = sys.argv[2]
#output_file = "/home/huang/文档/my_research/proceeding/practice/sim-Arthrobacter/test-refgenes-MC-comb-2.fna"
output_file = sys.argv[3]

#is_iteration = sys.argv[4]

#a defaultlist contains species names as keys, 
#and the list of gene indecies as value of each key
species_gene_dict = defaultdict(list)

try:
  with open(species_gene_file) as f:
     for line in f:
        if line.strip():
            ref_info = line.strip().split("\t")
            species_gene_dict[ref_info[0]].append(ref_info[1])
except Exception,e:  
  print "Error in processing species_genes_list_file!!\n"
  print e         
  exit()  
#print species_gene_dict

#According to species and their genes indecies list, iterate refgenes sequence file
#to extract a complete reference sequence (can be core-genomes etc.) for each species
#The sequence of a species will be yielded once all genes in the species list included
def refseq_iterator(refgenes_file, species_genes, index):
    print species_genes
    for species, genes in species_genes.iteritems():
        output_seq = ""
        handle = open(refgenes_file, "rU")
        for record in SeqIO.parse(handle, "fasta"):   
            if len(genes) == 0:
                species_genes[species].append('done')
                break        
            if record.id in genes:
                output_seq = output_seq + str(record.seq)
                genes.remove(record.id)
        yield SeqRecord(Seq(output_seq, generic_dna), id = str(species).replace(" ", "_"), description = str(index))
        print species_genes
        index += 1
        
def refseq_iterator2(refgenes_file, species_genes, index):
    
        output_seq = ""
        handle = open(refgenes_file, "rU")
        print species_genes
        for record in SeqIO.parse(handle, "fasta"): 
          for species, genes in species_genes.iteritems(): 
            if record.id in genes:
                output_seq = output_seq + str(record.seq)
                genes.remove(record.id)
            if len(genes) == 0:
                yield SeqRecord(Seq(output_seq, generic_dna), id = str(species).replace(" ", ""), description = str(index))
                species_genes[species].append('done')
                print species_genes
        
        index += 1

def refseq_list(refgenes_file, species_genes, index):
    for species, genes in species_genes.iteritems():
        output_list = list()
        output_seq = ""
        handle = open(refgenes_file, "rU")
        for record in SeqIO.parse(handle, "fasta"):   
            if len(genes) == 0:
                break        
            if record.id in genes:
                output_seq = output_seq + str(record.seq)
                genes.remove(record.id)
        output_list.append( SeqRecord(Seq(output_seq, generic_dna), id = str(species).replace(" ", ""), description = str(index)) )
        index += 1
    return output_list

try:
  output_handle = open(output_file, "w")
  #if is_iteration:
  SeqIO.write(refseq_iterator(refgenes_file, species_gene_dict, 1), output_handle, "fasta")
  #else:
  #SeqIO.write(refseq_list(refgenes_file, species_gene_dict, 1), output_handle, "fasta")
  output_handle.close()
except Exception,e:
  print "Error in processing all_genes_seqs_file to generate refseqs!!\n"
  print e
  exit()
