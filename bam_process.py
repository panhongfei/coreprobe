#!/usr/bin/python
import sys
import pysam
from collections import defaultdict

##
def nonredundant_list( old_list, new_item ):
    if new_item not in old_list:
        old_list.append( new_item )
    return old_list

if len( sys.argv ) != 3:
  print "Arguments needed!!\n" \
       +"usage: python bam_process.py [input_file] [output_file]\n" \
       +"where input_file is a bamfile with index"
  exit()

#here a sorted indexed bam/sam file is needed
#input_file = "/home/huang/Documents/my_research/proceeding/practice/sim-Arthrobacter/test6-sorted.bam"
input_file = sys.argv[1]
#output file
#output_file = "/home/huang/Documents/my_research/proceeding/practice/sim-Arthrobacter/test-c++/test-stat5-length"
output_file = sys.argv[2]

#the following two list contains non-redundant reference names and their lengths
refspecies_list = list()
refspecies_lengths = list()
#the reads_list contains the counts of each read
reads_list = defaultdict(int)
#the assignments_list contains the assignment of each read, 
#where the indecies are consistent with those in refspecies_list and reads_list
assignments_list = defaultdict( list )


try:
  samfile = pysam.AlignmentFile( input_file )
  handle = samfile.fetch(until_eof=True)
  print "start to process the input bam/sam file..."
  for alignment in handle:
    if not alignment.is_unmapped:
      reads_list[alignment.qname] += 1
      refname = samfile.getrname( alignment.reference_id )
      #nonredundant_list( refspecies_list, refname )
      if refname not in refspecies_list:
        refspecies_list.append( refname )
        refspecies_lengths.append( samfile.lengths[ alignment.reference_id ] )    
      assignments_list[ alignment.qname ].append( refspecies_list.index( refname ) )

  samfile.close()
except Exception,e:
  print "Error!!\n"
  print e
  exit()


refspecies_list = [ x for x in refspecies_list ]

print "start to write to the output file..."
with open( output_file, "w" ) as my_file:
    my_file.write( str(len(refspecies_list)) )
    
    my_file.write( "\n" )  
    for ref_name in refspecies_list:
      my_file.write( str(ref_name) + "\t" )
    
    my_file.write( "\n" )  
    for ref_length in refspecies_lengths:
      my_file.write( str(ref_length) + "\t" )
          
    my_file.write( "\n" )
    for read_cluster, read_cluster_assignments in assignments_list.iteritems():
        for assignment in read_cluster_assignments:
            my_file.write( #str(read_cluster) +"\t"+
                           str( assignment ) + "\t" )

    my_file.write( "\n" )
    for read_cluster, count in reads_list.iteritems():
        my_file.write( #str(read_cluster)+"\t"+
                       str( count ) + "\t" )
     
my_file.close()
print "process finished, please use PraBam to estimate relative abundance"
