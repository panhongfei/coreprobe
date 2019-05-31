#!/bin/bash
#$1 is the input samfile name
echo $1
#$2 is the output bamfile name
echo $2

samtools view -bS $1 -o $1.bam
samtools sort $1.bam $1-sorted
mv $1-sorted.bam $2

samtools index $2
