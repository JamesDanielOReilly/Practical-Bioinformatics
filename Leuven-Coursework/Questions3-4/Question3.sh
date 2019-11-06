#! /bin/bash

# This script is for question 3 of the first practical bioinformatics
# assignment at KU Leuven.

file='mmp_mut_strains.txt'     # file variable

# Question: How many chromosomes are there in the dataset?
# Solution: Return the number of unique elements in the chromosome column

no_chr=$( cut -f4 $file | sed "1 d"| sort | uniq | grep -v '[+a-z]' | wc -l )
echo 'Number of chromosomes: ' $no_chr

# Question: How many mutations are on introns and exons?
# Solution: Return the number of rows with "intron" and "coding_exon" in the feature column

mut_intron=$( awk '( $8=="intron" )' $file | wc -l )
echo 'Number of mutations on introns: ' $mut_intron

mut_exon=$( awk '( $8=="coding_exon" )' $file | wc -l)
echo 'Number of mutations on exons: ' $mut_exon

# Question: What type of non-synonymous mutations do you find in exons?
# Solution: First get only the rows with exons. Then return the unique entries in the
# effect column that or not synonymous.

non_syn_effects=$( awk '( $8=="coding_exon" && $11!="synonymous" )'\
                   $file | cut -f11 | sort | uniq )
echo 'Types of non-synonymous mutations: ' $non_syn_effects
