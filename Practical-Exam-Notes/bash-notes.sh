#! /bin/bash
# This file contains useful commands for the Practical Computing open book exam

# Navigating Files
pwd                 #show working path
ls -l               #show files as a list
ls -h               #show in human friendly format
ls *.txt            #* is a wildcard

# Manipulating Files
cat text1.txt text2.txt
cat text1.txt text2.txt > outputtext.txt
cat text1.txt text2.txt >> outputtext.txt
head outputtext.txt
# Unzip with gunzip and zip with gzip

# Text Processing
awk '{print $1}' mmp_mut_strains.txt | head -n 1 #returns column one (any amount of whitespace between)
awk '{print $1,3}' mmp_mut_strains.txt | head -n 1
awk '( $1=="coding_exon" && $11!="synonymous" )' mmp_mut_strains.txt
cut -f4 mmp_mut_strains.txt | head -n 10

# sort -k2 -n file

function names_numbers () {
  grep -rL "$1" media | grep -rh "^[0-9][0-9][0-9][0-9]\."
}

no_chr=$( cut -f4 $file | sed "1 d"| sort | uniq | grep -v '[+a-z]' | wc -l )
echo 'Number of chromosomes: ' $no_chr | tee -a script-outputs.txt

non_syn_effects=$( awk '( $8=="coding_exon" && $11!="synonymous" )'\
                   $file | cut -f11 | sort | uniq )
echo 'Types of non-synonymous mutations: ' $non_syn_effects
