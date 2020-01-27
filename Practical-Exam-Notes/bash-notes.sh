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
