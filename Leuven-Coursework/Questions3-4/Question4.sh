#! /bin/bash

# This script is for question 4 of the first practical bioinformatics
# assignment at KU Leuven.

# Question: Write a script that, given a substance as an input argument,
# lists the names and numbers of the media that don't contain this substance.

function names_numbers () {
  grep -rL "$1" media | grep -rh "^[0-9][0-9][0-9][0-9]\."
}

names_numbers "Yeast extract"

function names_numbers_spaces () {
  substance="'$*'"
  grep -rL "$substance" media | grep -rh "^[0-9][0-9][0-9][0-9]\."
}

# names_numbers_spaces Yeast extract
