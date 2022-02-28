#!/usr/bin/env bash

#arg 1: location to search of AA files
#arg 2: output name

find $1 -name "*_cycles.txt" | sort > scf.txt
find $1 -name "*_graph.txt" | sort > sgf.txt
cat scf.txt | rev | cut -d '/' -f 1 | cut -c12- | rev > san.txt
paste san.txt scf.txt sgf.txt > $2.input
rm san.txt scf.txt sgf.txt
