#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  3 09:26:58 2022

@author: rfaure
"""

from ete3 import NCBITaxa
import sys
import argparse
import os


def parse_args() :
    
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fastq", help="Fastq of the barcoded reads. The reads MUST be sorted by barcode")
    parser.add_argument("-c", "--classification", help="File containing classification of the reads, following the default kraken2 output format.\
                        The reads should be in the same order as in the fastq file")
    parser.add_argument("-o", "--output", help="Output file, in the default kraken2 output format")
    
    return parser.parse_args()


def main() :
    
    args = parse_args()
    file = args.classification
    fastq = args.fastq
    fileOut = args.output
    
    #now pre-treat all the data to attach barcode to kraken reads
    print("Pre-processing the files...")
    command = "awk '{if (NR%4==1) print;}' " + fastq + " > tmp889.txt"
    print("command : ", command)
    os.system(command)
    
    command = "paste " + file + " tmp889.txt > tmp887.txt"
    print("command : ", command)
    os.system(command)
    
    print("Now all reads have been attached to their barcodes")
    print("Moving on to loading NCBI taxonomy")
    
    ncbi = NCBITaxa()
    ncbi.update_taxonomy_database()
    
    print("Loaded taxonomy, let's improve the taxonomic assignment with barcodes now")
    
    # file = "/home/rfaure/Documents/stage_M2/datasets/ATCC/results/results_with_taxa.txt"
    # fileOut = "/home/rfaure/Documents/stage_M2/datasets/ATCC/results/results_with_improved_taxa.txt"
    
    args_command = parse_args()
    

    f = open("tmp887.txt", "r")
    fo = open(fileOut, "w")
    reads_nb = 0
    lineCount = 0
    oldTag = ''
    currentTag = ''
    cloudTaxa = set()
    reads = [] 
    nbtags = 0
    
    foundSpecies = 0
    foundStrains = 0
    
    for line in f :
        
        ls = line.split()
        currentTag = ""
        improvedTag = "0"
        for field in ls :
            if "BX:Z:" in field[:6] or "BC:Z:" in field[:6] :
                
                currentTag = field.split('-')[0] 
                improvedTag = field.split('-')[-1]
                #improvedTag = '1' #non-deconvoluted
                
        if currentTag != "" and line[0] == "C" :
        
            if currentTag != oldTag : #we're changing tag, it's time to see what happened        
            
            
                alreadyseen = set()
                #first build the taxonomic tree of this barcode
                descent = {} #list of children for each node with their improved tag
        
                for (leaf, improved) in cloudTaxa :
                  
                    if (leaf, improved) not in alreadyseen :
                        
                        alreadyseen.add((leaf, improved))
                        lineage = ncbi.get_lineage(leaf)[::-1]#go through the lineages from the leafs
                        
                        pastAncestor = leaf
                        if leaf not in descent :
                            descent[leaf] = {}
                        
                        for ancestor in lineage[1:] : 
                            
                            alreadyseen.add(ancestor)
                            
                            if ancestor not in descent :
                                descent[ancestor] = {}
                            
                            if improved not in descent[ancestor] :
                                descent[ancestor][improved] = set()
                                
                            descent[ancestor][improved].add(pastAncestor)
                            pastAncestor = ancestor
                
                #print(descent)
        
                #now start the rescuing of the reads
                for (read, improvedID) in reads :
                    predictedtaxon = int(read.split()[2]) #this is the taxon predicted by kraken 
                    lowertaxon = int(read.split()[2])
                    #improvedID = 1 #for testing convoluted barcodes
                    
                    cont = True
                    while cont :
                        if len(descent[lowertaxon]) == 1 : #only one improved barcode
                            for j in descent[lowertaxon].keys() :
                                if len(descent[lowertaxon][j])==1 :
                                    for k in descent[lowertaxon][j] : 
                                        lowertaxon = int(k)
                                else :
                                    cont = False
                                # print("lower taxon : ", lowertaxon, " ", descent[lowertaxon], " ", list(list(descent[lowertaxon])[0])[0])
                                # lowertaxon = int( list(list(descent[lowertaxon])[0])[0])
        
                        elif len(descent[lowertaxon]) == 0 : #no more descent
                            cont = False
                        elif improvedID in descent[lowertaxon] : #then we know where to look
                            if len(descent[lowertaxon][improvedID]) == 1 : #only one possibility in that barcode
                                lowertaxon =  int(list(descent[lowertaxon][improvedID])[0])
                            else : #several possibilities in that barcode too
                                cont = False
                        else : #hesitation between several non-improved descendents
                            #check if they all agree
                            candidateLower = set()
                            for barcode in descent[lowertaxon].keys() :
                                for l in descent[lowertaxon][barcode] :
                                    candidateLower.add(l)
                            if len(candidateLower) == 1 :
                                lowertaxon = list(candidateLower)[0]
                            else :
                                cont = False
                                
                    #now we know the new taxon, it's lowertaxon
                    
                    ls = read.split('\t')
                    newline = ls[0] + '\t' + ls[1] + '\t' + str(lowertaxon) + '\t' + ls[3] + '\t' + ls[4] + '\n'
                    
                    fo.write(newline)
                
                oldTag = currentTag    
                lineCount = 0
                cloudTaxa = set()
                reads = []
                nbtags += 1
                #print("Processed ", nbtags, " tag\r")
                # if nbtags == 5 :
                #     break;
                    
            else :
                
                reads += [(line, improvedTag)]
                taxon = int(line.split('\t')[2])
                cloudTaxa.add((taxon, improvedTag))
                lineCount += 1
        
        else : #no barcode information here, just output the kraken2 line
            ls = line.split()
            newline = ""
            append = True
            for field in ls :
                if field[0] == "@" :
                    append = False
                if append :
                    newline += field + "\t"
            newline += "\n"
            fo.write(newline)

        
    #print("I went down to the strain level for ", foundStrains, " read, and to species level for ", foundSpecies, " reads")
    os.system("rm tmp889.txt")
    os.system("rm tmp887.txt")
    
    print("Done !")

if __name__ == "__main__":
    main()
    
        