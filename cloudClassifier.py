#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May  3 09:26:58 2022

@author: rfaure
"""

from ete3 import NCBITaxa
import sys

ncbi = NCBITaxa()
ncbi.update_taxonomy_database()


file = "/home/rfaure/Documents/stage_M2/datasets/ATCC/results/results_with_taxa.txt"
fileOut = "/home/rfaure/Documents/stage_M2/datasets/ATCC/results/results_with_improved_taxa.txt"

f = open(file, "r")
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
    
    currentTag = line.split()[3].split('-')[0] 
    improvedTag = line.split()[3].split('-')[-1]
    improvedTag = '1' #non-deconvoluted
    #currentTag = line.split()[3] #deconvoluted
    
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
        for read in reads :
            predictedtaxon = int(read.split()[1]) #this is the taxon predicted by kraken 
            lowertaxon = int(read.split()[1])
            improvedID = read.split()[3].split('-')[-1]
            improvedID = 1 #for testing unconvoluted barcodes
            
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
                                
            if read.split()[0] == "SRR12283286.1.125015170" :
                print(descent)
                print("Read SRR12283286.1.125015170 : ", lowertaxon)
            
            fo.write(read.strip("\n")+"\tIT:"+str(lowertaxon)+"\n")
        
        oldTag = currentTag    
        lineCount = 0
        cloudTaxa = set()
        reads = []
        nbtags += 1
        #print("Processed ", nbtags, " tag\r")
        # if nbtags == 5 :
        #     break;
            
    else :
        
        reads += [line]
        taxon = int(line.split()[1])
        cloudTaxa.add((taxon, improvedTag))
        lineCount += 1
        
print("I went down to the strain level for ", foundStrains, " read, and to species level for ", foundSpecies, " reads")
    
        