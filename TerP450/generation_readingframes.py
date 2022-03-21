#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 16:33:09 2022

@author: friederike
"""

import re
import itertools

#listreadingframes:
minlengthreadingframe=20
maxlengthreadingframe=40
aa_in_core=7
aromatic_amino_acids_forward=('CAT',"CAC","TTT","TTC","TAT","TAC","TGG")
aromatic_amino_acids_reverse=("CCA","GTA","ATA","GAA","AAA","GTG","ATG")
filenamereadingframes="Output/readingframes.txt"
listreadingframe=[]

stopcodon1=('TGA', 'TAA', 'TAG')
startcodon1=('ATG',"GTG")
stopcodon2=("TCA","TTA","CTA")
startcodon2=("CAT","CAC")

#source=https://parts.igem.org/Ribosome_Binding_Sites/Catalog ,http://arep.med.harvard.edu/labgc/adnan/projects/Utilities/revcomp.html
ribosomalbindingsites1=["AAAGA[A-Z]{3}GA[A-Z]{3}","AACTAGAATCACCTCTTGCTTTTGGGTAAGAAAGAGGAGA","AACTAGAATCACCTCTTGGATTTGGGTATTAAAGAGGAGA","ATTAAAGAGGAGAAA","TTCACACAGGAAACC",'GTGTG','GTGTGTCTAG','TCACACAGGAAACCGGTTCGATG','TCACACAGGAAAGGCCTCGATG','TCACACAGGACGGCCGGATG','TCTCACGTGTGTCAAG','TCTCACGTGTGT','CATCCCT','TCACATCCCT','TCACATCCCTCC','ACTGCACGAGGTAACACAAG','TACGAGGAGGATGAAGAGTA','ACTTTACTTATGAGGGAGTA','ACGAAGACGGAGACTTCTAA','AACCCTCAGGAGGTAAACCA','AAGACATGGAGACACATTTA','ACTGCACGAGGTAACACAAG','GAGGAGGATGAAGAGTAATGTGAAGAGCTG','AGGAGGTCATC','GCAAGCTCTTTTTTCAGTTGTCTC','CTGATAGTTAAAATCACCAGCATGA','TAAAAACAAGAGGAAAACAA','TCTCCTCTTT','ACGGAGAAGCAGCGAA','GAGGTTGGGACAAG','TAAATGTATCCGTTTATAAGGACAGCCCGA','CTCTTAAGTGGGAGCGGCT','CTCTACCGGAGAAATT','CTCATCGTTAAAGAGCGACTAC','CTCAGCCTGTACCTGGAGAGCCTTTC','CTCAAGGAGG','GAGAGG','AGGAGGATTACAA','AAAGAGGAGAAA','TCACACAGGAAAG','GGAAGAGG','TTTCTCCTCTTTAAT','TCACACAGGAAAGGCCTCG','ATTAAAGAGGAGAAATTAAGC',' TCGTTTCTGAAAAATTTTCGTTTCTGAAAA','TGGCTAACATAGGGT','TGGCTAACTGAGGAT','TGGCTAACCCAGGGT','TGGCTAACTCAGGTG','TGGCTAACCCTGGTA','TGGCTAACTTGGGAC','TGGCTAACGCAGGTC','TGGCTAACATCGGTG','TTAATTAAGGAAAAGATCT','CAGAAGAGGATATTAATA','TTGATAAGGAATTGTA','TCAGAGGAGATAATTTA','TGACACGTTGAGCGGTATGA','ACAGATAACAGGAGTAAGTA','TAAAGGGAGAAAAAT','GAGTCTTGAGGTAACTAT','TCAGGAATATTAAAAACGCT','ATTTGAAGGAAAATATT','CAAAAACATACTGCAGGAAT','TGCCATTGCAAAGGAGAAGACT','AAGGGGGAATTCAAAT','AAGGGGTGCAGAAT','AGGTGGAATCACAG','ATAGATAAAAATGGTAACAAT','GGGATATAGCCTGAGGGGCCTGTA','CGGCAATAACAGAGGCGATTT','ATTAAAGAGGAGAAATA','TCACACAGGAAAGTA','AAAGGAGGTGT','AGAGGTGGTGT','AGGAGG','GAGG','TAAAGGAGGAA','AAAGGTGGTGAA','AGGAAACAGAACC','ATATTAAGAGGAGGAG','AGAGAACAAGGAGGGG','GATTGGGATAAATAAT','ATCAACCGGGGTACAT','TTTGGAGATTTTCAAC','AAAAAAGGTAATTCAA','CATAAGGTAATTCACA','ATAAGGAGTCTTAATC','GTTCCGGCTAAGTAAC','TAATGGAAACTTCCTC','TCGCTGGGGGTCAAAG','ATTTGAGGGGGATTCA','AATTTAGGTCAGAAG','AATCAATAGGAGAAATCAAT','TTAAAGAGGAGAAATACTAG']
ribosomalbindingsites2=["[A-Z]{3}TC[A-Z]{3}TCTTT","TCTCCTCTTTCTTACCCAAAAGCAAGAGGTGATTCTAGTT","TCTCCTCTTTAATACCCAAATCCAAGAGGTGATTCTAGTT","TTTCTCCTCTTTAAT",'GGTTTCCTGTGTGAA','CACAC','CTAGACACAC','CATCGAACCGGTTTCCTGTGTGA','CATCGAGGCCTTTCCTGTGTGA','CATCCGGCCGTCCTGTGTGA','CTTGACACACGTGAGA','ACACACGTGAGA','AGGGATG','AGGGATGTGA','GGAGGGATGTGA','CTTGTGTTACCTCGTGCAGT','TACTCTTCATCCTCCTCGTA','TACTCCCTCATAAGTAAAGT','TTAGAAGTCTCCGTCTTCGT','TGGTTTACCTCCTGAGGGTT','TAAATGTGTCTCCATGTCTT','CTTGTGTTACCTCGTGCAGT','CAGCTCTTCACATTACTCTTCATCCTCCTC','GATGACCTCCT','GAGACAACTGAAAAAAGAGCTTGC','TCATGCTGGTGATTTTAACTATCAG','TTGTTTTCCTCTTGTTTTTA','AAAGAGGAGA','TTCGCTGCTTCTCCGT','CTTGTCCCAACCTC','TCGGGCTGTCCTTATAAACGGATACATTTA','AGCCGCTCCCACTTAAGAG','AATTTCTCCGGTAGAG','GTAGTCGCTCTTTAACGATGAG','GAAAGGCTCTCCAGGTACAGGCTGAG','CCTCCTTGAG','CCTCTC','TTGTAATCCTCCT','TTTCTCCTCTTT','CTTTCCTGTGTGA','CCTCTTCC','ATTAAAGAGGAGAAA','CGAGGCCTTTCCTGTGTGA','GCTTAATTTCTCCTCTTTAAT','TTTTCAGAAACGAAAATTTTTCAGAAACGA','ACCCTATGTTAGCCA','ATCCTCAGTTAGCCA','ACCCTGGGTTAGCCA','CACCTGAGTTAGCCA','TACCAGGGTTAGCCA','GTCCCAAGTTAGCCA','GACCTGCGTTAGCCA','CACCGATGTTAGCCA','AGATCTTTTCCTTAATTAA','TATTAATATCCTCTTCTG','TACAATTCCTTATCAA','TAAATTATCTCCTCTGA','TCATACCGCTCAACGTGTCA','TACTTACTCCTGTTATCTGT','ATTTTTCTCCCTTTA','ATAGTTACCTCAAGACTC','AGCGTTTTTAATATTCCTGA','AATATTTTCCTTCAAAT','ATTCCTGCAGTATGTTTTTG','AGTCTTCTCCTTTGCAATGGCA','ATTTGAATTCCCCCTT','ATTCTGCACCCCTT','CTGTGATTCCACCT','ATTGTTACCATTTTTATCTAT','TACAGGCCCCTCAGGCTATATCCC','AAATCGCCTCTGTTATTGCCG','TATTTCTCCTCTTTAAT','TACTTTCCTGTGTGA','ACACCTCCTTT','ACACCACCTCT','CCTCCT','CCTC','TTCCTCCTTTA','TTCACCACCTTT','GGTTCTGTTTCCT','CTCCTCCTCTTAATAT','CCCCTCCTTGTTCTCT','ATTATTTATCCCAATC','ATGTACCCCGGTTGAT','GTTGAAAATCTCCAAA','TTGAATTACCTTTTTT','TGTGAATTACCTTATG','GATTAAGACTCCTTAT','GTTACTTAGCCGGAAC','GAGGAAGTTTCCATTA','CTTTGACCCCCAGCGA','TGAATCCCCCTCAAAT','CTTCTGACCTAAATT','ATTGATTTCTCCTATTGATT','CTAGTATTTCTCCTCTTTAA']

#create different shine-dalgarno sites:
dalgarno1=[]
dalgarno2=[]

#also look for common cores
commoncore=["TCCTGGTACATCTGGTACTGA","GCCTGGTACCTCTGGTACTGA","GCCTGGTACATCTGGTACTGA","GCCTGGTACCACTGGTACTAA"]
commoncorerev=["TCAGTACCAGATGTACCAGGA","TCAGTACCAGAGGTACCAGGC","TCAGTACCAGATGTACCAGGC","TTAGTACCAGTGGTACCAGGC"]
def generate_necessary_motifs(listreadingframe):
    forward_aa_list=list(itertools.product(aromatic_amino_acids_forward, repeat=2))
    forward_aa=[]
    for pair in forward_aa_list:
        forward_aa+=["".join(pair)]
    reverse_aa_list=itertools.product(aromatic_amino_acids_reverse, repeat=2)
    reverse_aa=[]
    for pair in reverse_aa_list:
        reverse_aa+=["".join(pair)]
    for i in range(minlengthreadingframe,maxlengthreadingframe):
        i=i*3
        for end in stopcodon1:
            for n in range (1,5):
                for start in startcodon1:
                    for forward in forward_aa:
                        listreadingframe=listreadingframe+[start+"[A-Z]{"+str(i)+"}"+"[A-Z]{"+str((4-n)*3)+"}"+str(forward)+"[A-Z]{"+str((n)*3)+"}"+str(end)]
    for i in range(minlengthreadingframe,maxlengthreadingframe):
        i=i*3
        for end in stopcodon2:
            for n in range (1,5):
                for start in startcodon2:
                    for reverse in reverse_aa:
                        listreadingframe=listreadingframe+[end+"[A-Z]{"+str((4-n)*3)+"}"+str(reverse)+"[A-Z]{"+str((n)*3)+"}"+"[A-Z]{"+str(i)+"}"+str(start)]
    return listreadingframe
listreadingframe=generate_necessary_motifs(listreadingframe)
#create reading frame patterns
for i in range(0,15):
    dalgarno1=dalgarno1+["AGGAGG[A-Z]{"+str(i)+"}"]
    dalgarno2=dalgarno2+["[A-Z]{"+str(i)+"}CCTCCT"]
ribosomalbindingsites1=ribosomalbindingsites1+dalgarno1
ribosomalbindingsites2=ribosomalbindingsites2+dalgarno2
for i in range(minlengthreadingframe,maxlengthreadingframe):
    i=i*3
    for end in stopcodon1:
        for n in range (4,8):
            for start in startcodon1:
                for rbs in ribosomalbindingsites1:

                    listreadingframe=listreadingframe+["(?<="+str(rbs)+"[A-Z]{"+str(n)+"}"+")"+start+"[A-Z]{"+str(i)+"}"+str(end)]
for i in range(minlengthreadingframe,maxlengthreadingframe):
    i=i*3
    for end in stopcodon2:
        for n in range (4,8):
            for start in startcodon2:
                for rbs in ribosomalbindingsites2:
                    newframe=end+"[A-Z]{"+str(i)+"}"+start+"(?=[A-Z]{"+str(n)+"}"+str(rbs)+")"
                    listreadingframe=listreadingframe+ [newframe]
for i in range(minlengthreadingframe-7,maxlengthreadingframe):
    i=i*3
    for end in commoncore:
        for start in startcodon1:
            listreadingframe=listreadingframe+[start+"[A-Z]{"+str(i)+"}"+str(end)]
for i in range(minlengthreadingframe-7,maxlengthreadingframe):
    i=i*3
    for end in commoncorerev:
        for start in startcodon2:
            newframe=end+"[A-Z]{"+str(i)+"}"+start
            listreadingframe=listreadingframe+ [newframe]
#safe file
with open(filenamereadingframes, 'w') as f:
    for s in listreadingframe:
        f.write(str(s) + '\n')