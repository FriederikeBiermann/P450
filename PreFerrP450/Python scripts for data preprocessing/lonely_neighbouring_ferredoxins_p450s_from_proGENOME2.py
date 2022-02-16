#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 15:34:28 2022

@author: friederike
#a tool to find all neighbouring P450/ferredoxins in fasta from ProGenome2 database
"""

import Bio
import re
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Entrez
import os
import glob

#fill in correct filenames for input/output files

path = '/home/friederike/Documents/Databases/proGenomes211012022/Raw data/'
filename_neighbouring_out="tablep450ferredoxin_neighbouring_progenome2_neighbouring=4.csv"
filename_lonely_out="tablep450ferredoxin_lonely_progenome2_neighbouring=4.csv"
output_ferredoxin_lonely="lonely_ferredoxins_progenome2_neighbouring=4.fasta"
output_p450_neighbouring="neighbouring_p450ferre_progenome2_neighbouring=4.fasta"
output_ferredoxin_neighbouring="neighbouring_ferredoxins_progenome2_neighbouring=4.fasta"
output_p450_lonely="lonely_p450ferre_progenome2_neighbouring=4.fasta"
gene_difference_neighbouring=4

#creates empty data frame
lonely_ferredoxins_p450s=pd.DataFrame(columns=["accession",'p450_id','p450_sequence','p450_description', 'ferredoxin_id','ferredoxin_sequence', 'ferredoxin_description'])
neighbouring_ferredoxins_p450s=pd.DataFrame(columns=["accession",'p450_id','p450_sequence','p450_description',  'ferredoxin_id','ferredoxin_sequence', 'ferredoxin_description'])
def match_ferredoxins_p450s (ferredoxins,p450s):
    try:
        global neighbouring_ferredoxins_p450s
        global lonely_ferredoxins_p450s
        for index_ferredoxin,ferredoxin in enumerate(ferredoxins):
            if index_ferredoxin != len(ferredoxins):
                for index_p450,p450 in enumerate(p450s):
                    if index_p450 != len(p450s):
                        if -gene_difference_neighbouring<ferredoxin[0]-p450[0]<gene_difference_neighbouring:
                            
                            new_row={"accession":ferredoxin[3], 'p450_id':p450[1], 'p450_sequence':p450[2],'p450_description':p450[4].description, 'ferredoxin_id':ferredoxin[1], 'ferredoxin_sequence':ferredoxin[2],'ferredoxin_description':ferredoxin[4].description}
                            neighbouring_ferredoxins_p450s = neighbouring_ferredoxins_p450s.append(new_row, ignore_index=True)
                            del p450s[index_p450]
                            del ferredoxins[index_ferredoxin]
                            list_neighbouring_p450.append(p450[4])
                            list_neighbouring_ferredox.append(ferredoxin[4])
                    
        if len(ferredoxins)==1 and len(p450s)>0:
            for ferredoxin in ferredoxins:
                for p450 in p450s:
                    new_row={"accession":ferredoxin[3],'p450_id':p450[1],'p450_sequence':p450[2],'p450_description':p450[4].description, 'ferredoxin_id':ferredoxin[1],'ferredoxin_sequence':ferredoxin[2],'ferredoxin_description':ferredoxin[4].description}
                    lonely_ferredoxins_p450s = lonely_ferredoxins_p450s.append(new_row, ignore_index=True)
                    list_lonely_p450.append(p450[4])
                    list_lonely_ferredox.append(ferredoxin[4])
    except: print (ferredoxins)

 
p450s=[]
ferredoxins=[]
list_lonely_p450=[]
list_lonely_ferredox=[]
list_neighbouring_p450=[]
list_neighbouring_ferredox=[]
old_sample_id=""
#process antismash file
for progenome_file in glob.glob(os.path.join(path, '*.fasta')):

        progenome_file_name=str(progenome_file.split("/")[-1])
        filename_neighbouring_out_n=str(progenome_file_name)+filename_neighbouring_out
        filename_lonely_out_n=str(progenome_file_name)+filename_lonely_out
        output_ferredoxin_lonely_n=str(progenome_file_name)+output_ferredoxin_lonely
        output_p450_neighbouring_n=str(progenome_file_name)+output_p450_neighbouring
        output_ferredoxin_neighbouring_n=str(progenome_file_name)+output_ferredoxin_neighbouring
        output_p450_lonely_n=str(progenome_file_name)+output_p450_lonely
        print(output_p450_lonely_n)
        #for each genome
        for progenome_record in SeqIO.parse(progenome_file, "fasta"):
            # seperate record identifiers by starting a new list, every time a new record identifier appears
            record_identifiers=progenome_record.id.split(".")
            if record_identifiers[1]==old_sample_id:
                    protein_index+=1
                    if 'ferredoxin"' in progenome_record.description or 'Ferredoxin,' in progenome_record.description or '4Fe-4S ferredoxin' in progenome_record.description :
                            ferredoxins.append([protein_index,record_identifiers[2],progenome_record.seq,record_identifiers[1],progenome_record])  
                    if " p450" in progenome_record.description or " P450" in progenome_record.description and 'ferredoxin"' not in progenome_record.description and 'Ferredoxin,' not in progenome_record.description:
                            p450s.append([protein_index,record_identifiers[2],progenome_record.seq,record_identifiers[1],progenome_record])

            else:   
                    match_ferredoxins_p450s(ferredoxins,p450s)
                    old_sample_id=record_identifiers[1]
                    protein_index=1
                    p450s=[]
                    ferredoxins=[]
                    if 'ferredoxin"' in progenome_record.description or 'Ferredoxin,' in progenome_record.description or '4Fe-4S ferredoxin' in progenome_record.description :
                            ferredoxins.append([protein_index,record_identifiers[2],progenome_record.seq,record_identifiers[1],progenome_record])  
                    if " p450" in progenome_record.description or " P450" in progenome_record.description and 'ferredoxin"' not in progenome_record.description and 'Ferredoxin,' not in progenome_record.description:
                            p450s.append([protein_index,record_identifiers[2],progenome_record.seq,record_identifiers[1],progenome_record])

        #safe files
        SeqIO.write(list_lonely_p450, output_p450_lonely_n, 'fasta')
        SeqIO.write(list_lonely_ferredox,output_ferredoxin_lonely_n, 'fasta')
        SeqIO.write(list_neighbouring_p450, output_p450_neighbouring_n, 'fasta')
        SeqIO.write(list_neighbouring_ferredox,output_ferredoxin_neighbouring_n, 'fasta')
        lonely_ferredoxins_p450s.to_csv(filename_lonely_out_n, index=False) 
        neighbouring_ferredoxins_p450s.to_csv(filename_neighbouring_out_n, index=False) 


