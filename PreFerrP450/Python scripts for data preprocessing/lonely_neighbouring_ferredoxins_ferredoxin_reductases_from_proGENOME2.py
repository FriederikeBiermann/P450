#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 15:34:28 2022

@author: friederike
#a tool to find all neighbouring ferredoxin_NAD_reductase/ferredoxins in fasta from ProGenome2 database
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

path = '/home/friederike/Documents/Databases/proGenomes211012022/Raw data'
filename_neighbouring_out="tableferredoxin_NAD_reductaseferredoxin_neighbouring_progenome2_neighbouring=4.csv"
filename_lonely_out="tableferredoxin_NAD_reductaseferredoxin_lonely_progenome2_neighbouring=4.csv"
output_ferredoxin_lonely="lonely_ferredoxins_progenome2_neighbouring=4.fasta"
output_ferredoxin_NAD_reductase_neighbouring="neighbouring_ferredoxin_NAD_reductaseferre_progenome2_neighbouring=4.fasta"
output_ferredoxin_neighbouring="neighbouring_ferredoxins_progenome2_neighbouring=4.fasta"
output_ferredoxin_NAD_reductase_lonely="lonely_ferredoxin--NAD_reductaseferre_progenome2_neighbouring=4.fasta"
gene_difference_neighbouring=10

#creates empty data frame
lonely_ferredoxins_ferredoxin_NAD_reductases=pd.DataFrame(columns=["accession",'ferredoxin_NAD_reductase_id','ferredoxin_NAD_reductase_sequence','ferredoxin_NAD_reductase_description', 'ferredoxin_id','ferredoxin_sequence', 'ferredoxin_description'])
neighbouring_ferredoxins_ferredoxin_NAD_reductases=pd.DataFrame(columns=["accession",'ferredoxin_NAD_reductase_id','ferredoxin_NAD_reductase_sequence','ferredoxin_NAD_reductase_description',  'ferredoxin_id','ferredoxin_sequence', 'ferredoxin_description'])
def match_ferredoxins_ferredoxin_NAD_reductases (ferredoxins,ferredoxin_NAD_reductases):
   # try:
        global neighbouring_ferredoxins_ferredoxin_NAD_reductases
        global lonely_ferredoxins_ferredoxin_NAD_reductases
        for index_ferredoxin,ferredoxin in enumerate(ferredoxins):
            if index_ferredoxin != len(ferredoxins):
                for index_ferredoxin_NAD_reductase,ferredoxin_NAD_reductase in enumerate(ferredoxin_NAD_reductases):
                    if index_ferredoxin_NAD_reductase != len(ferredoxin_NAD_reductases):
                        if -gene_difference_neighbouring<ferredoxin[0]-ferredoxin_NAD_reductase[0]<gene_difference_neighbouring and ferredoxin[0]!=ferredoxin_NAD_reductase[0]:
              
                            new_row={"accession":ferredoxin[3], 'ferredoxin_NAD_reductase_id':ferredoxin_NAD_reductase[1], 'ferredoxin_NAD_reductase sequence':ferredoxin_NAD_reductase[2],'ferredoxin_NAD_reductase_description':ferredoxin_NAD_reductase[4].description, 'ferredoxin_id':ferredoxin[1], 'ferredoxin_sequence':ferredoxin[2],'ferredoxin_description':ferredoxin[4].description}
                            neighbouring_ferredoxins_ferredoxin_NAD_reductases = neighbouring_ferredoxins_ferredoxin_NAD_reductases.append(new_row, ignore_index=True)
                            del ferredoxin_NAD_reductases[index_ferredoxin_NAD_reductase]
                            del ferredoxins[index_ferredoxin]
                            list_neighbouring_ferredoxin_NAD_reductase.append(ferredoxin_NAD_reductase[4])
                            list_neighbouring_ferredox.append(ferredoxin[4])
                    
        if len(ferredoxins)>0 and len(ferredoxin_NAD_reductases)==1:
            for ferredoxin in ferredoxins:
                for ferredoxin_NAD_reductase in ferredoxin_NAD_reductases:
                    new_row={"accession":ferredoxin[3],'ferredoxin_NAD_reductase_id':ferredoxin_NAD_reductase[1],'ferredoxin_NAD_reductase_sequence':ferredoxin_NAD_reductase[2],'ferredoxin_NAD_reductase_description':ferredoxin_NAD_reductase[4].description, 'ferredoxin_id':ferredoxin[1],'ferredoxin_sequence':ferredoxin[2],'ferredoxin_description':ferredoxin[4].description}
                    lonely_ferredoxins_ferredoxin_NAD_reductases = lonely_ferredoxins_ferredoxin_NAD_reductases.append(new_row, ignore_index=True)
                    list_lonely_ferredoxin_NAD_reductase.append(ferredoxin_NAD_reductase[4])
                    list_lonely_ferredox.append(ferredoxin[4])
    #except: print (ferredoxins)


ferredoxin_NAD_reductases=[]
ferredoxins=[]
list_lonely_ferredoxin_NAD_reductase=[]
list_lonely_ferredox=[]
list_neighbouring_ferredoxin_NAD_reductase=[]
list_neighbouring_ferredox=[]
old_sample_id=""
#process antismash file
for progenome_file in glob.glob(os.path.join(path, '*.fasta')):

        progenome_file_name=str(progenome_file.split("/")[-1])
        filename_neighbouring_out_n=str(progenome_file_name)+filename_neighbouring_out
        filename_lonely_out_n=str(progenome_file_name)+filename_lonely_out
        output_ferredoxin_lonely_n=str(progenome_file_name)+output_ferredoxin_lonely
        output_ferredoxin_NAD_reductase_neighbouring_n=str(progenome_file_name)+output_ferredoxin_NAD_reductase_neighbouring
        output_ferredoxin_neighbouring_n=str(progenome_file_name)+output_ferredoxin_neighbouring
        output_ferredoxin_NAD_reductase_lonely_n=str(progenome_file_name)+output_ferredoxin_NAD_reductase_lonely
        print(output_ferredoxin_NAD_reductase_lonely_n)
        #for each genome
        for progenome_record in SeqIO.parse(progenome_file, "fasta"):
            # seperate record identifiers by starting a new list, every time a new record identifier appears
            record_identifiers=progenome_record.id.split(".")

            if record_identifiers[0]==old_sample_id:
                    protein_index+=1
                    if 'ferredoxin' in progenome_record.description or 'Ferredoxin' in progenome_record.description or '4Fe-4S ferredoxin' in progenome_record.description :
                        if "ferredoxin--NAD(+) reductase" not in progenome_record.description and "Ferredoxin-NADP reductase" not in progenome_record.description  and "ferredoxin--NADP(+) reductase" not in progenome_record.description and "subunit" not in progenome_record.description :
                            
                            ferredoxins.append([protein_index,record_identifiers[1],progenome_record.seq,record_identifiers[1],progenome_record])  
                    if "ferredoxin--NAD(+) reductase" in progenome_record.description or "Ferredoxin-NADP reductase" in progenome_record.description  or "ferredoxin--NADP(+) reductase" in progenome_record.description:
                         
                            ferredoxin_NAD_reductases.append([protein_index,record_identifiers[1],progenome_record.seq,record_identifiers[1],progenome_record])

            else:   
                    match_ferredoxins_ferredoxin_NAD_reductases(ferredoxins,ferredoxin_NAD_reductases)
                    old_sample_id=record_identifiers[0]
                    protein_index=1
                    ferredoxin_NAD_reductases=[]
                    ferredoxins=[]
                    if 'ferredoxin' in progenome_record.description or 'Ferredoxin' in progenome_record.description or '4Fe-4S ferredoxin' in progenome_record.description :
                        if "ferredoxin--NAD(+) reductase" not in progenome_record.description and "Ferredoxin-NADP reductase" not in progenome_record.description  and "ferredoxin--NADP(+) reductase" not in progenome_record.description and "subunit" not in progenome_record.description :
                            ferredoxins.append([protein_index,record_identifiers[1],progenome_record.seq,record_identifiers[0],progenome_record])  
                    if "ferredoxin--NAD(+) reductase" in progenome_record.description or "Ferredoxin-NADP reductase" in progenome_record.description  or "ferredoxin--NADP(+) reductase" in progenome_record.description:

                        ferredoxin_NAD_reductases.append([protein_index,record_identifiers[2],progenome_record.seq,record_identifiers[0],progenome_record])

        #safe files
        SeqIO.write(list_lonely_ferredoxin_NAD_reductase, output_ferredoxin_NAD_reductase_lonely_n, 'fasta')
        SeqIO.write(list_lonely_ferredox,output_ferredoxin_lonely_n, 'fasta')
        SeqIO.write(list_neighbouring_ferredoxin_NAD_reductase, output_ferredoxin_NAD_reductase_neighbouring_n, 'fasta')
        SeqIO.write(list_neighbouring_ferredox,output_ferredoxin_neighbouring_n, 'fasta')
        lonely_ferredoxins_ferredoxin_NAD_reductases.to_csv(filename_lonely_out_n, index=False) 
        neighbouring_ferredoxins_ferredoxin_NAD_reductases.to_csv(filename_neighbouring_out_n, index=False) 


