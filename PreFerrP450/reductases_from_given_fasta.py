#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  2 17:00:41 2022
get reductases assoicated with given ferredoxins
@author: friederike
"""

import re
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Entrez
from Bio.SeqRecord import SeqRecord

filename_input="Datasets with different SSN thresholds/neighbouring=4_threshold=20/ferredoxins/whole_dataset_neighbouring_4_threshold_20_deduplicated_ferredoxin_5.fasta"
filename_output=filename_input+"info.csv"
outputferre=filename_input+"ferre.fasta"
outputreductase=filename_input+"reductase.fasta"
tableofproteins=pd.DataFrame(columns=['Protein Reference','Acessiongenome','Species',"Definition","proteins","blasts","products"])
Entrez.email = "friederike@biermann-erfurt.de"
listfasta=[]
liststrains=[]

def find_between( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""
def find_infrontof( s, last ):
    try:
        start = 0
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""
def find_after( s, first ):
    try:
        start =s.index( first ) + len( first )
        end = len(s)
        return s[start:end]
    except ValueError:
        return ""

listreductase=[]
listferredoxins=[]
countsafe=0

#actually finding in fastas
for record in SeqIO.parse(filename_input, "fasta"):
    
    #try:
            countsafe=countsafe+1
            if "ref" in record.id:
                gi=find_between(record.id, "ref|","|")
            if "gb" in record.id:
                gi=find_between(record.id, "gb|","|") 
            if "sp" in record.id:
                pi=find_between(record.id, "sp|","|") 
                handle = Entrez.elink(dbFrom = "protein", db = "protein",id=pi)
                recordtest = str(Entrez.read(handle))
                handle.close()
                for match in re.finditer("(?<=Id': ')[0-9]+", recordtest):
                            gi=match.group()
            if "ferredoxin" not in record.description:
                gi=record.id
            else: 
                gi=find_between(record.description,'protein_id="','" /')
            print(gi)
      
            ref=str(record.seq)
            # fetch all genomes associated with protein
            handle = Entrez.elink(dbFrom = "protein", db = "nucleotide",id=gi)
            record3 = str(Entrez.read(handle))
            handle.close()
            listgenes=re.findall("(?<=Id': ')[0-9]+",record3)
            listgenes=listgenes[:int(len(listgenes)/2)]
            for name1 in listgenes[:2]:
                dnahandle1 = Entrez.efetch(db="nucleotide", id=name1, rettype='gbwithparts',retmode="xml")
                record4 = Entrez.read(dnahandle1)
                dnahandle1.close()
                definition=str(record4[0]["GBSeq_definition"])
                species = record4[0]["GBSeq_organism"]
                listproteins=[]
                listsequences=[]
                listproducts=[]
                # find protein of interrest in genome and search neighbouring genes for "reductase"
                for match in re.finditer(ref, str(record4)):
                    proteinspan1=match.span()[0]
                    proteinspan2=int((proteinspan1))-50000
                    proteinspan3=int((proteinspan1))+52000
                    intervalrecord=str(record4)[proteinspan2:proteinspan3]                   
                    if "reductase" in intervalrecord or "Reductase" in intervalrecord:
                        for match in re.finditer("(?<={'GBQualifier_name': 'product', 'GBQualifier_value': ')...................", intervalrecord):
                            if "reductase" in match.group() or "Reductase" in match.group():
                                interval2= intervalrecord[match.span()[0]-1000:match.span()[0]+1000]
                                listproducts.append(match.group())
                                for match in re.finditer("(?<={'GBQualifier_name': 'translation', 'GBQualifier_value': ')\w+", interval2):
                                    listsequences.append(match.group())
                                    translation=match.group()
                                for match in re.finditer("(?<={'GBQualifier_name': 'locus_tag', 'GBQualifier_value': ')\w+", interval2):
                                    listproteins.append(match.group())
                                    locus=match.group()
                                listreductase.append(SeqRecord(Seq(translation),id=locus))
                                listferredoxins.append(record)
                        res = [] 
                        [res.append(x) for x in listproteins if x not in res]
                        
                        #safe results
                        new_row={'Protein Reference':gi,'Acessiongenome':name1,'Species':species,"Definition":definition,"proteins":res,"blasts":listsequences,"products":listproducts}
                        tableofproteins = tableofproteins.append(new_row, ignore_index=True)
                        
                        #safe every 100 entries
                        if countsafe% 10 == 0:
                            tableofproteins.to_csv(filename_output, index=False) 
                            SeqIO.write(listferredoxins, outputferre, 'fasta')
                            SeqIO.write(listreductase, outputreductase, 'fasta')
    #except: print ("3")
      
#final saving of files
tableofproteins.to_csv(filename_output, index=False) 
SeqIO.write(listferredoxins, outputferre, 'fasta')
SeqIO.write(listreductase, outputreductase, 'fasta')