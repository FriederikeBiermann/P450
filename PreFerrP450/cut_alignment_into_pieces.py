#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 18:48:47 2022
splits alignments at specific positions
@author: friederike

"""
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import glob
start=0
end=400
splitting_list=[["begin",start,92],["sbr1",93,192],["sbr2",193,275],["core",276,395],["end",396,end],["fes1",54,115],["fes2",302,401]]
input_alignments="/home/friederike/Documents/Coding/P450/PreFerrP450/whole_dataset_threshold_30/Alignments/"

def indexing_reference(record):
    list_reference=list(str(record.seq))

    index_aa=0
    index_mapping=[]
    for index,AA in enumerate(list_reference):
        if AA !="-":
            index_aa+=1
            index_mapping.append([index_aa,index])

    return (index_mapping)
def convert_splitting_list(splitting_list,index_reference):
    converted_splitting_list=[]
    for segment in splitting_list:
        converted_splitting_list.append([segment[0],index_reference[segment[1]][1],index_reference[segment[2]-1][1]])
    return converted_splitting_list
def split_alignment_and_put_to_fasta(alignment,segment,output_name):
    file_name=output_name+segment[0]+".fasta"
    start=segment[1]
    end=segment[2]
    fasta=[]
    if segment[0]=="begin":
        start=1
    if segment[0]!="end":
        for record in alignment:
            subsequence=str(record.seq)[start-1:end]
            fasta.append(SeqRecord(Seq(subsequence),id=record.id))
    else:
        for record in alignment:
            subsequence=str(record.seq)[start-1:]
            fasta.append(SeqRecord(Seq(subsequence),id=record.id))
    SeqIO.write(fasta, file_name, 'fasta')
for alignment_file in glob.glob(os.path.join(input_alignments, '* Alignment.fasta')):
    print (alignment_file)
    output_name=alignment_file[:-6]
    alignment = AlignIO.read(open(alignment_file), "fasta")
    for record in alignment:
        
        if record.id=="Reference_P450":
            print (record.seq)
            index_reference=indexing_reference(record)
            print(index_reference)
            converted_splitting_list=convert_splitting_list(splitting_list,index_reference)
            for segment in converted_splitting_list:
                split_alignment_and_put_to_fasta(alignment,segment,output_name)
            