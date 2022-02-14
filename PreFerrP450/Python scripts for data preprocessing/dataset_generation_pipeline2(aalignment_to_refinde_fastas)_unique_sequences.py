#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 19 18:48:47 2022
splits alignments at specific positions and Efficiently remove all duplicates+ remove all empty files+keep only files present in all clusters+ only keep files that are also present in specifgic fasta
@author: friederike

"""
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os
import glob
import re
start=0
end=400
without_duplicates=False
clusterlist= range(1,50)
splitting_list=[["begin",start,92],["sbr1",93,192],["sbr2",193,275],["core",276,395],["end",396,end],["fes1",54,115],["fes2",302,401]]
fragments=["begin","sbr1","sbr2","core","end","fes1","fes2"]
input_alignments="Datasets with different SSN thresholds/neighbouring=4_threshold=20/Alignment/"
master_fasta="/home/friederike/Documents/Databases/proGenomes211012022/extracted ferredoxins from all datasets neighbouring=4/trial of dereplication/dereplicated_p450s_by_CD_HIT_95/1644253029.fas.1"
list_headers_master_fasta=[]
for record in SeqIO.parse(master_fasta, "fasta"):
            list_headers_master_fasta=list_headers_master_fasta+[record.id]
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
            if "Reference_P450" not in record.id:
                subsequence=str(record.seq)[start-1:end]
                fasta.append(SeqRecord(Seq(subsequence),id=record.id))
    else:
        for record in alignment:
            if "Reference_P450" not in record.id:
                subsequence=str(record.seq)[start-1:]
                fasta.append(SeqRecord(Seq(subsequence),id=record.id))

    SeqIO.write(fasta, file_name, 'fasta')
def remove_duplicates_in_alignment(alignment):
        sequences=[]
        alignment_without_duplicates=[]
        for record in alignment:
            if str(record.seq).replace("-","") not in sequences or record.id == "Reference_P450":
                sequences.append(str(record.seq).replace("-",""))
                alignment_without_duplicates.append(record)
        print (len(alignment_without_duplicates))
        return alignment_without_duplicates


for alignment_file in glob.glob(os.path.join(input_alignments, '* Alignment.fasta')):

    output_name=alignment_file[:-6]
    alignment = AlignIO.read(open(alignment_file), "fasta")
    if without_duplicates==True:
        alignment= remove_duplicates_in_alignment(alignment)
    for record in alignment:
        
        if record.id=="Reference_P450":
            print ("x")
            index_reference=indexing_reference(record)
            converted_splitting_list=convert_splitting_list(splitting_list,index_reference)
            for segment in converted_splitting_list:
                split_alignment_and_put_to_fasta(alignment,segment,output_name)

for cluster in clusterlist:
    list_filenames=[]
    list_output_files=[]
    list_id_lists=[]
    final_id_list=[]
    list_output=[]
    for index_fragment,fragment in enumerate(fragments):
        list_id_lists.append([])
        list_output.append([])
        list_output_files.append("Datasets with different SSN thresholds/neighbouring=4_threshold=20/Alignment//whole_dataset_neighbouring_4_threshold_20_deduplicated_"+str(cluster)+ " Alignment_consens_"+fragment+".fasta")
        for record in SeqIO.parse("Datasets with different SSN thresholds/neighbouring=4_threshold=20/Alignment/whole_dataset_neighbouring_4_threshold_20_deduplicated_"+str(cluster)+ " Alignment"+fragment+".fasta", "fasta"):
            if record.id not in list_id_lists and re.findall(re.compile("_2$"), record.id)==[] and str(record.seq).replace("-","") !="":
                list_id_lists[index_fragment].append(record.id)
              
    for id in list_id_lists[0]:
        if id in list_id_lists[1] and  id in list_id_lists[2] and  id in list_id_lists[3] and id in list_id_lists[4] and  id in list_id_lists[5] and  id in list_id_lists[6] and  id in list_headers_master_fasta:
            final_id_list.append(id)
    for index_fragment,fragment in enumerate(fragments):
            list_records=[]
            for record in SeqIO.parse("Datasets with different SSN thresholds/neighbouring=4_threshold=20/Alignment/whole_dataset_neighbouring_4_threshold_20_deduplicated_"+str(cluster)+ " Alignment"+fragment+".fasta", "fasta"):
                    if record.id in final_id_list and record.id not in list_records:

                        list_records.append(record.id)
                        list_output[index_fragment].append(SeqRecord(Seq(str(record.seq).replace("-","")),id=record.id))
    for index_file,filename_out in enumerate(list_output_files):
        SeqIO.write(list_output[index_file], filename_out, "fasta")


            