#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 07:52:12 2022
Efficiently remove all duplicates+ remove all empty files+keep only files present in all clusters
@author: friederike
"""

import Bio
import re
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import os

fragments=["begin","sbr1","sbr2","core","end","fes1","fes2"]
for cluster in range(1,11):
    list_filenames=[]
    list_output_files=[]
    list_id_lists=[]
    final_id_list=[]
    list_output=[]
    for index_fragment,fragment in enumerate(fragments):
        list_id_lists.append([])
        list_output.append([])
        list_output_files.append("whole_dataset_threshold_30/Fragments/whole_dataset_threshold_30_withreference_"+str(cluster)+ " Alignment_consens_"+fragment+".fasta")
        for record in SeqIO.parse("whole_dataset_threshold_30/Fragments/whole_dataset_threshold_30_withreference_"+str(cluster)+ " Alignment"+fragment+".fasta", "fasta"):
            if record.id not in list_id_lists and re.findall(re.compile("_2$"), record.id)==[] and str(record.seq).replace("-","") !="":
                list_id_lists[index_fragment].append(record.id)
                print (record.id, cluster, "whole_dataset_threshold_30/Fragments/whole_dataset_threshold_30_withreference_"+str(cluster)+ " Alignment"+fragment+".fasta")
    for id in list_id_lists[0]:
        if id in list_id_lists[1] and  id in list_id_lists[2] and  id in list_id_lists[3] and id in list_id_lists[4] and  id in list_id_lists[5] and  id in list_id_lists[6]:
            final_id_list.append(id)
    for index_fragment,fragment in enumerate(fragments):
            list_records=[]
            for record in SeqIO.parse("whole_dataset_threshold_30/Fragments/whole_dataset_threshold_30_withreference_"+str(cluster)+ " Alignment"+fragment+".fasta", "fasta"):
                    if record.id in final_id_list and record.id not in list_records:

                        list_records.append(record.id)
                        list_output[index_fragment].append(SeqRecord(Seq(str(record.seq).replace("-","")),id=record.id))
    for index_file,filename_out in enumerate(list_output_files):
        SeqIO.write(list_output[index_file], filename_out, "fasta")
        
    