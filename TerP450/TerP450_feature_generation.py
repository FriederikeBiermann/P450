#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 08:26:52 2022

@author: friederike
"""

import Bio
import pandas as pd
from Bio import SeqIO
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
import re

fastas_aligned_before=False

include_charge_features=False
#fill in filenames here!
foldernameoutput="Output/Terpene"
name_450terpenes="Fasta/non_terpene_p450.faa" #insert either name of fasta or alignment in fasta format (alignment must include "Reference_P450")
name_450nonterpenes="Fasta/terpene_p450.faa" #insert either name of fasta or alignment in fasta format(alignment must include "Reference_P450")
filename_permutations=foldernameoutput+"/permutations.txt"
alignmentfa=(SeqRecord(Seq("MSAVALPRVSGGHDEHGHLEEFRTDPIGLMQRVRDECGDVGTFQLAGKQVVLLSGSHANEFFFRAGDDDLDQAKAYPFMTPIFGEGVVFDASPERRKEMLHNAALRGEQMKGHAATIEDQVRRMIADWGEAGEIDLLDFFAELTIYTSSACLIGKKFRDQLDGRFAKLYHELERGTDPLAYVDPYLPIESLRRRDEARNGLVALVADIMNGRIANPPTDKSDRDMLDVLIAVKAETGTPRFSADEITGMFISMMFAGHHTSSGTASWTLIELMRHRDAYAAVIDELDELYGDGRSVSFHLRQIPQLENVLKETLRLHPPLIILMRVAKGEFEVQGHRIHEGDLVAASPAISNRIPEDFPDPHDFVPARYEQPRQEDLLNRWTWIPFGAGRHRCVGAAFAIMQIKAIFSVLLREYEFEMAQPPESYRNDHSKMVVQLAQPACVRYRRRTGV"),id="Reference_P450"))
path_complete_feature_matrix=foldernameoutput+"/path_complete_feature_matrix.csv"
#S#14703] cytochrome P450 51 Cyp51 ([P#10800] [CYP51] cytochrome P450, family 51 )  from https://cyped.biocatnet.de/sequence/14703

start=0
end=400
splitting_list=[["begin",start,92],["sbr1",93,192],["sbr2",193,275],["core",276,395],["end",396,end],["fes1",54,115],["fes2",302,401]]
fragments=splitting_list[:][0]
with open(filename_permutations, 'r') as file:
    permutations = [line.rstrip('\n') for line in file]
def find_between( string, first, last ):
    try:
        start = string.index( first ) + len( first )
        end = string.index( last, start )
        return string[start:end]
    except ValueError:
        return ""
def calculate_charge(sequence):
    AACharge = {"C":-.045,"D":-.999,  "E":-.998,"H":.091,"K":1,"R":1,"Y":-.001}
    charge = -0.002
    seqstr=str(sequence)
    seqlist=list(seqstr)
    for aa in seqlist:
        if aa in AACharge:
            charge += AACharge[aa]
    return charge
def easysequence (sequence):
    #creates a string out of the sequence file, that only states if AA is acidic (a), basic (b), polar (p), neutral/unpolar (n),aromatic (r),Cystein (s) or a Prolin (t)
    seqstr=str(sequence)
    seqlist=list(seqstr)
    easylist=[]
    for i in seqlist:
        if i == 'E' or i== 'D':           
            easylist=easylist+['a']
        if i == 'K' or i=='R' or i=='H':      
            easylist=easylist+['b']
        if i == 'S' or i=='T' or i=='N' or i=='Q':
            easylist=easylist+['p']
        if i == 'F' or i=='Y' or i=='W':
            easylist=easylist+['r']
        if i == 'C':
            easylist=easylist+['s']
        if i == 'P':
            easylist=easylist+['t']
        if i == 'G' or i=='A' or i=='V' or i=='L' or i=='I' or i=='M':
            easylist=easylist+['n']
            
    seperator=''
    easysequence=seperator.join(easylist)
    return easysequence
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
    for fragment in splitting_list:
        converted_splitting_list.append([fragment[0],index_reference[fragment[1]][1],index_reference[fragment[2]-1][1]])
    return converted_splitting_list
def split_alignment(alignment,fragment):
    start=fragment[1]
    end=fragment[2]
    seqRecord_list_per_fragment=[]
    if fragment[0]=="begin":
        start=1
    if fragment[0]!="end":
        for record in alignment:
            subsequence=str(record.seq)[start-1:end]
            seqRecord_list_per_fragment.append([record.id,subsequence])
    else:
        for record in alignment:
            subsequence=str(record.seq)[start-1:]
            seqRecord_list_per_fragment.append([record.id,subsequence])
    return seqRecord_list_per_fragment
def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z       
def fragment_alignment(alignment,splitting_list):

    fragment_matrix=pd.DataFrame()
    if type(alignment)==class 'Bio.pairwise2.Alignment':
        print (x)
        seqa=find_between(str(alignments[0]),"seqA='","'")
        seqb=find_between(str(alignments[0]),"seqB='","'")
        index_reference=indexing_reference(seqa)
        converted_splitting_list=convert_splitting_list(splitting_list,index_reference)
        for fragment in converted_splitting_list:
                name_fragment=fragment[0]
                seqRecord_list_per_fragment=split_alignment(SeqRecord(Seq(seqb),id=find_between(str(alignments[0]),'seqB="ID:','\\'),fragment)
                fragment_matrix[name_fragment]=seqRecord_list_per_fragment[1]
                fragment_matrix.set_index(pd.Index(seqRecord_list_per_fragment[0]))
    else:
        for record in alignment:     
            print (record)
            if record.id=="Reference_P450":
                print (record.seq)
                index_reference=indexing_reference(record)
                print(index_reference)
                converted_splitting_list=convert_splitting_list(splitting_list,index_reference)
                for fragment in converted_splitting_list:
                    name_fragment=fragment[0]
                    seqRecord_list_per_fragment=split_alignment(alignment,fragment)
                    fragment_matrix[name_fragment]=seqRecord_list_per_fragment[1]
                    fragment_matrix.set_index(pd.Index(seqRecord_list_per_fragment[0]))
    return fragment_matrix
def featurize(fragment_matrix, permutations, fragments, include_charge_features):
    feature_matrix=pd.DataFrame()
    for index, row in fragment_matrix.iterrows(): 
        new_row={}
        for fragment in fragments:
            sequence_fragment=row[fragment]
            easysequence_fragment=easysequence(sequence_fragment)
            for motif in permutations:
                name_column=motif+fragment
                new_row =merge_two_dicts(new_row,{name_column:easysequence_fragment.count(motif)})
        if include_charge_features==True:
                acidic=fragment+"acidic"
                new_row =merge_two_dicts(new_row,{acidic:(easysequence_fragment.count("a")/len(easysequence_fragment)+1)})
                acidic_absolute=fragment+"acidic absolute"
                new_row =merge_two_dicts(new_row,{acidic_absolute:(easysequence_fragment.count("a"))})
                charge_name=fragment+"charge"
                new_row =merge_two_dicts(new_row,{charge_name:(calculate_charge(sequence_fragment))})
                basic=fragment+"basic"
                basic_absolute=fragment+"basic absolute"
                new_row =merge_two_dicts(new_row,{basic:(easysequence_fragment.count("b")/len(easysequence_fragment)+1)})
                new_row =merge_two_dicts(new_row,{basic_absolute:(easysequence_fragment.count("b"))})
        feature_matrix.append(new_row, ignore_index=True)
    feature_matrix.set_index(fragment_matrix.index)
    if include_charge_features==True:
        chargerows=[]
        acidicrows=[]
        basicrows=[]
        absacidicrows=[]
        absbasicrows=[]
        for fragment in fragments:
            chargerows.append(str(fragment)+"charge")
            acidicrows.append(str(fragment)+"acidic")
            basicrows.append(str(fragment)+"basic")
            absacidicrows.append(str(fragment)+"acidic absolute")
            absbasicrows.append(str(fragment)+"basic absolute")
        feature_matrix['complete charge']=feature_matrix[chargerows].sum(axis=1)
        feature_matrix['mean acidic']=feature_matrix[acidicrows].mean(axis=1)  
        feature_matrix['mean basic']=feature_matrix[basicrows].mean(axis=1)  
        feature_matrix['absolute acidic']=feature_matrix[absacidicrows].sum(axis=1)  
        feature_matrix['absolute basic']=feature_matrix[absbasicrows].sum(axis=1)
    return feature_matrix
complete_feature_matrix=pd.DataFrame()
for dataset in (name_450terpenes,name_450nonterpenes):
    if fastas_aligned_before==True:
        alignment = AlignIO.read(open(dataset), "fasta")
        fragment_matrix=fragment_alignment (alignment,splitting_list)
    if fastas_aligned_before==False:
        fragment_matrix=pd.DataFrame()
        seq_record_ids=[]
        for seq_record in SeqIO.parse(dataset, "fasta"):
             fewgaps = lambda x, y: -20 - y
             specificgaps = lambda x, y: (-2 - y)
             alignment = pairwise2.align.globalmc(alignmentfa,  seq_record, 1, -1, fewgaps, specificgaps,one_alignment_only=True)
             print (alignment[0], type(alignment[0]))
             fragment_matrix_for_record=fragment_alignment (alignment,splitting_list)
             fragment_matrix.append(fragment_matrix_for_record, ignore_index = True)
             seq_record_ids.append(seq_record.id)
        feature_matrix.set_index(pd.Index(seq_record_ids))     
    feature_matrix=featurize(fragment_matrix, permutations, fragments, include_charge_features)

    if dataset==name_450terpenes:
        feature_matrix["target"]=1
    if dataset==name_450nonterpenes: 
         feature_matrix["target"]=0
    complete_feature_matrix.append(fragment_matrix_for_record, ignore_index = True)
    
#modify table (drop accessionnumber, sequences)

complete_feature_matrix.to_csv(path_complete_feature_matrix, index=False)    