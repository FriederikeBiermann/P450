#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan  3 18:24:41 2022
GUI
@author: friederike
"""


import pandas as pd
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
from tkinter import *
import math
import itertools
import pickle
from tkinter import messagebox
fastas_aligned_before=False

include_charge_features=True
include_sequence_motifs=True
include_parts=True
#fill in filenames here!
foldernameoutput=""
filename_permutations=foldernameoutput+"permutations.txt"
alignmentfa=("MSAVALPRVSGGHDEHGHLEEFRTDPIGLMQRVRDECGDVGTFQLAGKQVVLLSGSHANEFFFRAGDDDLDQAKAYPFMTPIFGEGVVFDASPERRKEMLHNAALRGEQMKGHAATIEDQVRRMIADWGEAGEIDLLDFFAELTIYTSSACLIGKKFRDQLDGRFAKLYHELERGTDPLAYVDPYLPIESLRRRDEARNGLVALVADIMNGRIANPPTDKSDRDMLDVLIAVKAETGTPRFSADEITGMFISMMFAGHHTSSGTASWTLIELMRHRDAYAAVIDELDELYGDGRSVSFHLRQIPQLENVLKETLRLHPPLIILMRVAKGEFEVQGHRIHEGDLVAASPAISNRIPEDFPDPHDFVPARYEQPRQEDLLNRWTWIPFGAGRHRCVGAAFAIMQIKAIFSVLLREYEFEMAQPPESYRNDHSKMVVQLAQPACVRYRRRTGV")
used_classifier=foldernameoutput+'extratreeclassifier.sav'

start=0
end=400
splitting_list=[["begin",start,92],["sbr1",93,192],["sbr2",193,275],["core",276,395],["end",396,end],["fes1",54,115],["fes2",302,401]]
fragments=["begin","sbr1","sbr2","core","end","fes1","fes2"]
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
    if fastas_aligned_before==False:
        alignment=[alignment]
    seqRecord_list_per_fragment=np.array([["",""]])
    if fragment[0]=="begin":
        start=1
    if fragment[0]!="end":
        for record in alignment:
            subsequence=str(record.seq)[start-1:end]
     
            seqRecord_list_per_fragment=np.append(seqRecord_list_per_fragment,[[record.id,subsequence]],axis=0)
    else:
        for record in alignment:
            subsequence=str(record.seq)[start-1:]
            seqRecord_list_per_fragment=np.append(seqRecord_list_per_fragment,[[record.id,subsequence]],axis=0)
    seqRecord_list_per_fragment=np.delete(seqRecord_list_per_fragment, 0, axis=0)
    return seqRecord_list_per_fragment
def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z       
def fragment_alignment(alignment,splitting_list):

    fragment_matrix=pd.DataFrame()
    if fastas_aligned_before==False:
     
        seqa=alignment[0]
        seqb=alignment[1]
        index_reference=indexing_reference(SeqRecord(Seq(seqa),id=seqa))

        converted_splitting_list=convert_splitting_list(splitting_list,index_reference)
        for fragment in converted_splitting_list:
                name_fragment=fragment[0]
                seqRecord_list_per_fragment=split_alignment(SeqRecord(Seq(seqb),id=seqb),fragment)

                fragment_matrix[name_fragment]=seqRecord_list_per_fragment[:,1]
                fragment_matrix.set_index(pd.Index(seqRecord_list_per_fragment[:,0]))
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
                    fragment_matrix[name_fragment]=seqRecord_list_per_fragment[:,1]
                    fragment_matrix.set_index(pd.Index(seqRecord_list_per_fragment[:,0]))
    return fragment_matrix
def featurize(fragment_matrix, permutations, fragments, include_charge_features):
    feature_matrix=pd.DataFrame()
    for index, row in fragment_matrix.iterrows(): 
        new_row={}
        for fragment in fragments:
            sequence_fragment=row[fragment]
            easysequence_fragment=easysequence(sequence_fragment)
            if include_sequence_motifs==True:
                for motif in permutations:
                    name_column=motif+fragment
                    new_row =merge_two_dicts(new_row,{name_column:easysequence_fragment.count(motif)})
            if include_charge_features==True:
                new_row=append_charge_features(new_row,fragment,easysequence_fragment,sequence_fragment)
        feature_matrix=feature_matrix.append(new_row, ignore_index=True)
    if include_charge_features==True:
        feature_matrix=sum_charge_features(feature_matrix)
    print (feature_matrix)
    return feature_matrix

def append_charge_features(new_row,fragment,easysequence_fragment,sequence_fragment):
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

    if include_parts==True:
                if fragment in fragments:
                    how_many_parts=3
                    each_part_length = math.ceil(len(easysequence_fragment)/how_many_parts)
                    iterator = [iter(easysequence_fragment)] * each_part_length
                    parts = list(itertools.zip_longest(*iterator, fillvalue='X'))
                    for i in range (0,len(parts)):
                            part=fragment+"acidic"+ str(i)
                            par=parts[i]
                            new_row =merge_two_dicts(new_row,{part:(par.count("a")/(len(par)+1))})
                            part=fragment+"basic"+ str(i)
                            par=parts[i]
                            new_row =merge_two_dicts(new_row,{part:(par.count("b")/(len(par)+1))})
                    if len(parts)<3:
                        for dif in range ((len(parts)+1),4):
                            part=fragment+"acidic"+ str(dif)
                            new_row =merge_two_dicts(new_row,{part:0})
                            part=fragment+"basic"+ str(dif)
                            new_row =merge_two_dicts(new_row,{part:0})
    
    return new_row
def sum_charge_features(feature_matrix):
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
#getting Data for p450 associated with terpenes
def check_for_AA_sequence(ref_seq):
    if "M" in ref_seq or "S" in ref_seq or "m" in ref_seq or "s" in ref_seq:
        print ("")
    else:
        print (f'You might have inserted a DNA sequence!')
        return messagebox.showinfo('Warning', f'You might have inserted a DNA sequence!')
def ferredoxin_calculation (seq_record):
    seq_record =input_box.get("1.0", "end")
    check_for_AA_sequence(seq_record)
    seq_record=SeqRecord(Seq(seq_record),id="sequence of interest")
    if fastas_aligned_before==False:
        fragment_matrix=pd.DataFrame()
        seq_record_ids=[]
        fewgaps = lambda x, y: -20 - y
        specificgaps = lambda x, y: (-2 - y)
        alignment = pairwise2.align.globalmc(alignmentfa, seq_record.seq, 1, -1, fewgaps, specificgaps)
        fragment_matrix_for_record=fragment_alignment (alignment[0],splitting_list)
        fragment_matrix=fragment_matrix.append(fragment_matrix_for_record, ignore_index = True)
        print (fragment_matrix)
        seq_record_ids=seq_record_ids+[seq_record.id]
    feature_matrix=featurize(fragment_matrix, permutations, fragments, include_charge_features)
    print (feature_matrix)
    classifier = pickle.load(open(used_classifier, 'rb'))
    predicted_ferredoxin=classifier.predict(feature_matrix)
    y_score = classifier.predict_proba(feature_matrix)
    print (y_score)
    
    return messagebox.showinfo('Result', f'The matching family of ferredoxins is probably {predicted_ferredoxin}, the scores are {y_score}')

gui = Tk()
gui.title('PreFerrP450')
gui.geometry('500x500')

Label(gui, text='Enter amino acid sequence of CYP P450 & hit "Enter" Key').pack(pady=20)
input_box = Text(gui,height=2000, width= 450, yscrollcommand=TRUE )
input_box.bind('<Return>',ferredoxin_calculation)
input_box.pack()

gui.mainloop()
