#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 08:50:35 2022

@author: friederike
"""
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 08:50:35 2022

@author: friederike
"""


import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import itertools
import math
import numpy as np
#choose clusters from the sequence similarity network to use as classes of ferredoxins
clusterlist =(1,2,3,4,5,6,7,8,9,10)
# fill in names of files here!
foldernameoutput="Output"
foldernameinput="whole_dataset_threshold_30/Fragments"
fragments=("begin","sbr1","sbr2","core","end","fes1","fes2")
filenamescomplete=[]
include_parts=False
for cluster in clusterlist:
    filenames_per_cluster=[]  
    for fragment in fragments:
        output=foldernameinput+"/whole_dataset_threshold_30_withreference_"+str(cluster)+" Alignment_consens_"+fragment+".fasta"
        filenames_per_cluster.append(output)
    filenamescomplete.append(filenames_per_cluster)
filename_permutations=foldernameoutput+"/permutations.txt"
fragments=("begin","sbr1","sbr2","core","end","fes1","fes2")
with open(filename_permutations, 'r') as file:
    permutations = [line.rstrip('\n') for line in file]

def find_between( s, first, last ):
    try:
        start = s.index( first ) + len( first )
        end = s.index( last, start )
        return s[start:end]
    except ValueError:
        return ""
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
def calculate_charge(sequence):
    AACharge = {"C":-.045,"D":-.999,  "E":-.998,"H":.091,"K":1,"R":1,"Y":-.001}
    charge = -0.002
    seqstr=str(sequence)
    seqlist=list(seqstr)
    for aa in seqlist:
        if aa in AACharge:
            charge += AACharge[aa]
    return charge

def merge_two_dicts(x, y):
    z = x.copy()   # start with x's keys and values
    z.update(y)    # modifies z with y's keys and values & returns None
    return z        

# creates tables
tablep450n=pd.DataFrame()
for counter,listfiles in enumerate(filenamescomplete):
        table=pd.DataFrame()
        table['target']=counter+1

        for counter2, file in enumerate(listfiles):
            print (counter, file)
            
            
            tempframe=pd.DataFrame()
 
            fragment=fragments[counter2]
            easy_sequences=[]
            aa_sequences=[]
            for seq_record in SeqIO.parse(file, "fasta"):

                seq_record.seq=Seq(str(seq_record.seq).replace('-', ''))
                easy_sequences.append(easysequence(str(seq_record.seq)))
                aa_sequences.append(seq_record.seq) 

            for index_sequence,easy_sequence in  enumerate(easy_sequences):       
                new_row={}
                #add "motif" features
                for motiv in permutations:
                    #defines column of table where to add
                    name=motiv+fragment
                    new_row =merge_two_dicts(new_row,{name:easy_sequence.count(motiv)})
                #add general features
                acidic=fragment+"acidic"
                new_row =merge_two_dicts(new_row,{acidic:(easy_sequence.count("a")/len(easy_sequence)+1)})
                acidic_absolute=fragment+"acidic absolute"
                new_row =merge_two_dicts(new_row,{acidic_absolute:(easy_sequence.count("a"))})
                charge=fragment+"charge"
                new_row =merge_two_dicts(new_row,{charge:calculate_charge(aa_sequences[index_sequence])})
                basic=fragment+"basic"
                basic_absolute=fragment+"basic absolute"
                new_row =merge_two_dicts(new_row,{basic:(easy_sequence.count("b")/len(easy_sequence)+1)})
                new_row =merge_two_dicts(new_row,{basic_absolute:(easy_sequence.count("b"))})
                if include_parts==True:
                    if fragment in fragments:
                        how_many_parts=3
                        each_part_length = math.ceil(len(easy_sequence)/how_many_parts)
                        iterator = [iter(easy_sequence)] * each_part_length
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
                tempframe = tempframe.append(new_row, ignore_index=True)
            table=pd.concat((table,tempframe),axis=1)
        # if counter==0:
        #     table['target']=1
        # else:

        tablep450n=pd.concat((tablep450n,table),axis=0)
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

tablep450n['complete charge']=tablep450n[chargerows].sum(axis=1)
tablep450n['mean acidic']=tablep450n[acidicrows].mean(axis=1)  
tablep450n['mean basic']=tablep450n[basicrows].mean(axis=1)  
tablep450n['absolute acidic']=tablep450n[absacidicrows].sum(axis=1)  
tablep450n['absolute basic']=tablep450n[absbasicrows].sum(axis=1)
tablep450n = tablep450n.replace('nan', np.nan).fillna(0)
tablep450n.to_csv("whole_dataset_threshold_30/whole_dataset_threshold_30_neighbouring_feature_matrix.csv", index=False) 
