#!/usr/bin/env python
# coding: utf-8

# In[3]:


#!/usr/bin/env python
# coding: utf-8

# In[6]:


#!/usr/bin/env python
# coding: utf-8

# In[19]:


import subprocess
from Bio import SeqIO
import csv
import pickle

chrom_lens = {"chr1": 248956422, 
             "chr2": 242193529, 
             "chr3": 198295559,
             "chr4": 190214555,
             "chr5": 181538259,
             "chr6": 170805979,
             "chr7": 159345973,
             "chr8": 145138636,
             "chr9": 138394717,
             "chr10": 133797422,
             "chr11": 135086622,
             "chr12": 133275309,
             "chr13": 114364328,
             "chr14": 107043718,
             "chr15": 101991189,
             "chr16": 90338345,
             "chr17": 83257441,
             "chr18": 80373285,
             "chr19": 58617616,
             "chr20": 64444167,
             "chr21": 46709983,
             "chr22": 50818468,
             "chrX": 156040895,
             "chrY": 57227415}

# In[ ]:


# In[20]:
def tocsv(dictName):
    with open("overlap_pn_frags_mrna_seqs.tsv", "w") as outfile:
        writer = csv.writer(outfile, delimiter = "\t")
        for data in dictName:
            writer.writerow(data)
    return

def dna_to_mrna_fwd(input_fasta):
    trslt_fasta = ''
    l = list(input_fasta)  
    for i in range(len(input_fasta)):
        if(l[i]=='G'):
            l[i]='C'
  
        elif(l[i]=='C'):
            l[i]='G'
  
        elif (l[i] == 'T'):
            l[i] = 'A'
  
        elif (l[i] == 'A'):
            l[i] = 'U'
            
        elif(l[i]=='g'):
            l[i]='c'   
            
        elif(l[i]=='c'):
            l[i]='g'
  
        elif (l[i] == 't'):
            l[i] = 'a'
  
        elif (l[i] == 'a'):
            l[i] = 'u'
    for char in l:
        trslt_fasta += char
    return trslt_fasta

def dna_to_mrna_rev(input_fasta):
    trslt_fasta = ''
    l = list(input_fasta)  
    for i in range(len(input_fasta)-1, -1, -1):
        if(l[i]=='G'):
            l[i]='C'
  
        elif(l[i]=='C'):
            l[i]='G'
  
        elif (l[i] == 'T'):
            l[i] = 'A'
  
        elif (l[i] == 'A'):
            l[i] = 'U'
            
        elif(l[i]=='g'):
            l[i]='c'   
            
        elif(l[i]=='c'):
            l[i]='g'
  
        elif (l[i] == 't'):
            l[i] = 'a'
  
        elif (l[i] == 'a'):
            l[i] = 'u'
    for char in reversed(l):
        trslt_fasta += char
    return trslt_fasta

def pad_seq(g_seq, g_start, g_end):
    tholder = []
    g_start = int(g_start)
    g_end = int(g_end)
    if g_start < 350:
        g_start = 0
        g_end = g_end + 700-g_start
    elif g_end > chrom_lens[g_seq] - 350:
        g_end = chrom_lens[g_seq]
        g_start = g_start-700+chrom_lens[g_seq]-g_end
    else:
        g_start = g_start - 350
        g_end = g_end + 350
    tholder.append(g_seq)
    tholder.append(g_start)
    tholder.append(g_end)
    return tholder
# In[32]:


def file_to_mrna_pos(input_file):
    count =0
    with open(input_file) as f:
        for line in f:
            L = line.strip().split()
            if L[3] in frag_and_mrna_seq:
                count+=1
                continue
                
            else:
                pholder = pad_seq(L[0], L[1], L[2])

                given_seq = pholder[0]
                given_start = str(pholder[1])
                given_end = str(pholder[2])

                pos_strand = True
                if (str(L[5]) == "+"):
                    pos_strand = True
                elif (str(L[5]) == '-'):
                    pos_strand = False
                #ar3 = str(input_file) + '.fa'
                ar4 = '-seq=' + str(given_seq)
                ar5 = "-start=" + str(given_start)
                ar6 = "-end=" + str(given_end)

                ret_seq = subprocess.Popen(["./twoBitToFa", "-noMask", "hg38.2bit", "placeholder_file.fa", ar4, ar5, ar6])
                for seq_record in SeqIO.parse("placeholder_file.fa", "fasta"):
                    og_seq = str(seq_record.seq)
                    new_seq = ''
                    if (pos_strand == True):
                        new_seq = dna_to_mrna_fwd(og_seq)
                    else:
                        new_seq = dna_to_mrna_rev(og_seq)

                    frag_and_mrna_seq[L[3]] = new_seq 
                count+=1
                if (count%5000 ==0):
                    print(count)
    print(count)
    return 

def file_to_mrna_neg(input_file):
    count =0
    with open(input_file) as f:
        for line in f:
            L = line.strip().split()
            if L[3] in frag_and_mrna_seq:
                count+=1
                continue
                
            else:
                pholder = pad_seq(L[0], L[1], L[2])

                given_seq = pholder[0]
                given_start = str(pholder[1])
                given_end = str(pholder[2])

                pos_strand = True
                if (str(L[5]) == "+"):
                    pos_strand = True
                elif (str(L[5]) == '-'):
                    pos_strand = False
                #ar3 = str(input_file) + '.fa'
                ar4 = '-seq=' + str(given_seq)
                ar5 = "-start=" + str(given_start)
                ar6 = "-end=" + str(given_end)

                ret_seq = subprocess.Popen(["./twoBitToFa", "-noMask", "hg38.2bit", "placeholder_file.fa", ar4, ar5, ar6])
                for seq_record in SeqIO.parse("placeholder_file.fa", "fasta"):
                    og_seq = str(seq_record.seq)
                    new_seq = ''
                    if (pos_strand == True):
                        new_seq = dna_to_mrna_fwd(og_seq)
                    else:
                        new_seq = dna_to_mrna_rev(og_seq)

                    frag_and_mrna_seq[L[3]] = new_seq 
                count+=1
                if (count%5000 ==0):
                    print(count)
    print(count)
    return 
# In[26]:

frag_and_mrna_seq = {}
temp1 = file_to_mrna_pos("renamed_pos_overlaps.bed")
temp2 = file_to_mrna_neg("renamed_neg_overlaps.bed")

with open("tempfile.pickle", "wb") as f:
    pickle.dump(frag_and_mrna_seq, f)

with open('tempfile.pickle', 'rb') as handle:
    a = pickle.load(handle)

with open("overlap_pn_frags_mrna_seqs.tsv", "w") as outfile:
    writer = csv.writer(outfile, delimiter = "\t")
    for key, value in a.items():
        writer.writerow([key, value])

print(len(frag_and_mrna_seq.keys()))


# In[ ]:





# In[ ]:




