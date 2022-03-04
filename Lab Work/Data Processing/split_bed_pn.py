#!/usr/bin/env python
# coding: utf-8

# In[3]:


def check_op_strands():
    txt_path = 'all_overlaps_sf.bed'
    f = open(txt_path, 'r')
    f.readline()
    pos_bed_path = 'pos_overlaps.bed'
    neg_bed_path = 'neg_overlaps.bed'
    p_bed = open(pos_bed_path, 'w')
    n_bed = open(neg_bed_path, 'w')
    for line in f:
        line = line.strip().split('\t')    
        if line[4] == "+":
            p_bed.write(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + '\t' + line[4] + '\t' + line[5] + '\t' + line[6] + '\t'  + line[7] + '\t' + line[8] + '\t'+ line[9] + '\t'+ line[10] + '\t'+ line[11] + '\n')
        elif line[4] == '-':
            n_bed.write(line[0] + '\t' + line[1] + '\t' + line[2] + '\t' + line[3] + '\t' + line[4] + '\t' + line[5] + '\t' + line[6] + '\t'  + line[7] + '\t' + line[8] + '\t'+ line[9] + '\t'+ line[10] + '\t'+ line[11] + '\n')
    return 


# In[3]:


check_op_strands()


# In[ ]:





# In[ ]:




