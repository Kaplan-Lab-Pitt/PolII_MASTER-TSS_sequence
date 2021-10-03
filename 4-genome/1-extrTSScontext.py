### python 2 codes ###

### to extract 20nt TSS context (positions between -11 to +9) for every position in 401 windows of 5979 promoters version
# Input-1: [S288C_reference_sequence_R64-1-1_20110203.fsa], reference genome
# Input-2: [250Up_agMed_uniq_150Down_5979.gff], gff file for known 5979x401 promoter windows (Qiu et al., 2020)

# Output-1: [TSScontext_401x5979.txt], matrix w/ 5+1+401= 407 columns:
    # [chr]-[start]-[end]-[strand]-[attribute]-[PW_seq]-[n250]-...-[medTSS]-[p1]-...[p150]

import pandas as pd
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

## to store genome reference sequence for each chromosome into a dict
ref = open("S288C_reference_sequence_R64-1-1_20110203.fsa",'r').readlines()
chr_indx = [index for index,value in enumerate(ref) if value.startswith('>ref')]
chr_ref = dict() # to store ref sequence for each chromosome
for i in range(16): # 16 yeast chromosomes
    chr_name = 'chr' + ref[chr_indx[i]][(ref[chr_indx[i]].find('chromosome=')+11):ref[chr_indx[i]].find(']\n')]
    chr_seq = ''.join(ref[chr_indx[i]+1:chr_indx[i+1]]).replace('\n','')
    chr_ref[chr_name] = chr_seq

fout = open('TSScontext_401x5979.txt', 'w') # Output-1
print >> fout, '\t'.join(['chr','start','end','strand','attribute','401PW_seq'] + ['up'+str(x) for x in range(250,0,-1)] + ['medTSS'] + ['dn'+str(x) for x in range(1,151,1)])

pW_len = 401 # promoter window range
context_start = -11 # start from position -11
context_end = +9 # until position +9
PWs = pd.read_csv('250Up_agMed_uniq_150Down_5979.gff', delimiter='\t', header=None,\
                    names=['chr','source','feature','start','end','unk1','strand','unk2','attribute']) # gff for promoter windows
for _, row in PWs.iterrows():
    if row['strand'] == '+':
        full_seq = chr_ref[row['chr']][(row['start']-1+context_start):(row['start']-1+pW_len+(context_end-1))]
            # row['start']-1 : python is 0-based; GFF is 1-based
            # (row['start']-1)+context_start: need seq from context_start position of up250
            # (row['start']-1+pW_len)+context_end-1: need seq until context_end position of dn150, but dn150 itself is +1
    elif row['strand'] == '-':
        full_seq = str(Seq(chr_ref[row['chr']][(row['start']-1-(context_end-1)):(row['start']-1+pW_len+(-context_start))], IUPAC.unambiguous_dna).reverse_complement())
    
    out = [row['chr'], str(row['start']), str(row['end']), str(row['strand']), row['attribute']]
    out = out + [full_seq[(-context_start):(-context_start+pW_len)]] # add [PW_seq] that is 401nt
    out = out + [full_seq[(p+context_start):(p+context_end)] for p in range(-context_start, -context_start+pW_len)] # add each TSS context
    
    print >> fout, '\t'.join(out)
