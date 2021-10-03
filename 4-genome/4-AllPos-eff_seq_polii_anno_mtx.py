### python 3 codes ###
# remotely use Pitt HTC cluster
# python=python/anaconda3.8-2020.11 (Python 3.8.5.final.0); pandas=1.1.3; numpy=1.19.2;

### to generate matrix for each/all positions in 5979x401 window,
### including information of efficiency, sequence at positions -11 to +9, some annotations

# Input-1: [gnm_TSSseq-info_table.csv] containing updated information for all genome TSS-seq samples
# Input-2: [TSScontext_401x5979.txt], TSS context from positions -11 to +9 for each position within 401 windows, generated by <1-extrTSScontext.py>
    # 407 columns = [chr]-[start]-[end]-[strand]-[attribute]-[401PW_seq] - [up250]...[medTSS]...[dn150]
# Input-3: [Master_Annotation_Columns_uniq_5979.xlsx], annotation of 5979 promoter windows, from Kaplan lab
# Input-4s: [*-401x5979_effMtx.csv], 401x5979 efficiency matrix generated by <3-CntToEff.py>
# Input-5s: 401x5979 matrix of TSS-seq reads count, from Kaplan lab's previous study (Zhao et al., 2021)

# Output: [*-AllPos-eff_seqn11p9_anno-mtrx.csv], matrix with 1+20+10 columne: efficiency of every position from 401x5979, base at each -11 to +9 positions, Pol WT or mut, distance to medTSS (-250 to +150), TSS context, some other information
    # [eff] - ['n11','n10','n9','n8','n7','n6','n5','n4','n3','n2','n1','p1','p2','p3','p4','p5','p6','p7','p8','p9']x20 - [PolII]-[PWpos]-[distToMed]-[TSSn11p9]-[annotations]x4-[ttlRNAcnt]-[RNAcnt]

import datetime
job_start = datetime.datetime.now()
print(job_start)

import pandas as pd

cmn_folder = '/bgfs/ckaplan/Yunye/0-common_files/'
inp_folder_eff = '/bgfs/ckaplan/Yunye/3-TSS_sequence_library/9-gnm/3-CntToEff/' # Input-4s
inp_folder_cnt = '/bgfs/ckaplan/Yunye/3-TSS_sequence_library/9-gnm/0-cnt_mtx/' # Input-5s
out_folder = '/bgfs/ckaplan/Yunye/3-TSS_sequence_library/9-gnm/4-AllPos_mtx/'

smp_info = pd.read_csv(cmn_folder+'gnm_TSSseq-info_table.csv')
# to only process certain samples
smp_info = smp_info.loc[smp_info['super_id'].isin([25])]

context = pd.read_csv(cmn_folder+'TSScontext_401x5979.txt', sep = '\t', header = 0) # Input 2
anno = pd.read_excel(cmn_folder+'Master_Annotation_Columns_uniq_5979.xlsx', header=0) # Input 3

for _, smp in smp_info.iterrows():
    file_prefix = smp['Grad']+'_'+smp['PolII']+'_'+smp['mrgd_indv']+'-'
    print("Notes for processing sample: ", smp['Notes'])
    print("output prefix:", file_prefix)

    pos_eff_mtx = pd.read_csv(inp_folder_eff+smp['5979x401_eff_mtx'], header=None, index_col=None)
    pos_cnt_mtx = pd.read_csv(inp_folder_cnt+smp['5979x401_cnt_mtx'], header=None, index_col=None)
    
    mtrx = pd.DataFrame()
    for p in range(0,401):
        pMtrx = pd.DataFrame()
        pMtrx['eff'] = pos_eff_mtx.iloc[:, p]
        for j, pos in enumerate(['n11','n10','n9','n8','n7','n6','n5','n4','n3','n2','n1','p1','p2','p3','p4','p5','p6','p7','p8','p9']):
            pMtrx[pos] = context.iloc[:, (6+p)].str.slice(start=j, stop=j+1)
        pMtrx['PolII'] = smp['PolII']
        pMtrx['PWpos'] = context.columns[6+p]
        pMtrx['distToMed'] = p-250
        pMtrx['TSSn11p9'] = context.iloc[:, (6+p)]
        pMtrx['INDX6044'] = anno['INDX']
        pMtrx['gene'] = anno['gene']
        pMtrx['TATAclass'] = anno['TATA_class corrected']
        pMtrx['TAF1status'] = anno['TAF1_status']
        pMtrx['ttlRNAcnt'] = pos_cnt_mtx.sum(axis=1) # total TSS-seq reads cnt for each promoter window
        pMtrx['RNAcnt'] = pos_cnt_mtx.iloc[:, p] # TSS-seq reads cnt for a particular position
        
        mtrx = mtrx.append(pMtrx, ignore_index=True)

    mtrx.to_csv(out_folder+file_prefix+'AllPos-eff_seqn11p9_anno-mtrx.csv', index=False, header=True)
    
job_end = datetime.datetime.now()
print("job finished in", datetime.timedelta.total_seconds(job_end - job_start)/60, 'mins at', job_end)