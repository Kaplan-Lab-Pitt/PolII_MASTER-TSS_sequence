### python 3 codes ###
# remotely use Pitt HTC cluster
# python=python/anaconda3.8-2020.11 (Python 3.8.5.final.0); pandas=1.1.3; numpy=1.19.2;

### merge DNA-seq and TSS-seq data
# use "inner" merge, meaning only keep Barcodes that exist in both DNAseq and RNAseq
# seperately merge TSS-seq from 2 batches

## Before merging DNA-seq and RNA-seq data,
# DNA-seq:
    # 1) TSS & Barcode extraction
    # 2) TSS & Barcode individual correction using UMI-tools
    # 3) Barcode_multiple_TSS filter: >= 90%
    # 4) pool DNAseq matrix based on E.coli and 3 WT Yreps:
    #    A. keep TSS-Barcode variants that
    #       (a) >=5 reads in libraries
    #       (b) exist in >=2 Ecoli+3Yreps libraries, meaning missing in 0 or 1 or 2 libraries
    # ==> 6 columns, sep='\t': [TSS] - [Barcode_Ds] - [Ec] - [rep1] - [rep2] - [rep3]

# TSS-seq:
    # 1) 5'end_15nt_UMI & actual_TSS_region & Barcode extraction
    # 2) 5'end_15nt_UMI & actual_TSS_region & Barcode individual correction using UMI-tools
    # 3) RNA_Count deduplication: count each UMI_actlTSS_Bar variant as "1" RNA, regardless of RNA_Count
        # then Accumulate "real" RNA_Count for each actual_TSS_region + Barcode variant (i.e. matches to how many different 5'end_15nt_barcodes);
    # ==> 3 columns, sep='\t': [actlTSS] - [Barcode24] - [RNA_ddpCount]

# Input-1: [DTmerge-MASTER-all_final-info_table.csv] containing updated information for all MASTER samples
# Input-2s: [*-poolDs-5cnt_2EY_ct-TSS_Bar_EYCnt-mtx.txt], pooled DNAseq matrix
# Input-3s: [*-CXed-actlTss_Bar_ddpCnt-mtrx.txt], processed TSSseq matrix

# Output-1s: [*-pos_cnt_mtx.csv] TSS usage matrix for positions between -68 to +25, relative to designed +1 TSS as +1
# Output-2s: []*-pos_eff_mtx.csv] TSS efficiency matrix for positions between -68 to +25, relative to designed +1 TSS as +1
# Output-3: [1-merge_poolDs_TsB12-info.csv-info.csv] table collecting info/stat for each sample
# Other outputs are documented in scripts

import datetime
job_start = datetime.datetime.now()
print(job_start)

import pandas as pd

cmn_folder = '/bgfs/ckaplan/Yunye/0-common_files/'
inp_folder_Ds = '/bgfs/ckaplan/Yunye/3-TSS_sequence_library/1-WT-DNAseq-TxGen_18179Kap/8-poolDseq-5cnt_2EY_ct/'
inp_folder_TsB1 = '/bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/9-ddpUmi/'
inp_folder_TsB2 = '/bgfs/ckaplan/Yunye/3-TSS_sequence_library/6-reSeq_NovaSeq/2-extr_CX_ddp/9-ddpUmi/'
out_folder = '/bgfs/ckaplan/Yunye/3-TSS_sequence_library/7-poolDs_cmbTs/'

smp_info = pd.read_csv(cmn_folder+'DTmerge-MASTER-all_final-info_table.csv', na_filter= False) # w/o filling empty cells as NaN
# to process all final samples (3x3x5=45):
smp_info = smp_info[smp_info['final'] == 'yes']

# table to collect info for each sample
op_info = smp_info[['super_id','PolII','lib','rep']].set_index('super_id')

for _, smp in smp_info.iterrows():
    op_info.loc[smp['super_id'],'input_Dseq_mtx'] = smp['lib']+'-poolDs-5cnt_2EY_ct-TSS_Bar_EYCnt-mtx.txt'
    for b in range(int(smp['Tseq-ttlBth'])): # 27 samples have TSS-seq batch_2
        b = b+1 # batch#
        file_prefix = smp['lib']+'_'+smp['PolII']+'_'+smp['rep']+'_b'+str(b)+'-pDs-'
        op_info.loc[smp['super_id'],'Ts_b'+str(b)+'-input_Tseq_mtx'] = smp['lib']+'_'+smp['PolII']+'_'+smp['rep']+'_b'+str(b)+'-CXed-actlTss_Bar_ddpCnt-mtrx.txt'
        op_info.loc[smp['super_id'],'Ts_b'+str(b)+'-file_prefix'] = file_prefix

        DNAseq = pd.read_csv(inp_folder_Ds+smp['lib']+'-poolDs-5cnt_2EY_ct-TSS_Bar_EYCnt-mtx.txt',\
                             sep = '\t', header=0, index_col=False)[['TSS','Barcode_Ds']]
        if b==1: inp_folder_Ts=inp_folder_TsB1
        elif b==2: inp_folder_Ts=inp_folder_TsB2
        RNAseq = pd.read_csv(inp_folder_Ts+smp['lib']+'_'+smp['PolII']+'_'+smp['rep']+'_b'+str(b)+'-CXed-actlTss_Bar_ddpCnt-mtrx.txt',\
                             sep = '\t', header=0, index_col=False)

        ### inner merge to only keep barcodes that could be found in both DNAseq and TSSseq
        DNA_RNA_merged_inner = pd.merge(DNAseq, RNAseq, left_on = 'Barcode_Ds', right_on = 'Barcode24', how = 'inner')
            # each DNAseq TSS_Bar_DNA could match to multiple RNAseq actlTSS_Bar_RNA var
            # but each RNAseq actlTSS_Bar_RNA var should only match to 1 DNAseq TSS_Bar_DNA var

        op_info.loc[smp['super_id'],'Ts_b'+str(b)+'-n_Bar_in_bothDR'] = DNA_RNA_merged_inner.groupby('Barcode_Ds').count().shape[0]
        op_info.loc[smp['super_id'],'Ts_b'+str(b)+'-n_reads_Rseq'] = DNA_RNA_merged_inner['RNA_ddpCount'].sum()
        op_info.loc[smp['super_id'],'Ts_b'+str(b)+'-n_raw_TSS_var'] = DNA_RNA_merged_inner.groupby('TSS').count().shape[0]

        # to generate RNA_ddpCount distribution table after inner merge
        DNA_RNA_merged_inner.groupby('RNA_ddpCount').count()['Barcode24']\
                            .to_csv(out_folder+'1-stat/'+file_prefix+'mrgd-RNA_ddpCount-dist.csv', index=True)

        # to pool DNAseq info by TSS
            # keep 'count',lambda x: tuple(x) info of "Barcode_Ds"
        DNAseq_info = DNAseq[DNAseq['Barcode_Ds'].isin(DNA_RNA_merged_inner['Barcode_Ds'])]\
                           .groupby('TSS').agg(['count',lambda x: tuple(x)])

        # to pool TSSseq info by TSS 
            # keep 'sum','count',lambda x: tuple(x) info of "actlTSS"
        DNA_RNA_merged_inner_RNAseq_scl = DNA_RNA_merged_inner[['TSS','actlTSS','RNA_ddpCount']]\
                                                              .groupby(['TSS', 'actlTSS']).sum().reset_index()\
                                                              .sort_values(by='RNA_ddpCount',ascending=False)\
                                                              .groupby('TSS').agg(['sum','count',lambda x: tuple(x)])
        DNA_RNA_merged_inner_RNAseq_info = DNA_RNA_merged_inner_RNAseq_scl.drop(DNA_RNA_merged_inner_RNAseq_scl.columns[[0,1]], axis=1)

        # to merge pooled DNAseq & TSSseq info
        final_merge = pd.merge(DNAseq_info, DNA_RNA_merged_inner_RNAseq_info, left_index=True, right_index=True) # both _info mtx use TSS as index

        final_merge.columns = pd.Index([idx[0] + '-' + idx[1] for idx in final_merge.columns.tolist()])
            # to convert 2-levels of column indexes to 1-level, using combined level1-level2 name as new index name
        final_merge.rename(columns={'Barcode_Ds-<lambda_0>':'Barcode_Ds-list',\
                                    'actlTSS-<lambda_0>':'actlTSS-list', 'RNA_ddpCount-<lambda_0>':'RNA_ddpCount-list'}, inplace=True)

        final_merge.reset_index().to_csv(out_folder+'1-aggr_info/'+file_prefix+'mrgd-aggr_info.csv', index=False)

        ### to generate TSS usage matrix, containing positions between -68 to +25, relative to designed +1 TSS as +1
        pos_cnt_mtx = pd.DataFrame(columns=[pc for pc in range(-68,26) if pc != 0])
        rnas_w_mis = pd.DataFrame(columns=['TSS', 'RNAs_with_mismatch-list', 'RNAs_with_mismatch-Counts_list', 'RNAs_with_mismatch-Counts_sum', 'RNA_w/o_mismatch-Counts_sum'])
    
        for TSS, row in final_merge.iterrows():
            refer_seq = 'AAACTTCTTCCCTTTGTTACTTCTTCTTTAAAATTCAAATTTTTCTTTTGATTTTTTTTC'+TSS+'ACATTTTCAAAAGGCTAACATCAG' # [-68 to -9] - [TSS -8 to +1] - [FD]
            pos_reads = [0] * 93 # list where each item represents one position between -68 to +25
            mis_seq = [] # to save RNA seq w/ mismatch
            mis_count =[] # to save reads count of RNA seq w/ mismatch, same order as mis_seq list
        
            for j in range(row['RNA_ddpCount-count']): # row['RNA_ddpCount-count'] equals how many different actlTSS variants for this TSS, then go through each actlTSS
                seq = row['actlTSS-list'][j]
                count = row['RNA_ddpCount-list'][j]
                if seq == refer_seq[-len(seq):]: # perfect match
                    pos_reads[-len(seq)] += count # add count value to corresponding TSS position
                else: # contains slippage or regular mismatch
                    mis_seq.append(seq)
                    mis_count.append(count)
               
            if sum(pos_reads) != 0: # at least one perfect match
                pos_cnt_mtx.loc[TSS] = pos_reads
            if len(mis_seq)>0:
                mis_info = pd.Series({'TSS':TSS, 'RNAs_with_mismatch-list':mis_seq, 'RNAs_with_mismatch-Counts_list':mis_count,\
                                      'RNAs_with_mismatch-Counts_sum':sum(mis_count), 'RNA_w/o_mismatch-Counts_sum':sum(pos_reads)})
                rnas_w_mis = rnas_w_mis.append(mis_info, ignore_index=True)
    
        pos_cnt_mtx.to_csv(out_folder+'1-pos_cnt_mtx/'+file_prefix+'pos_cnt_mtx.csv')
        op_info.loc[smp['super_id'],'Ts_b'+str(b)+'-n_RNAddpCount_prf_mtch'] = pos_cnt_mtx.sum(axis=1).sum()
        rnas_w_mis.to_csv(out_folder+'1-aggr_info/'+file_prefix+'Reads_w_mismatch.csv', index=False)
    
        ### to generate TSS efficiency matrix
            # definition: yield of a given position / the sum of yield at this position with the yield at all downstream position(s)
            # Filter: to avoid the very last a few positions being highly efficient,
                # filter out positions with >=20% efficiency and with <=5 reads PolII flux
        pos_eff_mtx = pd.DataFrame(columns=[pe for pe in range(-68,26) if pe != 0])

        for TSS, pos_count in pos_cnt_mtx.iterrows():
            row_list = pos_count.tolist()
            pos_eff = [None] * 93
            for k in range(len(row_list)):
                if sum(row_list[:k]) != sum(row_list): # still polII flux remaining
                    eff = float(row_list[k])/float(sum(row_list[k:]))*100  # calculate TSS efficiency
                    if not ((eff >= 20) and (sum(row_list[k:])<=5)):
                        pos_eff[k] = eff
                    else: # filter out positions with >=20% efficiency with <=5 reads PolII flux
                        break
                else:
                    break
            pos_eff_mtx.loc[TSS] = pos_eff

        pos_eff_mtx.to_csv(out_folder+'1-pos_eff_mtx/'+file_prefix+'pos_eff_mtx.csv')

        pos_eff_mtx['mrgd_prf_ttlRNAcnt'] = pos_cnt_mtx.sum(axis=1)
        pos_eff_mtx.to_csv(out_folder+'1-pos_eff_mtx/'+file_prefix+'pos_eff_mtx-wTtl.csv')

op_info.reset_index().to_csv('1-merge_poolDs_TsB12-info.csv', sep=',', index=False, header=True, mode='w')

job_end = datetime.datetime.now()
print("job finished in", datetime.timedelta.total_seconds(job_end - job_start)/60, 'mins at', job_end)
