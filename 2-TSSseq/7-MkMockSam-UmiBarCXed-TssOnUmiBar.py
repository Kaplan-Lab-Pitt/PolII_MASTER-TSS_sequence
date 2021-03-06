### python 3 codes ###
# remotely use Pitt HTC cluster
# python=python/anaconda3.8-2020.11 (Python 3.8.5.final.0); pandas=1.1.3; numpy=1.19.2;

### To make mock SAM file for UMI-tools to correct actual_TSS sequence based on UMI-tools corrected 5'end_UMI+Barcode sequence

# Input-1: [DTmerge-MASTER-all_final-info_table.csv] containing information for all MASTER samples
# Input-2: [sam_header19.sam] containing first 19 lines of an example SAM file for S.c.
# Input-3s: UMI-tools output files [*-UmiCXed-mockBarCXonUmiTss-sorted-grouped.tsv] generated by <6-UMItools-UmiCXed-Bar_CX.slurm>
    # 9 columns
# Input-4s: 5'end_UMI-actlTSS Index file [*-UmiCXed-UmiTss_idx.txt] generated by <5-MkMockSam-UmiCXed-BarOnUmiTss.py>

# Output-1: 5'end_UMI+24bp_Barcode Index file [*-UmiBarCXed-UmiBar_idx.txt]
    # 2 columns, w/o titles: [Umi_Bar] \t [index]
# Output-2: mock SAM files [*-UmiBarCXed-mockTssCXonUmiBar.sam] for UMI-tools to correct actual_TSS sequence
# Output-3: [7-MkMockSam-UmiBarCXed-TssOnUmiBar-info.csv], details for processed samples

import datetime
job_start = datetime.datetime.now()
print(job_start)

import pandas as pd

# qname1 = 'Best wishes for all PhD students'
flag2 = 0 # "mapped to the forward strand"
rname3 = 'II'
# pos4 = 20160111 # 1st day in Kaplan lab :)
mapq5 = 60 # indicates uniquely mapped read
cigar6 = '21M'
rnext7 = '*'
pnext8 = 0
tlen9 = 0
seq10 = 'TCTATGGCTAGAACATATTAT'
qual11 = 'LikesSoccerScienceArt'

cmn_folder = '/bgfs/ckaplan/Yunye/0-common_files/'
inp_folder = '/bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/5_6-UmiCXed-BarOnUmiTss/'
out_folder = '/bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/7_8-UmiBarCXed-TssOnUmiBar/'

smp_info = pd.read_csv(cmn_folder+'DTmerge-MASTER-all_final-info_table.csv', na_filter= False) # w/o filling empty cells as NaN
# to process TSSseq_b1 of all final samples (3x3x5=45):
smp_info = smp_info[smp_info['final'] == 'yes']

# table to collect info for each sample
op_info = smp_info[['super_id','PolII','lib','rep']].set_index('super_id')

for _, smp in smp_info.iterrows():
    file_prefix = smp['lib']+'_'+smp['PolII']+'_'+smp['rep']+'_b1-'
    op_info.loc[smp['super_id'],'file_prefix'] = file_prefix
    op_info.loc[smp['super_id'],'input_mtx'] = file_prefix+'UmiCXed-mockBarCXonUmiTss-sorted-grouped.tsv'
    op_info.loc[smp['super_id'],'input_idx'] = file_prefix+'UmiCXed-UmiTss_idx.txt'

    ## Get 5'end_UMI-actlTSS infomation back based on index files
    UmiTss_idx = dict([idx_line.split() for idx_line in open(inp_folder+file_prefix+'UmiCXed-UmiTss_idx.txt', 'r')])
    idx_UmiTss = dict((v,k) for (k,v) in UmiTss_idx.items())

    grp = pd.read_csv(inp_folder+file_prefix+'UmiCXed-mockBarCXonUmiTss-sorted-grouped.tsv', sep='\t', index_col=False)\
            [['position','final_umi','final_umi_count']].drop_duplicates()
    grp['Umi_actlTSS'] = ((grp['position']+1).apply(str)).map(idx_UmiTss)
        # +1 because 1) Alignment position is not the start position of the read in the BAM file but the start of the read taking into account the read strand and cigar
        #            2) when making mock SAM file, "0" was put in FLAG field, which means "mapped to the forward strand"
        #            3) UMI-tools uses 0-indexed but SAM uses 1-indexed
    grp['CXed_Umi'] = grp['Umi_actlTSS'].str.slice(start=0, stop=15)
    grp['actlTSS'] = grp['Umi_actlTSS'].str.slice(start=15)
    grp['CXed_Bar'] = grp['final_umi']
    grp['RNA_Count'] = grp['final_umi_count']
 
    ## To make and save Umi_Bar-Index file
    UmiBar = grp[['CXed_Umi', 'CXed_Bar']].drop_duplicates().reset_index().drop('index', axis=1)
    UmiBar.index = UmiBar.index + 1 # SAM uses 1-indexed
    UmiBar['CXedUmi_CXedBar'] = UmiBar['CXed_Umi'] + UmiBar['CXed_Bar']
    UmiBar.reset_index()[['CXedUmi_CXedBar','index']].to_csv(out_folder+file_prefix+'UmiBarCXed-UmiBar_idx.txt', sep='\t', index=False, header=False)

    ## To make mock SAM file
    UmiBar_idx = dict([idx_line2.split() for idx_line2 in open(out_folder+file_prefix+'UmiBarCXed-UmiBar_idx.txt', 'r')])

    fout = open(out_folder+file_prefix+'UmiBarCXed-mockTssCXonUmiBar.sam', 'w')
    
    samHeader = [sam_line.rstrip() for sam_line in open(cmn_folder+'sam_header19.sam', 'r')]
    print("\n".join(samHeader), file = fout)

    row_num = 0 # to record # of Umi-actlTSS-Barcode variants
    readIn = 0
    readOut = 0
    grp = grp[['CXed_Umi','actlTSS','CXed_Bar','RNA_Count']]
    for index, row in grp.iterrows():
        # row['CXed_Umi']; row['actlTSS']; row['CXed_Bar']; row['RNA_Count']
        
        row_num += 1
        
        # UMI-tools requires same length of "UMI" to be CXed, but actual_TSS has various lengthes,
        # so need to add "stuff" sequence 'Z'
        
        # Maxiumn length of sequencing reads:
            # TSSseq_b1: WT-151, E1103G-151, F1086S-155, G1097D-151, H1085Q-151
            # TSSseq_b2 (NovaSeq): 201
        # Maximum length of "actlTSS" = max_read_len -15-24-27, meaning
            # [5End_15ntUMI] + [actlTSS] + [24bpBarcode_v2] + [ATGTC+CKO2191] = 15 + (1~max) + 24 + 27
            # Therefore,
                # TSSseq_b1: theoretical max= 155-15-24-27= 89; actual max= 89
                # TSSseq_b2: theoretical max= 201-15-24-27= 135; actual max= 113
        
        max_len = 89
        stuffed_Z_num = max_len-len(row['actlTSS'])
        actlTSS_stuffed = row['actlTSS'] + 'Z'*stuffed_Z_num # add 'Z's at the end to get equal length

        readIn += int(row['RNA_Count'])
        for _ in range(int(row['RNA_Count'])):
            readOut += 1
            qname1 = 'Mock'+ str(readOut)+ 'Z'+ str(stuffed_Z_num) + '_' + actlTSS_stuffed
            print("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s" % \
                  (qname1, flag2,rname3, int(UmiBar_idx[row['CXed_Umi']+row['CXed_Bar']]), mapq5,cigar6,rnext7,pnext8,tlen9,seq10,qual11), file = fout)

    op_info.loc[smp['super_id'],'n_reads_input'] = readIn
    op_info.loc[smp['super_id'],'n_reads_output'] = readOut
    op_info.loc[smp['super_id'],'n_Umi_actlTSS_Bar_variants_input'] = row_num

op_info.reset_index().to_csv('7-MkMockSam-UmiBarCXed-TssOnUmiBar-info.csv', sep=',', index=False, header=True, mode='w')
    
job_end = datetime.datetime.now()
print("job finished in", datetime.timedelta.total_seconds(job_end - job_start)/60, 'mins at', job_end)
