### python 3 codes ###
# remotely use Pitt HTC cluster
# python=python/anaconda3.8-2020.11 (Python 3.8.5.final.0); pandas=1.1.3; numpy=1.19.2;

### To make mock SAM file for UMI-tools to correct 5'end_UMI sequence based on actual_TSS+Barcode_v2 sequence

# Input-1: [DTmerge-MASTER-all_final-info_table.csv] containing information for all MASTER samples
# Input-2: [sam_header19.sam] containing first 19 lines of an example SAM file for S.c.
# Input-3s: matrix files [*-Umi_actlTss_Bar_Cnt-Mtx.txt] generated by <2-Umi_actlTss_Bar_Cnt-Mtx.py>
    # 4 columns, sep = '\t': [5UMI]-[actlTSS]-[Barcode24]-[RNA_Count]

# Output-1s: actlTSS-Barcode Index file [*-TssBar_idx.txt]
    # 2 columns, w/o titles, sep = '\t': [actlTSS_Bar]-[index]
# Output-2s: [*-mockUmiCXonTssBar.sam], mock SAM files for UMI-tools to correct 5'end_15bp_UMI based on actual_TSS_region+24bp_Barcode
# Output-3: [3-MkMockSam-UmiOnTssBar-info.csv], details for processed samples

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
inp_folder = '/bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/2-Umi_actlTss_Bar_Cnt-Mtx/'
out_folder = '/bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/3_4-UmiOnTssBar/'

smp_info = pd.read_csv(cmn_folder+'DTmerge-MASTER-all_final-info_table.csv', na_filter= False) # w/o filling empty cells as NaN
# to process TSSseq_b1 of all final samples (3x3x5=45):
smp_info = smp_info[smp_info['final'] == 'yes']

# table to collect info for each sample
op_info = smp_info[['super_id','PolII','lib','rep']].set_index('super_id')

for _, smp in smp_info.iterrows():
    file_prefix = smp['lib']+'_'+smp['PolII']+'_'+smp['rep']+'_b1-'
    op_info.loc[smp['super_id'],'file_prefix'] = file_prefix
    op_info.loc[smp['super_id'],'input_mtx'] = file_prefix+'Umi_actlTss_Bar_Cnt-Mtx.txt'

    ## To make and save actual_TSS_region-24bp_Bar Index file
    ActlTSSBar = pd.read_csv(inp_folder+file_prefix+'Umi_actlTss_Bar_Cnt-Mtx.txt', sep='\t', index_col=False)\
                   [['actlTSS', 'Barcode24']]\
                   .drop_duplicates().reset_index().drop('index', axis=1)
    ActlTSSBar.index = ActlTSSBar.index + 1
    ActlTSSBar['actlTSS_Bar'] = ActlTSSBar['actlTSS'] + ActlTSSBar['Barcode24']
    ActlTSSBar.reset_index()[['actlTSS_Bar','index']]\
              .to_csv(out_folder+file_prefix+'TssBar_idx.txt', sep='\t', index=False, header=False)
    
    ## To make mock SAM file
    ActlTSSBar_idx = dict([idx_line.split() for idx_line in open(out_folder+file_prefix+'TssBar_idx.txt', 'r')])
    
    fout = open(out_folder+file_prefix+'mockUmiCXonTssBar.sam', 'w')
    
    samHeader = [sam_line.rstrip() for sam_line in open(cmn_folder+'sam_header19.sam', 'r')]
    print("\n".join(samHeader), file = fout)
    
    line_num = 0
    readIn = 0
    readOut = 0
    for line in open(inp_folder+file_prefix+'Umi_actlTss_Bar_Cnt-Mtx.txt', 'r'):
        line_num += 1
        if line_num > 1: # first line is titles
            mtx = line.split('\t') # mtx[0] - 5End_15ntUMI; mtx[1] - Actual_TSS_Region; mtx[2] - 24bp_Barcode_v2; mtx[3] - RNA_Count
            readIn += int(mtx[3])
            for _ in range(0, int(mtx[3])):
                readOut += 1
                qname1 = 'Mock'+ str(readOut) + '_' + mtx[0]
                print("%s\t%d\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%s\t%s" % \
                     (qname1, flag2,rname3, int(ActlTSSBar_idx[mtx[1]+mtx[2]]), mapq5,cigar6,rnext7,pnext8,tlen9,seq10,qual11), file = fout)

    op_info.loc[smp['super_id'],'n_reads_input'] = readIn
    op_info.loc[smp['super_id'],'n_reads_output'] = readOut
    op_info.loc[smp['super_id'],'n_Umi_actlTSS_Bar_variants_input'] = line_num-1

op_info.reset_index().to_csv('3-MkMockSam-UmiOnTssBar-info.csv', sep=',', index=False, header=True, mode='w')
    
job_end = datetime.datetime.now()
print("job finished in", datetime.timedelta.total_seconds(job_end - job_start)/60, 'mins at', job_end)
