### python 3 codes ###
# remotely use Pitt HTC cluster
# python=python/anaconda3.8-2020.11 (Python 3.8.5.final.0); pandas=1.1.3; numpy=1.19.2;

# After combining TSS-seq batch 1&2 datasets
# 1. get matched #RNA reads info from 3 yeast reps for each promoter variant
# 2. calculate stat (exit in #rep, Sum, SD, Average, CV)

import pandas as pd

cmn_folder = '/bgfs/ckaplan/Yunye/0-common_files/'
inp_folder = '/bgfs/ckaplan/Yunye/3-TSS_sequence_library/7-poolDs_cmbTs/'
out_folder = '/bgfs/ckaplan/Yunye/3-TSS_sequence_library/7-poolDs_cmbTs/2p-reps_RNAstat/'

smp_info = pd.read_csv(cmn_folder+'DTmerge-MASTER-all_final-info_table.csv', na_filter= False) # w/o filling empty cells as NaN
# to process all final samples (3x3x5=45):
smp_info = smp_info[smp_info['final'] == 'yes']

for lib in ['AYR','BYR','ARY']:
    for pol in ['WT','E1103G','F1086S','G1097D','H1085Q']:
        mtx = pd.read_csv(cmn_folder+lib+'_TSS.csv', header=0, index_col=0)
        for _, smp in smp_info[(smp_info['lib']==lib) & (smp_info['PolII']==pol)].iterrows():
            file_prefix = smp['lib']+'_'+smp['PolII']+'_'+smp['rep']
            mtx[smp['rep']] = pd.read_csv(inp_folder+'2-B12Cmb-pos_cnt_mtx/'+file_prefix+'-pDs-B12Cmb-pos_cnt_mtx.csv',\
                                          header=0, index_col=0).sum(axis=1)
        mtx['n_exit']= 3-mtx.iloc[:,0:3].isnull().sum(axis=1)
        mtx['Sum'] = mtx.iloc[:,0:3].sum(axis=1)
        mtx['SD'] = mtx.iloc[:,0:3].std(axis=1)
        mtx['Average'] = mtx.iloc[:,0:3].mean(axis=1)
        mtx['CV'] = mtx['SD'] / mtx['Average']
                
        mtx.reset_index().to_csv(out_folder+lib+'_'+pol+'-pDs_B12Cmb-reps_RNAstat.csv', index=False, header=True)