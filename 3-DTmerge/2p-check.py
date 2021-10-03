### python 3 codes ###
# remotely use Pitt HTC cluster
# python=python/anaconda3.8-2020.11 (Python 3.8.5.final.0); pandas=1.1.3; numpy=1.19.2;

### to make mtx to check reps reproducibility using TSS efficiency of major TSS of each library
# row = TSS variants (promoter variants)
# column = samples (15 = reps x WT/mutant)
# cell = TSS efficiency of major TSS

import pandas as pd

cmn_folder = '/bgfs/ckaplan/Yunye/0-common_files/'
inp_folder = '/bgfs/ckaplan/Yunye/3-TSS_sequence_library/7-poolDs_cmbTs/'
out_folder = '/bgfs/ckaplan/Yunye/3-TSS_sequence_library/7-poolDs_cmbTs/2p-check/'

smp_info = pd.read_csv(cmn_folder+'DTmerge-MASTER-all_final-info_table.csv', na_filter= False) # w/o filling empty cells as NaN
# to process all final samples (3x3x5=45):
smp_info = smp_info[smp_info['final'] == 'yes']

for (lib, tss) in [('AYR','1'), ('BYR','1'), ('ARY','2')]:
    mtx = pd.read_csv(cmn_folder+lib+'_TSS.csv', header=0, index_col=0)

    for _, smp in smp_info[(smp_info['lib']==lib)].iterrows():
        file_prefix = smp['lib']+'_'+smp['PolII']+'_'+smp['rep']
        mtx[file_prefix] = pd.read_csv(inp_folder+'2-B12Cmb-pos_eff_mtx/'+file_prefix+'-pDs-B12Cmb-pos_eff_mtx.csv',\
                                       header=0, index_col=0)[tss]

    mtx.reset_index().to_csv(out_folder+lib+'_final-reps_pDs_B12Cmb-TSS_'+tss+'eff.csv', index=False, header=True)
