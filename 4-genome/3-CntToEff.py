### python 3 codes ###
# remotely use Kaplan lab iMac
# python=3.8.6; pandas=1.1.3;

### to generate efficiency matrix for 401x5979 matrix
# definition: yield of a given position / the sum of yield at this position with the yield at all downstream position(s)
# Filter: to avoid the very last a few positions being highly efficient,
    # filter out positions with >=20% efficiency and with <=5 reads PolII flux
   
# Input-1: [gnm_TSSseq-info_table.csv] containing updated information for all genomic TSS-seq samples
# Input-2s: 401x5979 matrix of TSS-seq reads count, from Kaplan lab's previous study (Zhao et al., 2021)
    # row= 5979 yeast promoters
    # column= individual position within known promoter windows ("median" TSS, 250nt up and 150nt dn from median TSS)
    # cell= TSS-seq reads count

# Output: [*-401x5979_effMtx.csv], 401x5979 matrix of TSS efficiency

import datetime
job_start = datetime.datetime.now()
print(job_start)

import pandas as pd

smp_info_file = '/Users/Informatics/Informatics/Yunye/0-common_files/gnm_TSSseq-info_table.csv'
smp_info = pd.read_csv(smp_info_file)
# to only process certain samples
smp_info = smp_info.loc[smp_info['super_id'].isin([25,26,27])]

cnt_mtx_folder = '0-cnt_mtx/'
eff_mtx_folder = '3-CntToEff/'
for _, smp in smp_info.iterrows():
    file_prefix = smp['Grad']+'_'+smp['PolII']+'_'+smp['mrgd_indv']+'-'
    print("5979x401 matrix:", smp['5979x401_cnt_mtx'])
    print("Notes for processing sample: ", smp['Notes'])
    print("output prefix:", file_prefix)
    
    pos_cnt_mtx = pd.read_csv(cnt_mtx_folder+smp['5979x401_cnt_mtx'], header=None, index_col=None)

    pos_eff_mtx = pd.DataFrame(index=range(pos_cnt_mtx.shape[0]), columns=range(pos_cnt_mtx.shape[1]))
    for idx, pos_count in pos_cnt_mtx.iterrows():
        row_list = pos_count.tolist()
        pos_eff = [None] * pos_cnt_mtx.shape[1]
        for k in range(pos_cnt_mtx.shape[1]):
            if sum(row_list[:k]) != sum(row_list): # still polII flux remaining
                eff = float(row_list[k])/float(sum(row_list[k:]))*100  # calculate TSS efficiency
                if not ((eff >= 20) and (sum(row_list[k:])<=5)):
                    pos_eff[k] = eff
                else: # filter out positions with >=20% efficiency with <=5 reads PolII flux
                    break
            else:
                break
        pos_eff_mtx.iloc[idx] = pos_eff

    pos_eff_mtx.to_csv(eff_mtx_folder+file_prefix+'401x5979_effMtx.csv', index=False, header=False)

job_end = datetime.datetime.now()
print("job finished in", datetime.timedelta.total_seconds(job_end - job_start)/60, 'mins at', job_end)
