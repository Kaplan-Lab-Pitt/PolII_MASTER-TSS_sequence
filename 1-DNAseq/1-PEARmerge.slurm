#!/bin/bash
#SBATCH --job-name=1-PEARmerge
#SBATCH --output=1-PEARmerge.out
#SBATCH -t 3-00:00
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=xxxxxx@pitt.edu
#SBATCH --mail-type=BEGIN,END,FAIL

### To merge paired-end reads of DNA-seq using PEAR (0.9.11)
# remotely use Pitt HTC cluster

# Inputs: raw sequencing files for R1 and R2 [*_R1_001.fastq.gz] and [*_R2_001.fastq.gz]
#   two data sets in two folders:
#       fastq_18179Kap_N18162_L002
#       fastq_18179Kap_N18164_L002

# Outputs: [*-Ds-R1R2mrgd*]
#   file prefix naming: library_mutant_replicate-Ds-R1R2mrgd

module purge
module load pear/0.9.11

# process files in fastq_18179Kap_N18162_L002 folder
cd /bgfs/ckaplan/Yunye/3-TSS_sequence_library/1-WT-DNAseq-TxGen_18179Kap/fastq_18179Kap_N18162_L002/

pear -f ARY_E_coli_v1_S8_L002_R1_001.fastq.gz -r ARY_E_coli_v1_S8_L002_R2_001.fastq.gz -o ARY__Ec-Ds-R1R2mrgd &
pear -f ARY_Yeast_v1_1_Rep_1_S5_L002_R1_001.fastq.gz -r ARY_Yeast_v1_1_Rep_1_S5_L002_R2_001.fastq.gz -o ARY_WT_rep1-Ds-R1R2mrgd &
pear -f ARY_Yeast_v1_1_Rep_2_S6_L002_R1_001.fastq.gz -r ARY_Yeast_v1_1_Rep_2_S6_L002_R2_001.fastq.gz -o ARY_WT_rep2-Ds-R1R2mrgd &
pear -f ARY_Yeast_v1_1_Rep_3_S7_L002_R1_001.fastq.gz -r ARY_Yeast_v1_1_Rep_3_S7_L002_R2_001.fastq.gz -o ARY_WT_rep3-Ds-R1R2mrgd &
pear -f AYR_E_coli_v4_S4_L002_R1_001.fastq.gz -r AYR_E_coli_v4_S4_L002_R2_001.fastq.gz -o AYR__Ec-Ds-R1R2mrgd &
pear -f AYR_Yeast_v4_1_Rep_1_S1_L002_R1_001.fastq.gz -r AYR_Yeast_v4_1_Rep_1_S1_L002_R2_001.fastq.gz -o AYR_WT_rep1-Ds-R1R2mrgd &
pear -f AYR_Yeast_v4_1_Rep_2_S2_L002_R1_001.fastq.gz -r AYR_Yeast_v4_1_Rep_2_S2_L002_R2_001.fastq.gz -o AYR_WT_rep2-Ds-R1R2mrgd &
pear -f AYR_Yeast_v4_1_Rep_3_S3_L002_R1_001.fastq.gz -r AYR_Yeast_v4_1_Rep_3_S3_L002_R2_001.fastq.gz -o AYR_WT_rep3-Ds-R1R2mrgd &
pear -f BRY_E_coli_v1_S12_L002_R1_001.fastq.gz -r BRY_E_coli_v1_S12_L002_R2_001.fastq.gz -o BYR__Ec-Ds-R1R2mrgd &
pear -f BYR_Yeast_v1_1_Rep_1_S9_L002_R1_001.fastq.gz -r BYR_Yeast_v1_1_Rep_1_S9_L002_R2_001.fastq.gz -o BYR_WT_rep1-Ds-R1R2mrgd &
pear -f BYR_Yeast_v1_1_Rep_2_S10_L002_R1_001.fastq.gz -r BYR_Yeast_v1_1_Rep_2_S10_L002_R2_001.fastq.gz -o BYR_WT_rep2-Ds-R1R2mrgd &
pear -f BYR_Yeast_v1_1_Rep_3_S11_L002_R1_001.fastq.gz -r BYR_Yeast_v1_1_Rep_3_S11_L002_R2_001.fastq.gz -o BYR_WT_rep3-Ds-R1R2mrgd &
wait

# process files in fastq_18179Kap_N18164_L002 folder
cd /bgfs/ckaplan/Yunye/3-TSS_sequence_library/1-WT-DNAseq-TxGen_18179Kap/fastq_18179Kap_N18164_L002/

pear -f ARY_E_coli_v1_S8_L002_R1_001.fastq.gz -r ARY_E_coli_v1_S8_L002_R2_001.fastq.gz -o ARY__Ec-Ds-R1R2mrgd &
pear -f ARY_Yeast_v1_1_Rep_1_S5_L002_R1_001.fastq.gz -r ARY_Yeast_v1_1_Rep_1_S5_L002_R2_001.fastq.gz -o ARY_WT_rep1-Ds-R1R2mrgd &
pear -f ARY_Yeast_v1_1_Rep_2_S6_L002_R1_001.fastq.gz -r ARY_Yeast_v1_1_Rep_2_S6_L002_R2_001.fastq.gz -o ARY_WT_rep2-Ds-R1R2mrgd &
pear -f ARY_Yeast_v1_1_Rep_3_S7_L002_R1_001.fastq.gz -r ARY_Yeast_v1_1_Rep_3_S7_L002_R2_001.fastq.gz -o ARY_WT_rep3-Ds-R1R2mrgd &
pear -f AYR_E_coli_v4_S4_L002_R1_001.fastq.gz -r AYR_E_coli_v4_S4_L002_R2_001.fastq.gz -o AYR__Ec-Ds-R1R2mrgd &
pear -f AYR_Yeast_v4_1_Rep_1_S1_L002_R1_001.fastq.gz -r AYR_Yeast_v4_1_Rep_1_S1_L002_R2_001.fastq.gz -o AYR_WT_rep1-Ds-R1R2mrgd &
pear -f AYR_Yeast_v4_1_Rep_2_S2_L002_R1_001.fastq.gz -r AYR_Yeast_v4_1_Rep_2_S2_L002_R2_001.fastq.gz -o AYR_WT_rep2-Ds-R1R2mrgd &
pear -f AYR_Yeast_v4_1_Rep_3_S3_L002_R1_001.fastq.gz -r AYR_Yeast_v4_1_Rep_3_S3_L002_R2_001.fastq.gz -o AYR_WT_rep3-Ds-R1R2mrgd &
pear -f BRY_E_coli_v1_S12_L002_R1_001.fastq.gz -r BRY_E_coli_v1_S12_L002_R2_001.fastq.gz -o BYR__Ec-Ds-R1R2mrgd &
pear -f BYR_Yeast_v1_1_Rep_1_S9_L002_R1_001.fastq.gz -r BYR_Yeast_v1_1_Rep_1_S9_L002_R2_001.fastq.gz -o BYR_WT_rep1-Ds-R1R2mrgd &
pear -f BYR_Yeast_v1_1_Rep_2_S10_L002_R1_001.fastq.gz -r BYR_Yeast_v1_1_Rep_2_S10_L002_R2_001.fastq.gz -o BYR_WT_rep2-Ds-R1R2mrgd &
pear -f BYR_Yeast_v1_1_Rep_3_S11_L002_R1_001.fastq.gz -r BYR_Yeast_v1_1_Rep_3_S11_L002_R2_001.fastq.gz -o BYR_WT_rep3-Ds-R1R2mrgd &
wait

echo "yeah! I am finished!"