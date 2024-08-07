#!/bin/bash
#SBATCH --job-name=6-UMItools-TssCXed-Bar_CX
#SBATCH --output=6-UMItools-TssCXed-Bar_CX.out
#SBATCH -t 3-00:00
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=xxxxxx@pitt.edu
#SBATCH --mail-type=BEGIN,END,FAIL

### To correct Barcode_v2 sequence based on UMI-tools corrected TSS sequence using UMI-tools
# remotely use Pitt HTC cluster

# Inputs: mock SAM files [*-TssCXed-mockBarCXonTss.sam] generated by <5-MkMockSam-TssCXed-BarOnTss.py>

# Outputs: standard output of UMI-tools, file prefix naming: library_mutant_replicate-Ds-TssCXed-mockBarCXonTss-sorted

module purge
module load gcc/8.2.0
module load umi-tools/1.0.0
module load samtools/1.10

cd /bgfs/ckaplan/Yunye/3-TSS_sequence_library/1-WT-DNAseq-TxGen_18179Kap/5_6-TssCXed-BarOnTss/

# 1. Convert SAM file to BAM file using samtools
samtools view -bS AYR__Ec-Ds-TssCXed-mockBarCXonTss.sam > AYR__Ec-Ds-TssCXed-mockBarCXonTss.bam &
samtools view -bS BYR__Ec-Ds-TssCXed-mockBarCXonTss.sam > BYR__Ec-Ds-TssCXed-mockBarCXonTss.bam &
samtools view -bS ARY__Ec-Ds-TssCXed-mockBarCXonTss.sam > ARY__Ec-Ds-TssCXed-mockBarCXonTss.bam &
samtools view -bS AYR_WT_rep1-Ds-TssCXed-mockBarCXonTss.sam > AYR_WT_rep1-Ds-TssCXed-mockBarCXonTss.bam &
samtools view -bS AYR_WT_rep2-Ds-TssCXed-mockBarCXonTss.sam > AYR_WT_rep2-Ds-TssCXed-mockBarCXonTss.bam &
samtools view -bS AYR_WT_rep3-Ds-TssCXed-mockBarCXonTss.sam > AYR_WT_rep3-Ds-TssCXed-mockBarCXonTss.bam &
samtools view -bS BYR_WT_rep1-Ds-TssCXed-mockBarCXonTss.sam > BYR_WT_rep1-Ds-TssCXed-mockBarCXonTss.bam &
samtools view -bS BYR_WT_rep2-Ds-TssCXed-mockBarCXonTss.sam > BYR_WT_rep2-Ds-TssCXed-mockBarCXonTss.bam &
samtools view -bS BYR_WT_rep3-Ds-TssCXed-mockBarCXonTss.sam > BYR_WT_rep3-Ds-TssCXed-mockBarCXonTss.bam &
samtools view -bS ARY_WT_rep1-Ds-TssCXed-mockBarCXonTss.sam > ARY_WT_rep1-Ds-TssCXed-mockBarCXonTss.bam &
samtools view -bS ARY_WT_rep2-Ds-TssCXed-mockBarCXonTss.sam > ARY_WT_rep2-Ds-TssCXed-mockBarCXonTss.bam &
samtools view -bS ARY_WT_rep3-Ds-TssCXed-mockBarCXonTss.sam > ARY_WT_rep3-Ds-TssCXed-mockBarCXonTss.bam &
wait

# 2. sort BAM file
samtools sort AYR__Ec-Ds-TssCXed-mockBarCXonTss.bam -o AYR__Ec-Ds-TssCXed-mockBarCXonTss-sorted.bam &
samtools sort BYR__Ec-Ds-TssCXed-mockBarCXonTss.bam -o BYR__Ec-Ds-TssCXed-mockBarCXonTss-sorted.bam &
samtools sort ARY__Ec-Ds-TssCXed-mockBarCXonTss.bam -o ARY__Ec-Ds-TssCXed-mockBarCXonTss-sorted.bam &
samtools sort AYR_WT_rep1-Ds-TssCXed-mockBarCXonTss.bam -o AYR_WT_rep1-Ds-TssCXed-mockBarCXonTss-sorted.bam &
samtools sort AYR_WT_rep2-Ds-TssCXed-mockBarCXonTss.bam -o AYR_WT_rep2-Ds-TssCXed-mockBarCXonTss-sorted.bam &
samtools sort AYR_WT_rep3-Ds-TssCXed-mockBarCXonTss.bam -o AYR_WT_rep3-Ds-TssCXed-mockBarCXonTss-sorted.bam &
samtools sort BYR_WT_rep1-Ds-TssCXed-mockBarCXonTss.bam -o BYR_WT_rep1-Ds-TssCXed-mockBarCXonTss-sorted.bam &
samtools sort BYR_WT_rep2-Ds-TssCXed-mockBarCXonTss.bam -o BYR_WT_rep2-Ds-TssCXed-mockBarCXonTss-sorted.bam &
samtools sort BYR_WT_rep3-Ds-TssCXed-mockBarCXonTss.bam -o BYR_WT_rep3-Ds-TssCXed-mockBarCXonTss-sorted.bam &
samtools sort ARY_WT_rep1-Ds-TssCXed-mockBarCXonTss.bam -o ARY_WT_rep1-Ds-TssCXed-mockBarCXonTss-sorted.bam &
samtools sort ARY_WT_rep2-Ds-TssCXed-mockBarCXonTss.bam -o ARY_WT_rep2-Ds-TssCXed-mockBarCXonTss-sorted.bam &
samtools sort ARY_WT_rep3-Ds-TssCXed-mockBarCXonTss.bam -o ARY_WT_rep3-Ds-TssCXed-mockBarCXonTss-sorted.bam &
wait

# 3. index sorted BAM file
samtools index AYR__Ec-Ds-TssCXed-mockBarCXonTss-sorted.bam &
samtools index BYR__Ec-Ds-TssCXed-mockBarCXonTss-sorted.bam &
samtools index ARY__Ec-Ds-TssCXed-mockBarCXonTss-sorted.bam &
samtools index AYR_WT_rep1-Ds-TssCXed-mockBarCXonTss-sorted.bam &
samtools index AYR_WT_rep2-Ds-TssCXed-mockBarCXonTss-sorted.bam &
samtools index AYR_WT_rep3-Ds-TssCXed-mockBarCXonTss-sorted.bam &
samtools index BYR_WT_rep1-Ds-TssCXed-mockBarCXonTss-sorted.bam &
samtools index BYR_WT_rep2-Ds-TssCXed-mockBarCXonTss-sorted.bam &
samtools index BYR_WT_rep3-Ds-TssCXed-mockBarCXonTss-sorted.bam &
samtools index ARY_WT_rep1-Ds-TssCXed-mockBarCXonTss-sorted.bam &
samtools index ARY_WT_rep2-Ds-TssCXed-mockBarCXonTss-sorted.bam &
samtools index ARY_WT_rep3-Ds-TssCXed-mockBarCXonTss-sorted.bam &
wait

# 4. group = error correct
umi_tools group -I AYR__Ec-Ds-TssCXed-mockBarCXonTss-sorted.bam --group-out=AYR__Ec-Ds-TssCXed-mockBarCXonTss-sorted-grouped.tsv --log=AYR__Ec-Ds-TssCXed-mockBarCXonTss-sorted-group.log &
umi_tools group -I BYR__Ec-Ds-TssCXed-mockBarCXonTss-sorted.bam --group-out=BYR__Ec-Ds-TssCXed-mockBarCXonTss-sorted-grouped.tsv --log=BYR__Ec-Ds-TssCXed-mockBarCXonTss-sorted-group.log &
umi_tools group -I ARY__Ec-Ds-TssCXed-mockBarCXonTss-sorted.bam --group-out=ARY__Ec-Ds-TssCXed-mockBarCXonTss-sorted-grouped.tsv --log=ARY__Ec-Ds-TssCXed-mockBarCXonTss-sorted-group.log &
umi_tools group -I AYR_WT_rep1-Ds-TssCXed-mockBarCXonTss-sorted.bam --group-out=AYR_WT_rep1-Ds-TssCXed-mockBarCXonTss-sorted-grouped.tsv --log=AYR_WT_rep1-Ds-TssCXed-mockBarCXonTss-sorted-group.log &
umi_tools group -I AYR_WT_rep2-Ds-TssCXed-mockBarCXonTss-sorted.bam --group-out=AYR_WT_rep2-Ds-TssCXed-mockBarCXonTss-sorted-grouped.tsv --log=AYR_WT_rep2-Ds-TssCXed-mockBarCXonTss-sorted-group.log &
umi_tools group -I AYR_WT_rep3-Ds-TssCXed-mockBarCXonTss-sorted.bam --group-out=AYR_WT_rep3-Ds-TssCXed-mockBarCXonTss-sorted-grouped.tsv --log=AYR_WT_rep3-Ds-TssCXed-mockBarCXonTss-sorted-group.log &
umi_tools group -I BYR_WT_rep1-Ds-TssCXed-mockBarCXonTss-sorted.bam --group-out=BYR_WT_rep1-Ds-TssCXed-mockBarCXonTss-sorted-grouped.tsv --log=BYR_WT_rep1-Ds-TssCXed-mockBarCXonTss-sorted-group.log &
umi_tools group -I BYR_WT_rep2-Ds-TssCXed-mockBarCXonTss-sorted.bam --group-out=BYR_WT_rep2-Ds-TssCXed-mockBarCXonTss-sorted-grouped.tsv --log=BYR_WT_rep2-Ds-TssCXed-mockBarCXonTss-sorted-group.log &
umi_tools group -I BYR_WT_rep3-Ds-TssCXed-mockBarCXonTss-sorted.bam --group-out=BYR_WT_rep3-Ds-TssCXed-mockBarCXonTss-sorted-grouped.tsv --log=BYR_WT_rep3-Ds-TssCXed-mockBarCXonTss-sorted-group.log &
umi_tools group -I ARY_WT_rep1-Ds-TssCXed-mockBarCXonTss-sorted.bam --group-out=ARY_WT_rep1-Ds-TssCXed-mockBarCXonTss-sorted-grouped.tsv --log=ARY_WT_rep1-Ds-TssCXed-mockBarCXonTss-sorted-group.log &
umi_tools group -I ARY_WT_rep2-Ds-TssCXed-mockBarCXonTss-sorted.bam --group-out=ARY_WT_rep2-Ds-TssCXed-mockBarCXonTss-sorted-grouped.tsv --log=ARY_WT_rep2-Ds-TssCXed-mockBarCXonTss-sorted-group.log &
umi_tools group -I ARY_WT_rep3-Ds-TssCXed-mockBarCXonTss-sorted.bam --group-out=ARY_WT_rep3-Ds-TssCXed-mockBarCXonTss-sorted-grouped.tsv --log=ARY_WT_rep3-Ds-TssCXed-mockBarCXonTss-sorted-group.log &
wait

date
echo "yeah! I am finished!"