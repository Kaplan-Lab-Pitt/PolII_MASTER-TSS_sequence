#!/bin/bash
#SBATCH --job-name=1-decompress
#SBATCH --output=1-decompress.out
#SBATCH -t 3-00:00
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=xxxxxx@pitt.edu
#SBATCH --mail-type=BEGIN,END,FAIL

### to decompress .fastq.gz files, with original .gz kept
# remotely use Pitt HTC cluster
# gzip=1.5;

# 1-WT-TSSseq-Vv1399_1407:
cd /bgfs/ckaplan/Yunye/3-TSS_sequence_library/1-WT-TSSseq-Vv1399_1407/
gzip -d -c VV1399_S3_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/VV1399_S3_R1_001.fastq &
gzip -d -c VV1400_S4_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/VV1400_S4_R1_001.fastq &
gzip -d -c VV1401_S5_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/VV1401_S5_R1_001.fastq &
gzip -d -c VV1402_S6_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/VV1402_S6_R1_001.fastq &
gzip -d -c VV1403_S7_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/VV1403_S7_R1_001.fastq &
gzip -d -c VV1404_S8_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/VV1404_S8_R1_001.fastq &
gzip -d -c VV1405_S9_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/VV1405_S9_R1_001.fastq &
gzip -d -c VV1406_S10_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/VV1406_S10_R1_001.fastq &
gzip -d -c VV1407_S11_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/VV1407_S11_R1_001.fastq &
wait

# 2-E1103G-TSSseq-Vv1477_1485:
cd /bgfs/ckaplan/Yunye/3-TSS_sequence_library/2-E1103G-TSSseq-Vv1477_1485/
gzip -d -c VV1477_S4_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/VV1477_S4_R1_001.fastq &
gzip -d -c VV1478_S5_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/VV1478_S5_R1_001.fastq &
gzip -d -c VV1479_S6_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/VV1479_S6_R1_001.fastq &
gzip -d -c VV1480_S7_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/VV1480_S7_R1_001.fastq &
gzip -d -c VV1481_S8_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/VV1481_S8_R1_001.fastq &
gzip -d -c VV1482_S9_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/VV1482_S9_R1_001.fastq &
gzip -d -c VV1483_S10_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/VV1483_S10_R1_001.fastq &
gzip -d -c VV1484_S11_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/VV1484_S11_R1_001.fastq &
gzip -d -c VV1485_S12_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/VV1485_S12_R1_001.fastq &
wait

# 3-F1086S-TSSseq-Vv1486_1494:
cd /bgfs/ckaplan/Yunye/3-TSS_sequence_library/3-F1086S-TSSseq-Vv1486_1494/
gzip -d -c VV1486_S4_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/VV1486_S4_R1_001.fastq &
gzip -d -c VV1487_S5_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/VV1487_S5_R1_001.fastq &
gzip -d -c VV1488_S6_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/VV1488_S6_R1_001.fastq &
gzip -d -c VV1489_S7_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/VV1489_S7_R1_001.fastq &
gzip -d -c VV1490_S8_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/VV1490_S8_R1_001.fastq &
gzip -d -c VV1491_S9_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/VV1491_S9_R1_001.fastq &
gzip -d -c VV1492_S10_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/VV1492_S10_R1_001.fastq &
gzip -d -c VV1493_S11_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/VV1493_S11_R1_001.fastq &
gzip -d -c VV1494_S12_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/VV1494_S12_R1_001.fastq &
wait

# 4-G1097D-TSSseq-Vv1589_1599:
cd /bgfs/ckaplan/Yunye/3-TSS_sequence_library/4-G1097D-TSSseq-Vv1589_1599/
gzip -d -c Vv1589_S1_L001_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1589_S1_L001_R1_001.fastq &
gzip -d -c Vv1589_S1_L002_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1589_S1_L002_R1_001.fastq &
gzip -d -c Vv1589_S1_L003_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1589_S1_L003_R1_001.fastq &
gzip -d -c Vv1589_S1_L004_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1589_S1_L004_R1_001.fastq &
gzip -d -c Vv1590_S2_L001_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1590_S2_L001_R1_001.fastq &
gzip -d -c Vv1590_S2_L002_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1590_S2_L002_R1_001.fastq &
gzip -d -c Vv1590_S2_L003_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1590_S2_L003_R1_001.fastq &
gzip -d -c Vv1590_S2_L004_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1590_S2_L004_R1_001.fastq &
gzip -d -c Vv1591_S3_L001_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1591_S3_L001_R1_001.fastq &
gzip -d -c Vv1591_S3_L002_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1591_S3_L002_R1_001.fastq &
gzip -d -c Vv1591_S3_L003_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1591_S3_L003_R1_001.fastq &
gzip -d -c Vv1591_S3_L004_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1591_S3_L004_R1_001.fastq &
gzip -d -c Vv1592_S4_L001_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1592_S4_L001_R1_001.fastq &
gzip -d -c Vv1592_S4_L002_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1592_S4_L002_R1_001.fastq &
gzip -d -c Vv1592_S4_L003_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1592_S4_L003_R1_001.fastq &
gzip -d -c Vv1592_S4_L004_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1592_S4_L004_R1_001.fastq &
gzip -d -c Vv1593_S5_L001_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1593_S5_L001_R1_001.fastq &
gzip -d -c Vv1593_S5_L002_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1593_S5_L002_R1_001.fastq &
gzip -d -c Vv1593_S5_L003_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1593_S5_L003_R1_001.fastq &
gzip -d -c Vv1593_S5_L004_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1593_S5_L004_R1_001.fastq &
gzip -d -c Vv1594_S6_L001_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1594_S6_L001_R1_001.fastq &
gzip -d -c Vv1594_S6_L002_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1594_S6_L002_R1_001.fastq &
gzip -d -c Vv1594_S6_L003_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1594_S6_L003_R1_001.fastq &
gzip -d -c Vv1594_S6_L004_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1594_S6_L004_R1_001.fastq &
gzip -d -c Vv1595_S7_L001_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1595_S7_L001_R1_001.fastq &
gzip -d -c Vv1595_S7_L002_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1595_S7_L002_R1_001.fastq &
gzip -d -c Vv1595_S7_L003_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1595_S7_L003_R1_001.fastq &
gzip -d -c Vv1595_S7_L004_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1595_S7_L004_R1_001.fastq &
gzip -d -c Vv1596_S8_L001_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1596_S8_L001_R1_001.fastq &
gzip -d -c Vv1596_S8_L002_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1596_S8_L002_R1_001.fastq &
gzip -d -c Vv1596_S8_L003_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1596_S8_L003_R1_001.fastq &
gzip -d -c Vv1596_S8_L004_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1596_S8_L004_R1_001.fastq &
gzip -d -c Vv1597_S9_L001_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1597_S9_L001_R1_001.fastq &
gzip -d -c Vv1597_S9_L002_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1597_S9_L002_R1_001.fastq &
gzip -d -c Vv1597_S9_L003_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1597_S9_L003_R1_001.fastq &
gzip -d -c Vv1597_S9_L004_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1597_S9_L004_R1_001.fastq &
gzip -d -c Vv1598_S10_L001_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1598_S10_L001_R1_001.fastq &
gzip -d -c Vv1598_S10_L002_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1598_S10_L002_R1_001.fastq &
gzip -d -c Vv1598_S10_L003_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1598_S10_L003_R1_001.fastq &
gzip -d -c Vv1598_S10_L004_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1598_S10_L004_R1_001.fastq &
gzip -d -c Vv1599_S11_L001_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1599_S11_L001_R1_001.fastq &
gzip -d -c Vv1599_S11_L002_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1599_S11_L002_R1_001.fastq &
gzip -d -c Vv1599_S11_L003_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1599_S11_L003_R1_001.fastq &
gzip -d -c Vv1599_S11_L004_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1599_S11_L004_R1_001.fastq &
wait

# 5-H1085Q-TSSseq-Vv1580_1588:
cd /bgfs/ckaplan/Yunye/3-TSS_sequence_library/5-H1085Q-TSSseq-Vv1580_1588/
gzip -d -c Vv1580_S2_L001_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1580_S2_L001_R1_001.fastq &
gzip -d -c Vv1580_S2_L002_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1580_S2_L002_R1_001.fastq &
gzip -d -c Vv1580_S2_L003_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1580_S2_L003_R1_001.fastq &
gzip -d -c Vv1580_S2_L004_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1580_S2_L004_R1_001.fastq &
gzip -d -c Vv1581_S3_L001_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1581_S3_L001_R1_001.fastq &
gzip -d -c Vv1581_S3_L002_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1581_S3_L002_R1_001.fastq &
gzip -d -c Vv1581_S3_L003_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1581_S3_L003_R1_001.fastq &
gzip -d -c Vv1581_S3_L004_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1581_S3_L004_R1_001.fastq &
gzip -d -c Vv1582_S4_L001_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1582_S4_L001_R1_001.fastq &
gzip -d -c Vv1582_S4_L002_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1582_S4_L002_R1_001.fastq &
gzip -d -c Vv1582_S4_L003_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1582_S4_L003_R1_001.fastq &
gzip -d -c Vv1582_S4_L004_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1582_S4_L004_R1_001.fastq &
gzip -d -c Vv1583_S5_L001_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1583_S5_L001_R1_001.fastq &
gzip -d -c Vv1583_S5_L002_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1583_S5_L002_R1_001.fastq &
gzip -d -c Vv1583_S5_L003_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1583_S5_L003_R1_001.fastq &
gzip -d -c Vv1583_S5_L004_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1583_S5_L004_R1_001.fastq &
gzip -d -c Vv1584_S6_L001_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1584_S6_L001_R1_001.fastq &
gzip -d -c Vv1584_S6_L002_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1584_S6_L002_R1_001.fastq &
gzip -d -c Vv1584_S6_L003_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1584_S6_L003_R1_001.fastq &
gzip -d -c Vv1584_S6_L004_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1584_S6_L004_R1_001.fastq &
gzip -d -c Vv1585_S7_L001_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1585_S7_L001_R1_001.fastq &
gzip -d -c Vv1585_S7_L002_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1585_S7_L002_R1_001.fastq &
gzip -d -c Vv1585_S7_L003_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1585_S7_L003_R1_001.fastq &
gzip -d -c Vv1585_S7_L004_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1585_S7_L004_R1_001.fastq &
gzip -d -c Vv1586_S8_L001_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1586_S8_L001_R1_001.fastq &
gzip -d -c Vv1586_S8_L002_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1586_S8_L002_R1_001.fastq &
gzip -d -c Vv1586_S8_L003_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1586_S8_L003_R1_001.fastq &
gzip -d -c Vv1586_S8_L004_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1586_S8_L004_R1_001.fastq &
gzip -d -c Vv1587_S9_L001_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1587_S9_L001_R1_001.fastq &
gzip -d -c Vv1587_S9_L002_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1587_S9_L002_R1_001.fastq &
gzip -d -c Vv1587_S9_L003_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1587_S9_L003_R1_001.fastq &
gzip -d -c Vv1587_S9_L004_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1587_S9_L004_R1_001.fastq &
gzip -d -c Vv1588_S10_L001_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1588_S10_L001_R1_001.fastq &
gzip -d -c Vv1588_S10_L002_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1588_S10_L002_R1_001.fastq &
gzip -d -c Vv1588_S10_L003_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1588_S10_L003_R1_001.fastq &
gzip -d -c Vv1588_S10_L004_R1_001.fastq.gz > /bgfs/ckaplan/Yunye/3-TSS_sequence_library/12345-Ts_b1-extr_CX_ddp/1-fastq/Vv1588_S10_L004_R1_001.fastq &
wait

date
echo "yeah! I am finished!"