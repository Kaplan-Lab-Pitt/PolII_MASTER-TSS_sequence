# remotely use Pitt HTC cluster
# R= 4.0.0
setwd("/bgfs/ckaplan/Yunye/3-TSS_sequence_library/9-gnm")

library(dplyr) # for separate(), %>%, bind_rows(), top_n()
library(tidyr) # for separate(), %>%

#### 0. functions ####
# 0.3 function for centralization of columns "col"/rows "row"/dataframe "df"
cen_df <- function(df, how="df"){
  if (how == 'df'){
    df <- cen_df(df, how='row')
    df_cen <- cen_df(df, how='col')
  } else {
    df_cen <- df[]
    for (r in seq(nrow(df))){
      for (c in seq(ncol(df))){
        if (how == 'row'){
          df_cen[r, c] = df[r, c] - mean(as.matrix(df[r,]))
        } else if (how == 'col'){
          df_cen[r, c] = df[r, c] - mean(as.matrix(df[,c]))
        }
      }
    }
  }
  return(df_cen)
}

#### 1. dataset setup ####
### get all data in 401*5979 windows
pol = 'WT'
grad = 'TZ'
mrdg_indv = 'mrgd'
smp = paste0(grad,'_',pol,'_',mrdg_indv)

mtx_folder = "/bgfs/ckaplan/Yunye/3-TSS_sequence_library/9-gnm/4-AllPos_mtx/"
master_mtx <- read.csv(paste0(mtx_folder,smp,'-AllPos-eff_seqn11p9_anno-mtrx.csv'),
                       colClasses=c("numeric",rep("factor",22),"numeric","factor","numeric",rep("factor",3),rep("numeric",2)))
  # 31 columns: [eff] - [pos-11to+9]*20 - [annotations]*10
  # "eff", "n11"-"p9", "PolII","PWpos","distToMed","TSSn11p9","INDX6044","gene","TATAclass","TAF1status","ttlRNAcnt","RNAcnt"
  # 401*5979 = 2,397,579

master_mtx <- master_mtx[!(is.na(master_mtx$eff)),] # only keep rows w/o NA for eff
master_mtx <- subset(master_mtx, n11 != '' & p9 !='') # filter out a few of missing data caused by 1) length of some promoter windows <401, or 2) promoters at the ends of chromosomes
for (p in seq(2,21,1)){
  master_mtx[,p] <- factor(master_mtx[,p], levels=c("A", "G", "C", "T")) # to re-ORDER levels
}

dim(master_mtx)
rm(p)

#### 6. median-based analysis of subgroups #####
### "relative efficiency" at positions -11 to +9 of "median TSS"
n11p9_poss <- c('n11','n10','n9','n8','n7','n6','n5','n4','n3','n2','n1','p1','p2','p3','p4','p5','p6','p7','p8','p9')

sub_mtx <- subset(master_mtx, PWpos == 'medTSS')

df <- data.frame(base=c('A','G','C','T'))
for (pos in n11p9_poss){
  temp = aggregate(eff ~ eval(as.symbol(pos)), data=sub_mtx, FUN=median)
  df <- merge(df,temp,by.x = 'base', by.y = 'eval(as.symbol(pos))', all=TRUE)
  colnames(df)[ncol(df)] = pos
}
rownames(df) <- df$base
out_folder = "/bgfs/ckaplan/Yunye/3-TSS_sequence_library/9-gnm/5-gnm/"
cen_df(subset(df, select = -base), how="col") %>%
  t() %>%
  subset(select=c('A','G','C','T')) %>% # re-order columns
  write.csv(paste0(out_folder,smp,"-medTSS-median_centered.csv"), row.names = TRUE)

rm(df,temp,pos, sub_mtx)  

#### 8. output TSSvar for weblogo #####
out_folder = out_folder = "/bgfs/ckaplan/Yunye/3-TSS_sequence_library/9-gnm/5-gnm/"

# top d*10 % expressed Median TSSs, based on TSS-seq reads count for median TSSs
sub_mtx <- subset(master_mtx, PWpos == 'medTSS')
d <- 2 # decile
write.table(sub_mtx[(sub_mtx$RNAcnt > quantile(sub_mtx$RNAcnt, prob=1-d*10/100)),]["TSSn11p9"],
            paste0(out_folder,smp,"-top",d*10,".csv"), row.names = FALSE, col.names = FALSE)
# subgroup: -9/-8: A vs non-A
write.table(subset(sub_mtx[(sub_mtx$RNAcnt > quantile(sub_mtx$RNAcnt, prob=1-d*10/100)),], n9 != 'A')["TSSn11p9"],
            paste0(out_folder,smp,"-top",d*10,"-n9B",".csv"), row.names = FALSE, col.names = FALSE)

rm(d)
