# remotely use Pitt HTC cluster
# R= 4.0.0

### 1. To train the final model with selected robust features, save parameters, and visualization
### 2. To perform PCA analysis on wt & mut models, and visualization
### 3. To predict TSS efficiency of positions in genome (within known promoter windows), and visualizaiton

setwd("/bgfs/ckaplan/Yunye/3-TSS_sequence_library/8-modeling-ct4")

library(LSD) # for heatscatter()
library(broom) # for tidy()
library(dplyr) # for separate(), %>%, bind_rows()
library(tidyr) # for separate(), %>%
library(factoextra) # for PCA visualization, fviz_eig(), fviz_pca_biplot()

##### 0. functions #####
# 0.1 function to evalidate model
eval_results <- function(actual, predicted, df){
  SSE <- sum((actual - predicted)^2)
  SST <- sum((actual - mean(actual))^2) # may need to use mean(training) as a baseline model?
  R_square <- 1 - SSE/SST
  RMSE = sqrt(SSE/nrow(df))
  
  data.frame(
    RMSE = RMSE,
    Rsquare = R_square
  )
}
# 0.2 function for visualization
heatscatter_vis <- function(x, y, main, xlab, ylab, xlim=NULL, ylim=NULL,cex.main=2,
                            tick=TRUE,
                            diagonal=TRUE, vline=FALSE, hline=FALSE){
  par(mar=c(4.5, 5, 3, 1), lwd=4, pty="s")
  heatscatter(x=x, y=y,
              colpal="bl2gr2rd",
              xlim=xlim, ylim=ylim, xaxt="n", yaxt="n",
              main="", xlab=xlab, ylab=ylab,
              cex.lab=2.5, cex.axis=2)
  title(main=main, cex.main=cex.main, line=0.5)
  axis(side=1, tick=tick, cex.axis = 2)
  axis(side=2, tick=tick, cex.axis = 2)
  if(diagonal){abline(coef = c(0,1), col="grey", lwd=3, lty=2)}
  if(vline){abline(v= vline, col= "red", lty = 2, lwd= 3)}
  if(hline){abline(h= hline, col= "red", lty = 2, lwd= 3)}
}
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

##### 1. dataset setup #####
### 1.1. get all datasets including (WT+all mutatns) x (-8 to +4 TSSs)
mtx_folder = '/bgfs/ckaplan/Yunye/3-TSS_sequence_library/7-poolDs_cmbTs/'
master_mtx <- read.csv(paste0(mtx_folder,'all-pDs_B12RepsCmb-ct4-eff-84_seqn11p9_info_mtx.csv'),
                       colClasses=c("numeric", rep("factor",24)))
  # 25 columns: [eff] - [pos-11to+9]*20 - [PolII] - [lib] - [TSS] - [TSSvar]
master_mtx$o_mtf811 <- ifelse(grepl("A[ACTG]{6}[CT][AG]",master_mtx$TSSvar, perl=TRUE), "AYR",
                              ifelse(grepl("T[ACTG]{6}[CT][AG]",master_mtx$TSSvar, perl=TRUE), "TYR",
                                     ifelse(grepl("C[ACTG]{6}[CT][AG]",master_mtx$TSSvar, perl=TRUE), "CYR",
                                            ifelse(grepl("G[ACTG]{6}[CT][AG]",master_mtx$TSSvar, perl=TRUE), "GYR",
                                                   ifelse(grepl("A[ACTG]{6}[AG][CT]",master_mtx$TSSvar, perl=TRUE), "ARY",
                                                          NA    )))))
master_mtx <- master_mtx[!(is.na(master_mtx$eff)),] # only keep rows w/o NA for eff

### 1.2. use datasets of (-8 to +1, +4 TSSs)x(AYR+BYR+ARY)+ (+2 TSS)x(ARY)
master_mtx <- subset(master_mtx, TSS != '3') # 4322382      26
master_mtx <- master_mtx[!(grepl("[ACTG]YR",master_mtx$o_mtf811, perl=TRUE) & master_mtx$TSS == '2'), ] # 4010354

pol = 'WT'
sub_mtx <- subset(master_mtx, PolII == pol)

### 1.3. training 80%, test 20%
set.seed(42) # sets the random seed for reproducibility of results during pipline generation
perc <- 80/100

index = sample(1:nrow(sub_mtx), perc*nrow(sub_mtx))
train <- sub_mtx[index,]
test <- sub_mtx[-index,]

train$eff01 <- train$eff/100 # transform eff from [0,100] to [0,1]
test$eff01 <- test$eff/100 # transform eff from [0,100] to [0,1]

dim(sub_mtx)
dim(train)
dim(test)
summary(sub_mtx)
rm(perc,index)

##### 2. model training & prediction on test set & parameters saving #####
### 2.1 train with 9 additive + 1 itr ----
lg <- glm(eff01 ~ 1+ n9+n8+n7+n4+n3+n2+n1+p1+p2+ n9*n8,
          data = train, family = "quasibinomial",
          y = FALSE, # prevent the response vector from being returned
          model = FALSE) # prevent the model.frame from being returned

## save models
#pol = 'WT'
LR_para = 'lg9i1' # logistic regression, 9 additive terms, 1 interaction
models_folder = "/bgfs/ckaplan/Yunye/3-TSS_sequence_library/8-modeling-ct4/2-LR_models/"
saveRDS(lg, paste0(models_folder,"LR_",LR_para,"-",pol,".RData")) # save model

## load models
pol = 'WT'
LR_para = 'lg9i1'
models_folder = "/bgfs/ckaplan/Yunye/3-TSS_sequence_library/8-modeling-ct4/2-LR_models/"
lg <- readRDS(paste0(models_folder,"LR_",LR_para,"-",pol,".RData"))

ll.null <- lg$null.deviance/(-2)
ll.proposed <- lg$deviance/(-2)
(ll.null - ll.proposed) / ll.null

summary(lg)
rm(ll.null,ll.proposed)

### 2.2 predict and evaluate the model on test data ----
test$pred01 <- predict(lg, newdata = test, type="response")
test$pred <- test$pred01*100

eval_results(test$eff01, test$pred01, test)
cor(test$eff01, test$pred01, method="pearson")
cor(test$eff01, test$pred01, method="spearman")

### 2.3 save parameters for trained model into dataframe / interaction matrix ----
output_folder = '/bgfs/ckaplan/Yunye/3-TSS_sequence_library/8-modeling-ct4/2-LR_lg9i1-para/'
## 2.3.1 save additives only (w/o interactions);
data.frame(tidy(lg)[c("term","estimate")]) %>%
  subset((term != "(Intercept)") & (!grepl(":",term))) %>%
  separate(col=term, into=c("pos","base"), sep=c(-1), remove=TRUE) %>%
  spread(base,estimate) %>%
  write.csv(paste0(output_folder,"LR_",LR_para,"-",pol,"-para_idv.csv"), row.names = FALSE)

## 2.3.2 save interaction (n9/n8) only
data.frame(tidy(lg)[c("term","estimate")]) %>%
  subset(grepl(":n8", term)) %>%
  separate(col=term, into=c("itr1_base","itr2_base"), sep=":", remove=TRUE) %>%
  spread(key=itr2_base, value=estimate) %>%
  write.csv(paste0(output_folder,"LR_",LR_para,"-",pol,"-para_n9n8itr.csv"), row.names = FALSE)

## 2.3.3 save itr2(n8) under different itr1(n9) (=centered (itr2 coef+interaction))
itr <- data.frame(tidy(lg)[c("term","estimate")]) %>%
  subset(grepl(":n8", term)) %>%
  separate(col=term, into=c("itr1_base","itr2_base"), sep=":", remove=TRUE) %>%
  spread(key=itr2_base, value=estimate)
rownames(itr) <- itr$itr1_base
itr$itr1_base <- NULL
itr$n8A <- 0
itr["n9A",] <- 0

n8_coef <- data.frame(tidy(lg)[c("term","estimate")]) %>%
  subset((grepl("n8",term)) & (!grepl(":",term))) %>%
  separate(col=term, into=c("pos","base"), sep=c(-1), remove=TRUE) %>%
  spread(base,estimate)
rownames(n8_coef) <- n8_coef$pos # set_index
n8_coef$pos <- NULL # remove pos column
n8_coef[,"A"] = 0

n9_itr_n8 = itr[]
for (base2 in c('A','G','C','T')){
  for (base1 in c('A','G','C','T')){
    n9_itr_n8[paste0("n9",base1),paste0("n8",base2)] = 
      itr[paste0("n9",base1),paste0("n8",base2)] + n8_coef["n8",base2]
  }
}
cen_df(n9_itr_n8, how = 'row')[paste0("n9",c('A','G','C','T')), paste0("n8",c('A','G','C','T'))] %>%
  write.csv(paste0(output_folder,"LR_",LR_para,"-",pol,"-n8para_under_n9.csv"), row.names = TRUE)

rm(itr, n8_coef, n9_itr_n8,base1,base2)


##### 3. visualization #####
### 3.1 measurement vs prediction ----
output_folder = '/bgfs/ckaplan/Yunye/3-TSS_sequence_library/8-modeling-ct4/2-LR_lg9i1-figs/'

# test set
pdf(file=paste0(output_folder,pol,"-test.pdf")) # vector image
heatscatter_vis(x = test$pred, y = test$eff,
                xlim= c(0,100), ylim=c(0,100),
                main=paste0(pol," - test set"), xlab="prediction",ylab="measurement",
                tick=TRUE,diagonal=TRUE, vline=FALSE, hline=FALSE)
dev.off()

# test set, w/ 5% eff cut-off
pdf(file=paste0(output_folder,pol,"-test_5cnt.pdf")) # vector image
heatscatter_vis(x = subset(test, eff>5 & pred>5)$pred, y = subset(test, eff>5 & pred>5)$eff,
                xlim= c(0,100), ylim=c(0,100),
                main=paste0(pol," - test set (cut-off = 5)"), xlab="prediction",ylab="measurement",
                tick=TRUE,diagonal=TRUE, vline=FALSE, hline=FALSE)
dev.off()

#### 4. PCA analysis on wt & mut models ####
### 4.1 get all parameters of all WT/mut models
LR_para = 'lg9i1'
models_folder = "/bgfs/ckaplan/Yunye/3-TSS_sequence_library/8-modeling-ct4/2-LR_models/"
para_folder = '/bgfs/ckaplan/Yunye/3-TSS_sequence_library/8-modeling-ct4/2-LR_lg9i1-para/'

for(pol in c('G1097D','E1103G','WT','F1086S','H1085Q')){
  temp <- readRDS(paste0(models_folder,"LR_",LR_para,"-",pol,".RData")) %>%
    tidy() %>%
    data.frame() %>%
    subset(select=c('term','estimate'))
  colnames(temp)[ncol(temp)] <- pol
  
  if (exists("PCA_mtx")){
    PCA_mtx <- merge(PCA_mtx,temp, by = 'term', all=TRUE)
  } else {
    PCA_mtx <- temp[]
  }
}
write.csv(PCA_mtx, paste0(para_folder,"LR_",LR_para,"-AllPol_para.csv"), row.names = FALSE)
rm(temp)

rownames(PCA_mtx) <- PCA_mtx$term
PCA_mtx_t <- PCA_mtx[, c('G1097D','E1103G','WT','F1086S','H1085Q')] %>% t()

### 4.2 PCA analaysis & visualization
pca.out = prcomp(PCA_mtx_t, scale = FALSE, center = TRUE)
summary(pca.out)

# scree plot of eigenvalues, showing the percentage of variances explained by each PC
fviz_eig(pca.out)
# Biplot of individuals and variables; save as 6.5 x 6 in
pdf(file=paste0(para_folder,"LR_",LR_para,"-Allmut_wt_para-PCA_15.pdf"))
fviz_pca_biplot(pca.out,
                pointsize = 3,
                col.var = "contrib", # Color by contributions to the PC
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                select.var = list(contrib = 15),
                repel = TRUE # Avoid text overlapping
                )# + xlim(-0.4, 0.4) + ylim (-0.4, 0.4)
dev.off()

rm(PCA_mtx,PCA_mtx_t,pca.out)

#### 5. prediction on genomic data #####
### 5.1 get genomic dataset in 401*5979 windows ----
pol = 'WT'
grad = 'TZ'
mrdg_indv = 'mrgd'
smp = paste0(grad,'_',pol,'_',mrdg_indv)

mtx_folder = "/bgfs/ckaplan/Yunye/3-TSS_sequence_library/9-gnm/4-AllPos_mtx/"
vivo <- read.csv(paste0(mtx_folder,smp,'-AllPos-eff_seqn11p9_anno-mtrx.csv'),
                 colClasses=c("numeric",rep("factor",22),"numeric","factor","numeric",rep("factor",3),rep("numeric",2)))
  # 31 columns: [eff] - [pos-11to+9]*20 - [annotations]*10
  # "eff", "n11"-"p9", "PolII","PWpos","distToMed","TSSn11p9","INDX6044","gene","TATAclass","TAF1status","ttlRNAcnt","RNAcnt"
  # 401*5979 = 2,397,579

vivo <- vivo[!(is.na(vivo$eff)),] %>% # only keep rows w/o NA for eff
  subset(n11 != '' & p9 !='') # filter out a few of missing data caused by 1) length of some promoter windows <401, or 2) promoters at the ends of chromosomes

### 5.2 LR prediction ----
LR_para = 'lg9i1'
models_folder = "/bgfs/ckaplan/Yunye/3-TSS_sequence_library/8-modeling-ct4/2-LR_models/"
lg <- readRDS(paste0(models_folder,"LR_",LR_para,"-",pol,".RData")) # load trained model
ll.null <- lg$null.deviance/(-2)
ll.proposed <- lg$deviance/(-2)
(ll.null - ll.proposed) / ll.null
rm(ll.null, ll.proposed)

vivo$eff01 <- vivo$eff/100 # transform eff from [0,100] to [0,1]

vivo$pred01 <- predict(lg, newdata = vivo, type="response")
vivo$pred <- vivo$pred01*100

eval_results(vivo$eff01, vivo$pred01, vivo)
cor(vivo$eff01, vivo$pred01, method="pearson")
cor(vivo$eff01, vivo$pred01, method="spearman")

### 5.3 visualization: measurement vs prediction ----
output_folder = "/bgfs/ckaplan/Yunye/3-TSS_sequence_library/8-modeling-ct4/2-LR_lg9i1-gnm/"

## 5.3.1 all positions within promoter windows or median TSSs only ----
# all positions
pdf(file=paste0(output_folder, smp,"-",LR_para,".pdf"))
heatscatter_vis(x = vivo$pred, y = vivo$eff,
                xlim=c(0,100), ylim=c(0,100),
                main=paste0(grad,'_',pol),
                xlab="prediction", ylab="measurement")
dev.off()
cor(vivo$pred, vivo$eff, method="pearson")

# all positions, w/ 5% eff cut-off
pdf(file=paste0(output_folder, smp,"-",LR_para,"_eff5cut.pdf"))
heatscatter_vis(x = subset(vivo, eff>5 & pred>5)$pred, y = subset(vivo, eff>5 & pred>5)$eff,
                xlim=c(0,100), ylim=c(0,100),
                main= paste0(grad,'_',pol," (cut-off = 5)"),
                xlab= "prediction", ylab= "measurement")
dev.off()
cor(subset(vivo, eff>5 & pred>5)$pred, subset(vivo, eff>5 & pred>5)$eff, method="pearson")

# median TSSs only
pdf(file=paste0(output_folder, smp,"-",LR_para,"-medTSS.pdf"))
heatscatter_vis(x= subset(vivo, distToMed == 0)$pred, 
                y= subset(vivo, distToMed == 0)$eff,
                xlim=c(0,100), ylim=c(0,100),
                main= paste0(grad,'_',pol," - medTSS"),
                xlab= "prediction", ylab= "measurement")
dev.off()
cor(subset(vivo, distToMed == 0)$pred, subset(vivo, distToMed == 0)$eff, method="pearson")

## 5.3.2 median TSSs within subgroups based on TAF1-status and exp levels ----
# output figures for all combinations
TAF1_info <- data.frame(abbr = c('TAF1dep', 'TAF1enr', 'TAF1unk'),
                        TAF1status = c("TAF1-depleted", "TAF1-enriched", "-"))
exp_info <- data.frame(abbr = c('low', 'med', 'high'),
                       HJ_WT_ttlRNA_l = c(0,200,1000),
                       HJ_WT_ttlRNA_h = c(200,1000,3400000),
                       TZ_WT_ttlRNA_l = c(0,200,1000),
                       TZ_WT_ttlRNA_h = c(200,1000,2200000))
for(i in c(1,2,3)){
  TAF1 <- as.character(TAF1_info[i,"TAF1status"])
  for(j in c(1,2,3)){
    ttlRNA_l = exp_info[j, paste0(grad,"_",pol,"_ttlRNA_l")]
    ttlRNA_h = exp_info[j, paste0(grad,"_",pol,"_ttlRNA_h")]
    vis_temp <- subset(vivo, 
                       distToMed == 0
                       & TAF1status == TAF1
                       & ttlRNAcnt >= ttlRNA_l
                       & ttlRNAcnt < ttlRNA_h 
    )
    pdf(file=paste0(output_folder, smp,"-",LR_para,"-medTSS_",as.character(TAF1_info[i,"abbr"]),"_",as.character(exp_info[j,"abbr"]),".pdf"))
    heatscatter_vis(x = vis_temp$pred, y = vis_temp$eff,
                    xlim=c(0,100), ylim=c(0,100),
                    main=paste(grad,'_',pol," - medTSS\n",TAF1," & ttlRNAcnt=[",ttlRNA_l,",",ttlRNA_h,")", sep=""),
                    xlab="prediction", ylab="measurement",
                    cex.main=1.3)
    dev.off()
    
    print(paste0(as.character(TAF1_info[i,"abbr"]),"-",as.character(exp_info[j,"abbr"]),": ",
                 "N= ",dim(vis_temp)[1], "; pearson-r = ", cor(vis_temp$pred, vis_temp$eff,method="pearson")))
  }
}