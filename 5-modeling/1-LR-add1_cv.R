# remotely use Pitt HTC cluster
# R= 4.0.0

### To select robust features,
### employ a forward stepwise strategy with a 5-fold Cross-Validation (CV) in two major stages:
### for additive terms and for interactions

setwd("/bgfs/ckaplan/Yunye/3-TSS_sequence_library/8-modeling-ct4")

library(caret) # for train()

##### 0. functions #####
# 0.1 function to evalidate model
eval_results <- function(actual, predicted, df){
  SSE <- sum((actual - predicted)^2)
  SST <- sum((actual - mean(actual))^2)
  R_square <- 1 - SSE/SST
  RMSE = sqrt(SSE/nrow(df))
  
  data.frame(
    RMSE = RMSE,
    Rsquare = R_square
  )
}

##### 1. dataset setup #####
## 1.1. get all datasets including (WT+all mutatns) x (-8 to +4 TSSs)
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

## 1.2. use datasets of (-8 to +1, +4 TSSs)x(AYR+BYR+ARY)+ (+2 TSS)x(ARY)
master_mtx <- subset(master_mtx, TSS != '3') # 4322382      26
master_mtx <- master_mtx[!(grepl("[ACTG]YR",master_mtx$o_mtf811, perl=TRUE) & master_mtx$TSS == '2'), ] # 4010354

pol = 'WT'
sub_mtx <- subset(master_mtx, PolII == pol)

## 1.3. training 80%, test 20%
set.seed(42) # sets the random seed for reproducibility of results during pipline generation
perc <- 80/100

index = sample(1:nrow(sub_mtx), perc*nrow(sub_mtx))
train <- sub_mtx[index,]
test <- sub_mtx[-index,]

train$eff01 <- train$eff/100 # transform eff from [0,100] to [0,1]
test$eff01 <- test$eff/100 # transform eff from [0,100] to [0,1]

rm(perc,index)

##### 2. stepwise model selection, w/ cv, using caret #####
# below is for 2 stages, for additive terms and for interactions.
# when processing 2.2 additive, change if(FALSE) in 2.2 to if(TRUE) and keep if(FALSE) in 2.3
# when processing 2.3 interaction, change if(FALSE) in 2.3 to if(TRUE) and keep if(FALSE) in 2.2

## 2.2 additive parameters; for loop ----
if (FALSE){
set.seed(42) # can all use seeds= argument in trainControl()
additive <- c("n11","n10","n9","n8","n7","n6","n5","n4","n3","n2","n1","p1","p2","p3","p4","p5","p6","p7","p8","p9")
R2_best <- data.frame(row.names = additive) # to save R_square of “optimal” model, which is the best one within cv
R2_train <- data.frame(row.names = additive) # Rsquare of complete train set predicted by optimal model
R2_test <- data.frame(row.names = additive) # Rsquare of complete test set predicted by optimal model

cv <- 5
added_var <- c()
to_add_var <- c("n11","n10","n9","n8","n7","n6","n5","n4","n3","n2","n1","p1","p2","p3","p4","p5","p6","p7","p8","p9")

out_folder <- "/bgfs/ckaplan/Yunye/3-TSS_sequence_library/8-modeling-ct4/1-LR-add1_cv/"
while (length(to_add_var)>0){
  var_num <- as.character(length(added_var)+1)
  R2_cv <- data.frame(row.names = additive)
  
  for (add_var in to_add_var){
    if (length(added_var)){ # >0 added var, not the 1st round
      lg_temp <- train(as.formula(paste("eff01 ~ 1",
                                        paste(added_var, collapse="+"),
                                        add_var, sep = "+ ")),
                       data = train,
                       trControl = trainControl(method = "cv", number = cv), # define training control
                       method = "glm",
                       family="quasibinomial")
    } else { # =0 added var, 1st round, nothing in added_var
      lg_temp <- train(as.formula(paste("eff01 ~ 1",
                                        add_var, sep = "+ ")),
                       data = train,
                       trControl = trainControl(method = "cv", number = cv), # define training control
                       method = "glm",
                       family="quasibinomial")
    }
    
    R2_best[add_var, var_num] <- lg_temp$results$Rsquared
    R2_cv[add_var, 1:cv] <- lg_temp$resample$Rsquared
    R2_train[add_var, var_num] <- eval_results(train$eff01, predict(lg_temp, newdata = train), train)[1,"Rsquare"]
    R2_test[add_var, var_num] <- eval_results(test$eff01, predict(lg_temp, newdata = test), test)[1,"Rsquare"]
  }
  write.csv(R2_cv, paste0(out_folder,"LR-add1_",var_num,"-cv",as.character(cv),"-R2_idv.csv"))
  
  # optional: output temp results, kind of monitor process:
  write.csv(R2_best, paste0(out_folder,"LR-add1_",var_num,"-cv",as.character(cv),"-R2_best.csv"))
  write.csv(R2_train, paste0(out_folder,"LR-add1_",var_num,"-cv",as.character(cv),"-R2_train.csv"))
  write.csv(R2_test, paste0(out_folder,"LR-add1_",var_num,"-cv",as.character(cv),"-R2_test.csv"))
  
  best_para <- rownames(R2_best)[which.max(R2_best[,var_num])] # variable that has highest R2
  added_var <- c(added_var, best_para) # add best-variable into added_var
  to_add_var <- to_add_var[to_add_var != best_para] # remove best-variable from to_add_var
}
print(added_var) # order of added variable
write.csv(R2_best, paste0(out_folder,"LR-add1_20idv-cv",as.character(cv),"-R2_best.csv"))
write.csv(R2_train, paste0(out_folder,"LR-add1_20idv-cv",as.character(cv),"-R2_train.csv"))
write.csv(R2_test, paste0(out_folder,"LR-add1_20idv-cv",as.character(cv),"-R2_test.csv"))

} # 2.2 end

## 2.3 (9 additive +)interaction parameters; for loop ----
if (FALSE){
set.seed(42) # can all use seeds= argument in trainControl()
additive <- c("n9","n8","n7","n4","n3","n2","n1","p1","p2") # selected robust additive terms
# generate a list of all possible interactions
for (x in seq(length(additive)-1)){
  for (y in seq(from=x+1, to=length(additive))){
    if (exists("itrs")){
      itrs <- c(itrs, paste0(additive[x],"*",additive[y]))
    } else {
      itrs <- c(paste0(additive[x],"*",additive[y]))
    }
  }
}
R2_best <- data.frame(row.names = c(additive,itrs)) # to save R_square of “optimal” model, which is the best one within cv
R2_train <- data.frame(row.names = c(additive,itrs)) # Rsquare of complete train set predicted by optimal model
R2_test <- data.frame(row.names = c(additive,itrs)) # Rsquare of complete test set predicted by optimal model

cv <- 5
added_var <- additive[]
to_add_var <- itrs[] # candidates pool

out_folder <- "/bgfs/ckaplan/Yunye/3-TSS_sequence_library/8-modeling-ct4/1-LR-add1_cv/lg9i/"
while (length(to_add_var)>0){
  var_num <- as.character(length(added_var)+1)
  R2_cv <- data.frame(row.names = c(additive,itrs))
  
  for (add_var in to_add_var){
    lg_temp <- train(as.formula(paste("eff01 ~ 1", 
                                      paste(added_var, collapse="+"),
                                      add_var, sep = "+ ")),
                     data = train,
                     trControl = trainControl(method = "cv", number = cv), # define training control
                     method = "glm",
                     family="quasibinomial")
    R2_best[add_var, var_num] <- lg_temp$results$Rsquared
    R2_cv[add_var, 1:cv] <- lg_temp$resample$Rsquared
    R2_train[add_var, var_num] <- eval_results(train$eff01, predict(lg_temp, newdata = train), train)[1,"Rsquare"]
    R2_test[add_var, var_num] <- eval_results(test$eff01, predict(lg_temp, newdata = test), test)[1,"Rsquare"]
  }
  write.csv(R2_cv, paste0(out_folder,"LR-9i_add1_",var_num,"-cv",as.character(cv),"-R2_idv.csv"))
  
  # optional: output temp results, kind of monitor process:
  write.csv(R2_best, paste0(out_folder,"LR-9i_add1_",var_num,"-cv",as.character(cv),"-R2_best.csv"))
  write.csv(R2_train, paste0(out_folder,"LR-9i_add1_",var_num,"-cv",as.character(cv),"-R2_train.csv"))
  write.csv(R2_test, paste0(out_folder,"LR-9i_add1_",var_num,"-cv",as.character(cv),"-R2_test.csv"))
  
  best_para <- rownames(R2_best)[which.max(R2_best[,var_num])] # variable that has highest R2
  added_var <- c(added_var, best_para) # add best-variable into added_var
  to_add_var <- to_add_var[to_add_var != best_para] # remove best-variable from to_add_var
}
print(added_var) # order of added variable
write.csv(R2_best, paste0(out_folder,"LR-9i_add1_36itr-cv",as.character(cv),"-R2_best.csv"))
write.csv(R2_train, paste0(out_folder,"LR-9i_add1_36itr-cv",as.character(cv),"-R2_train.csv"))
write.csv(R2_test, paste0(out_folder,"LR-9i_add1_36itr-cv",as.character(cv),"-R2_test.csv"))

} # 2.3 end