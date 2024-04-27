library(data.table)

UKBB_Xtest <- fread("Data_Revision/UKBB (plasma) X_test.csv")
UKBB_Xtrain  <- fread("Data_Revision/UKBB (plasma) X_train.csv")
UKBB_ytest  <- fread("Data_Revision/UKBB (plasma) y_test.csv")
UKBB_ytrain  <- fread("Data_Revision/UKBB (plasma) y_train.csv")


################################
################################


#LASSO in HMP2 dataset

library(readr)
library(tidyverse)
library(tidymodels)
library(vip)
library(finalfit)
library(workflowsets)
library(finetune)
library(patchwork)
library(glmnet)
library(future.apply)
set.seed(132)

mainDir <-"/Users/bravol/Desktop/Animesh_Serena"
NameRun <- 'LASSO_Revision_UKBB'
subDir <- paste0(sub('\\..*', '', NameRun), format(Sys.time(), '_%Y%m%d_%H%M'))
dir.create(file.path(mainDir, subDir))
IncludeFigHere <- file.path(mainDir, subDir)

############### PARALLELIZATION
library(doParallel)
cl <- makeCluster(detectCores())
registerDoParallel(cl)
plan(multisession) # prev multiprocess
############### PARALLELIZATION



#no missing values
set.seed(45)
#TrackSplit <- initial_split(LASSO, strata = Label)
train <- UKBB_Xtrain %>%
  select(-c(V1)) %>%
  add_column(Label = UKBB_ytrain$Label) #%>%
#mutate(Label = ifelse(Label == "IBD", 1, 0))

test <- UKBB_Xtest %>%
  select(-c(V1)) %>%
  add_column(Label = UKBB_ytest$Label) #%>%
#mutate(Label = ifelse(Label == "IBD", 1, 0))

Bootstraps <- 400

boots <- bootstraps(train, times = Bootstraps) #apparent = FALSE

#source("/Users/lxb732/Desktop/Bladder/PlotsScript.R")
source("/Users/lxb732/Desktop/Bladder/ClassificationScript.R")
#source("/Users/lxb732/Desktop/1KIP_files/Regression_Model.R")

Plots <- function(LASSOModel,boots, Title, Minimum, IncludeFigHere) {
  #added boots!
  #Thresh > 10

  LASSOModelPlot <- LASSOModel %>%
    group_by(name) %>%
    add_tally() %>%
    mutate(nn= n/dim(boots)[1]) %>%
    filter(name != "(Intercept)") %>%
    ungroup() %>%
    group_by(Boots) %>%
    add_tally(name = "Num_Variables") %>%
    ungroup()

  Threshold <- LASSOModelPlot %>%
    select(name, n) %>%
    unique() %>%
    mutate(Quantile = median(c(quantile(n)[4], quantile(n)[5]))) #prev mean

  print(Threshold$Quantile[1])

  LASSOModelPlot <- LASSOModelPlot %>%
    mutate(Quantile = Threshold$Quantile[1]) %>%
    mutate(Thresh = ifelse(n > Threshold$Quantile[1], 1, 0))


  #write.csv(LASSOModelPlot, paste0(IncludeFigHere, "/Dataset2_",Title, ".csv" ))

  FeatureSelectPlot <- ggplot(LASSOModelPlot %>%
                                ungroup() %>%
                                filter(n>Minimum) %>% #chhaaaaange
                                select(name, n,Thresh) %>%
                                unique(), aes(reorder(name,+n),n,fill= as.factor(Thresh)))+
    geom_col()+ #geom_bar(stat="identity")
    theme(axis.text.x = element_text(size=8),legend.position="none")+
    labs(x="Features",y="Frequency")+
    theme(panel.background = element_blank(),panel.grid.minor = element_line(colour="black"),axis.line = element_line(colour = "black"))+
    coord_flip() + scale_fill_manual(values=c("#999999", "#E69F00"))

  Betas <- LASSOModelPlot %>%
    filter(Thresh == 1 )

  BetasPlot <- ggplot(Betas, aes(name, coefficient )) +
    geom_hline(yintercept=0, linetype="dashed",
               color = "red", size=0.5) +
    geom_boxplot() +
    # geom_point() + #toomany!
    coord_flip() +
    theme_bw() +
    labs(x="Features",y="Coefficients")


  RMSE <- LASSOModelPlot %>%
    select(Boots, performance, lambda.1se,Num_Variables) %>%
    unique() %>%
    ggplot(., aes(lambda.1se, performance, colour = Num_Variables)) + geom_point() + theme_bw() + labs(x = "lambda.1se", y = "AUC")


  pdf(paste0(IncludeFigHere,"/",Title, ".pdf"), 13, 9)
  print(FeatureSelectPlot + BetasPlot/RMSE)
  dev.off()

  return(LASSOModelPlot)

}




RunAll_Class <- function(boots,Vars,Minimum,IncludeFigHere){

  #Vars <- "Label"
  #Minimum <- 0.3*Bootstraps

  ClassModel <- future_lapply(1:dim(boots)[1], function(x) {ClassLASSO(boots$splits[[x]],1, Vars,IncludeFigHere)}, future.seed = TRUE) %>%
    bind_rows(.id = "Boots")


  #ClassModel <- lapply(1:dim(boots)[1], function(x) {ClassLASSO(boots$splits[[x]],1, Vars,IncludeFigHere)}) %>%
  #  bind_rows(.id = "Boots")

  #Plots(ClassModel, paste0("Class_LASSO_", Vars), Threshold)

  Data <- Plots(ClassModel, boots, paste0("Class_LASSO_", Vars), Minimum, IncludeFigHere)

  return(Data)

}


Data <- RunAll_Class(boots, "Label", 0.6*Bootstraps, IncludeFigHere)

save(Data, train, test, file = paste0(IncludeFigHere, "/AllData.RData"))

##############
library(data.table)

#IncludeFigHere <- "/Users/lxb732/Desktop/Bladder/Trial1000norm_20220615_0116"
#file <- list.files(path = paste0(IncludeFigHere,"/"), pattern = "*.csv")
#tbl_fread <- fread(paste0(IncludeFigHere, "/",file))


Select <- Data %>%
  select(name, Thresh) %>%
  filter(Thresh == 1)

Vars <- c(unique(Select$name), "Label")

Vars <-  gsub("`", "", Vars)

#problem if more than 1 model appears with same performance
Vars2 <- Data %>%
  filter(performance == max(performance)) %>%
  select(name)

Vars2 <- c(unique(Vars2$name), "Label")

#########

#store and generate booostraps - maybe decision trees?

library(cutpointr)


FinalCheck <- function(Vars, train, test,alpha){

  rec <- recipe(Label ~., data=train) %>% #[,1:40]
    #step_normalize(all_predictors()) %>%
    #step_nzv(all_predictors()) %>%
    step_select(all_of(Vars))

  CalValuesPrep <- rec %>% prep()
  ModMat <- juice(CalValuesPrep) #%>%

  X <- model.matrix(Label ~ ., data=ModMat)[, -1] #whats

  cv.fit <- cv.glmnet(X, ModMat$Label,  family = "binomial", alpha = alpha,standardize = TRUE)

  XTest2 <- bake(CalValuesPrep, test ) #%>%

  XTest <- model.matrix(Label ~ ., data= XTest2)[, -1] #whats

  ypred <- predict(cv.fit, newx =  XTest , s = "lambda.1se", type = "response")

  MSE <- Metrics::auc((XTest2$Label), as.numeric(ypred))

  Results1se <- coeff2dt(cv.fit, s = "lambda.1se") %>%
    add_column(lambda.1se = cv.fit$lambda.1se) %>%
    add_column(performance = MSE)

  #Confusion table

  CT <- data.frame(truth = (XTest2$Label), pred = as.numeric(ypred))

  ############


  return(list(CT,Results1se))

}

Mean <- FinalCheck(Vars, train, test,1)
Metrics <- Results(Mean[[1]], "LASSOMean") #default 1000 bootstraps
cutpoint <- filter(Metrics, metric == "optimal_cutpoint" )$Mean[1]
#obtain optimal cutpoint!

save(Mean, cutpoint, file =  paste0(IncludeFigHere,"/ModelCoeffs.RData"))

Multiply <- Mean[[2]]

############## Confusion Tablelibrary(ggplot2)

library(yardstick)
CT <- Mean[[1]] %>%
  mutate(predicted_response =ifelse(pred > cutpoint, 1, 0))

outcomes <- table(CT$predicted_response, CT$truth)
confusion <- conf_mat(outcomes)

pdf(paste0(IncludeFigHere,"/HeatmapLASSO.pdf"), 7, 5)
autoplot(confusion, type = "heatmap")
dev.off()

