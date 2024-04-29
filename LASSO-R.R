library(data.table)

UKBB_Xtest <- fread("Data_Revision/UKBB (plasma) X_test.csv")
UKBB_Xtrain  <- fread("Data_Revision/UKBB (plasma) X_train.csv")
UKBB_ytest  <- fread("Data_Revision/UKBB (plasma) y_test.csv")
UKBB_ytrain  <- fread("Data_Revision/UKBB (plasma) y_train.csv")

HMP2_Xtest <- fread("Data_Revision/HMP2 (feces) X_test.csv")
HMP2_Xtrain  <- fread("Data_Revision/HMP2 (feces) X_train.csv")
HMP2_ytest  <- fread("Data_Revision/HMP2 (feces) y_test.csv")
HMP2_ytrain  <- fread("Data_Revision/HMP2 (feces) y_train.csv")


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

mainDir <-"/Users/bravol/Desktop/Animesh_Serena" #change to yours
NameRun <- 'Together'
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

source("path/to/ClassificationScript.R") ## EDIT PATH HERE

Vars <- "Label"
Minimum <- 0.6*Bootstraps


ClassModel <- future_lapply(1:dim(boots)[1], function(x) {ClassLASSO(boots$splits[[x]],1, Vars,IncludeFigHere)}, future.seed = TRUE) %>%
  bind_rows(.id = "Boots")


#Thresh > 10

LASSOModelPlot <- ClassModel %>%
  bind_rows %>%
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
  theme(axis.text.x = element_text(size=12),legend.position="none")+
  theme(axis.text.y = element_text(size=12),legend.position="none")+
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
  theme(axis.text.x = element_text(size=12))+
  theme(axis.text.y = element_text(size=12))+
 theme(panel.background = element_blank()) +
  labs(x="Features",y="Coefficients")


RMSE <- LASSOModelPlot %>%
  select(Boots, performance, lambda.1se,Num_Variables) %>%
  unique() %>%
  ggplot(., aes(lambda.1se, performance, colour = Num_Variables)) + geom_point() + theme_bw() + labs(x = "lambda.1se", y = "AUC") +
  theme(axis.text.x = element_text(size=12))+
  theme(axis.text.y = element_text(size=12))+
  theme(panel.background = element_blank())

Title <- "LASSO"
LASSO_Fig <- FeatureSelectPlot + BetasPlot/RMSE

#pdf(paste0(IncludeFigHere,"/",Title, ".pdf"), 11, 9)
#print(FeatureSelectPlot + BetasPlot/RMSE)
#dev.off()
#
#png(paste0(IncludeFigHere,"/",Title, ".png"), 11, 9)
print(FeatureSelectPlot + BetasPlot/RMSE)
#dev.off()
#
ggsave(filename = paste0(IncludeFigHere,"/",Title, ".png"), width = 11, height = 6, device='png', dpi=700)

################################
################################
################################
################################




#no missing values
set.seed(45)
#TrackSplit <- initial_split(LASSO, strata = Label)
train <- HMP2_Xtrain %>%
  select(-c(V1)) %>%
  add_column(Label = HMP2_ytrain$Label) %>%
  select(-NUA)#%>%
#mutate(Label = ifelse(Label == "IBD", 1, 0))

test <- HMP2_Xtest %>%
  select(-c(V1)) %>%
  add_column(Label = HMP2_ytest$Label)  %>%
  select(-NUA)#%>%#%>%
#mutate(Label = ifelse(Label == "IBD", 1, 0))

Bootstraps <- 400

boots <- bootstraps(train, times = Bootstraps) #apparent = FALSE

#source("/Users/lxb732/Desktop/Bladder/PlotsScript.R")
source("~/Desktop/Animesh_Serena/Data_Revision/Figures/ClassificationScript.R") ## add
#source("/Users/lxb732/Desktop/1KIP_files/Regression_Model.R")


Vars <- "Label"
Minimum <- 0.6*Bootstraps


ClassModel <- future_lapply(1:dim(boots)[1], function(x) {ClassLASSO(boots$splits[[x]],1, Vars,IncludeFigHere)}, future.seed = TRUE) %>%
  bind_rows(.id = "Boots")


#Thresh > 10

LASSOModelPlot <- ClassModel %>%
  bind_rows %>%
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
  theme(axis.text.x = element_text(size=12),legend.position="none")+
  theme(axis.text.y = element_text(size=12),legend.position="none")+
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
  theme(axis.text.x = element_text(size=12))+
  theme(axis.text.y = element_text(size=12))+
  theme(panel.background = element_blank()) +
  labs(x="Features",y="Coefficients")


RMSE <- LASSOModelPlot %>%
  select(Boots, performance, lambda.1se,Num_Variables) %>%
  unique() %>%
  ggplot(., aes(lambda.1se, performance, colour = Num_Variables)) + geom_point() + theme_bw() + labs(x = "lambda.1se", y = "AUC") +
  theme(axis.text.x = element_text(size=12))+
  theme(axis.text.y = element_text(size=12))+
  theme(panel.background = element_blank())

Title <- "HMP2_NoNUA"

#pdf(paste0(IncludeFigHere,"/",Title, ".pdf"), 11, 9)
#print(FeatureSelectPlot + BetasPlot/RMSE)
#dev.off()
#
#png(paste0(IncludeFigHere,"/",Title, ".png"), 11, 9)
print(FeatureSelectPlot + BetasPlot/RMSE)
HMP2_NoNUA_Fig <- FeatureSelectPlot + BetasPlot/RMSE
#dev.off()
#
ggsave(filename = paste0(IncludeFigHere,"/",Title, ".png"), width = 11, height = 6, device='png', dpi=700)

################################
################################
################################
################################

#no missing values
set.seed(45)
#TrackSplit <- initial_split(LASSO, strata = Label)
train <- HMP2_Xtrain %>%
  select(-c(V1)) %>%
  add_column(Label = HMP2_ytrain$Label) #%>%
#mutate(Label = ifelse(Label == "IBD", 1, 0))

test <- HMP2_Xtest %>%
  select(-c(V1)) %>%
  add_column(Label = HMP2_ytest$Label) #%>%
#mutate(Label = ifelse(Label == "IBD", 1, 0))

Bootstraps <- 400

boots <- bootstraps(train, times = Bootstraps) #apparent = FALSE

Vars <- "Label"
Minimum <- 0.6*Bootstraps


ClassModel <- future_lapply(1:dim(boots)[1], function(x) {ClassLASSO(boots$splits[[x]],1, Vars,IncludeFigHere)}, future.seed = TRUE) %>%
  bind_rows(.id = "Boots")


#Thresh > 10

LASSOModelPlot <- ClassModel %>%
  bind_rows %>%
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
  theme(axis.text.x = element_text(size=12),legend.position="none")+
  theme(axis.text.y = element_text(size=12),legend.position="none")+
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
  theme(axis.text.x = element_text(size=12))+
  theme(axis.text.y = element_text(size=12))+
  theme(panel.background = element_blank()) +
  labs(x="Features",y="Coefficients")


RMSE <- LASSOModelPlot %>%
  select(Boots, performance, lambda.1se,Num_Variables) %>%
  unique() %>%
  ggplot(., aes(lambda.1se, performance, colour = Num_Variables)) + geom_point() + theme_bw() + labs(x = "lambda.1se", y = "AUC") +
  theme(axis.text.x = element_text(size=12))+
  theme(axis.text.y = element_text(size=12))+
  theme(panel.background = element_blank())

Title <- "HMP2"

#pdf(paste0(IncludeFigHere,"/",Title, ".pdf"), 11, 9)
#print(FeatureSelectPlot + BetasPlot/RMSE)
#dev.off()
#
#png(paste0(IncludeFigHere,"/",Title, ".png"), 11, 9)
print(FeatureSelectPlot + BetasPlot/RMSE)
HMP2Fig <- FeatureSelectPlot + BetasPlot/RMSE
#dev.off()
#
ggsave(filename = paste0(IncludeFigHere,"/",Title, ".png"), width = 11, height = 6, device='png', dpi=700)


#Fig1
#print(d / d & plot_annotation(tag_levels = 'A'))

library(ggpubr)
print(ggarrange(LASSO_Fig, HMP2Fig,
         labels = c("A", "B"),
         ncol = 1, nrow = 2,  common.legend = TRUE))

Title <- "Fig1Final"
ggsave(filename = paste0(IncludeFigHere,"/",Title, ".png"), width = 11, height = 12, device='png', dpi=700)





