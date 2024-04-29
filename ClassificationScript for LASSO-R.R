coeff2dt <- function(fitobject, s) {

  #Coefficient values

  coeffs <- coef(fitobject, s)
  coeffs.dt <- data.frame(name = coeffs@Dimnames[[1]][coeffs@i + 1], coefficient = coeffs@x)
  return(coeffs.dt[order(coeffs.dt$coefficient, decreasing = T),])
}


ClassLASSO <- function(bootsSingle, alpha, Vars,IncludeFigHere){

  #bootsSingle <- boots$splits[[1]]
  #Vars <- "Label"
  #alpha <- 1
  #Everything LASSO
  #FactorLabel <- function(x){
  #  as.numeric(as.factor(x))-1
  #}


  rec <- recipe(formula(paste(Vars, "~ .")), data=analysis(bootsSingle)) %>% #[,1:40]
    #step_normalize(all_predictors()) %>%
    step_nzv(all_predictors())

  CalValuesPrep <- rec %>% prep()
  ModMat <- juice(CalValuesPrep) #%>%
  #mutate_at(vars(starts_with("chr")), as.numeric)

  #X <- model.matrix(formula(paste(Vars, "~ .")), ~ ., data=ModMat)[, -1] #whats
  X <- model.matrix(formula(paste(Vars, "~ .")), data=ModMat)[, -1] #whats


  print("here")

  #X <- model.matrix(Diseases2 ~ ., data=analysis(bootsSingle)[,1:20])[, -1] #whats up with the hot encoding?
  cv.fit <- cv.glmnet(X, analysis(bootsSingle)[[Vars]],  family = "binomial", alpha = alpha,standardize = TRUE)
  #plot(cv.fit)

  #XTest <- model.matrix(diseases2 ~ ., data=assessment(bootsSingle)[,1:20])[, -1]

  #CalValuesPrep2 <- arec %>% prep(assessment(bootsSingle)) #[,1:40]
  XTest2 <- bake(CalValuesPrep, assessment(bootsSingle) ) #%>%
  #mutate_at(vars(starts_with("chr")), as.numeric)
  #XTest2 <- juice(CalValuesPrep2)
  #XTest <- model.matrix(formula(paste(Vars, "~ .")), ~ ., data= XTest2)[, -1] #whats
  XTest <- model.matrix(formula(paste(Vars, "~ .")), data= XTest2)[, -1] #whats

  ypred <- predict(cv.fit, newx =  XTest , s = "lambda.1se", type = "response")
  table(ypred)
  #MSE <- Metrics::auc(FactorLabel((assessment(bootsSingle)[[Vars]])), as.numeric(ypred))
  MSE <- Metrics::auc(((assessment(bootsSingle)[[Vars]])), as.numeric(ypred))
  #R2 <- Metrics::R2(as.numeric(assessment(bootsSingle)$BMI), #as.numeric(ypred))
  #MAE<- Metrics::mae(as.numeric(assessment(bootsSingle)$BMI), #as.numeric(ypred))

  Results1se <- coeff2dt(cv.fit, s = "lambda.1se") %>%
    add_column(lambda.1se = cv.fit$lambda.1se) %>%
    add_column(performance = MSE)
  #lambda.min


  return(list(Results1se))


}



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


Results <- function(Mean, Title){

  cp <- cutpointr(Mean, pred, truth,
                  method = maximize_metric, metric = sum_sens_spec, boot_runs = 1000, na.rm = TRUE) #youden

  Metrics <- c("AUC", "acc", "sensitivity", "specificity", "optimal_cutpoint")

  make_df <- function(cp, metric){

    a <- boot_ci(cp, !!metric, in_bag = FALSE, alpha = 0.1) %>% #in bag or outbag!
      add_column(metric = eval(metric)) %>%
      add_column(Mean = cp[[!!metric]])
    print(a)
    return(a)
  }

  dfs <- lapply(Metrics, function(x) {make_df(cp,x )}) %>%
    bind_rows()

  write.csv(dfs, paste0(IncludeFigHere, "/Metrics_",Title, ".csv" ))

  pdf(paste0(IncludeFigHere,"/AUC_",Title, ".pdf"), 7, 5)
  print(plot_cutpointr(cp, xvar = fpr, yvar = tpr, aspect_ratio = 1, conf_lvl = 0.95) + theme_bw() )
  dev.off()

  return(dfs)

}
