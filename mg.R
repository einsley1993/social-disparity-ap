

# ============== Mediational g formula for survival outcome ============== #
### convert BMI, income to quartiles
library(dplyr)
dat2 <- dat1[,c("id",
                "cvd",
                "survt",
                "end.int",
                "age",
                "minority",
                "sex",
                "educ",
                "PM_1",
                "PM_2",
                "PM_3",
                "PM_4",
                "PM_5" )]
dat2$BMI_q <- ntile(dat1$BMI, 4)
dat2$urban_1_q <- ntile(dat1$urban_1, 5)
dat2$urban_2_q <- ntile(dat1$urban_2, 5)
dat2$urban_3_q <- ntile(dat1$urban_3, 5)
dat2$urban_4_q <- ntile(dat1$urban_4, 5)
dat2$urban_5_q <- ntile(dat1$urban_5, 5)

dat2$income_1_q <- ntile(dat1$income_1, 5)
dat2$income_2_q <- ntile(dat1$income_2, 5)
dat2$income_3_q <- ntile(dat1$income_3, 5)
dat2$income_4_q <- ntile(dat1$income_4, 5)
dat2$income_5_q <- ntile(dat1$income_5, 5)

library(nlme)
library(lme4)
require(MASS)
require(Hmisc)

N <- 5000

indiv_pt_fun <- function(data){
  urban_pt <- PM_pt <- income_pt <- list()
  
  urban <- data[, c("urban_1_q",
                     "urban_2_q",
                     "urban_3_q",
                     "urban_4_q",
                     "urban_5_q")]
  income <- data[, c("income_1_q",
                     "income_2_q",
                     "income_3_q",
                     "income_4_q",
                     "income_5_q")]
  PM <- data[, c("PM_1",
                 "PM_2",
                 "PM_3",
                 "PM_4",
                 "PM_5")]
  
  for(i in 1:N){
    urban_pt[[i]] <- rep(NA, nrow = data$end.int[i])
    PM_pt[[i]] <- rep(NA, nrow = data$end.int[i])
    income_pt[[i]] <- rep(NA, nrow = data$end.int[i])
    for(j in 1:T){
      urban_pt[[i]][j] <- urban[i,j]
      PM_pt[[i]][j] <- PM[i,j]
      income_pt[[i]][j] <- income[i,j]
    }
  }
  
  indiv <- list()
  for(i in 1:N){
    indiv[[i]] <- as.data.frame(matrix(NA, nrow = T, ncol = 14))
    indiv[[i]][, c(1:9)] <- data[i, c("id", "cvd", "survt", "end.int",
                                      "minority", "sex", "educ", "age", "BMI_q")] 
    indiv[[i]][, 10] <- 0
    indiv[[i]][, 10][dat2$end.int[i]] <- ifelse(data$cvd[i] == 1, 1, 0) # CVD outcome at each time point
    indiv[[i]][, 11] <- 0:(T - 1) # time point
    indiv[[i]][, c(12:14)] <- c(urban_pt[[i]], income_pt[[i]], PM_pt[[i]])
    colnames(indiv[[i]]) <- c("id", "cvd", "survt", "end.int",
                              "minority", "sex", "educ", "age", "BMI_q",
                              "event_at_t", "t", 
                              "urban_q", "income_q", "PM")
  }
  
  # person-time pt format data
  indiv_pt <- do.call(rbind.data.frame, indiv)
  return(list(indiv_pt = indiv_pt,
              indiv = indiv))
}

indiv_pt <- indiv_pt_fun(data = dat2)[[1]]
indiv_ls <- indiv_pt_fun(data = dat2)[[2]]


############################### STEP 1 ###############################
# 1(a). For t > 1, the joint density of the observed confounders, exposures, and mediators at t, 
# given past covariate history and survival to t‐1.

## Exposure models: no need to model exposure so far
# temp1 <- polr(factor(income_1_q) ~ age, data = dat2, Hess = TRUE)
# myGuess <- c(temp1$coefficients, 0, temp1$zeta)
# exp_mod_1 <- polr(factor(income_1_q) ~ 
#                     age  + 
#                     # + factor(BMI_q) 
#                     + factor(minority) + factor(sex) 
#                     + factor(educ)
#                     + factor(urban_1_q) # Error in optim(s0, fmin, gmin, method = "BFGS", ...) : initial value in 'vmmin' is not finite
#                   , data = dat2, Hess = TRUE)
# pred_exp_1 <- predict(exp_mod_1, dat2, type = "class")
# 
# exp_mod_2 <- polr(factor(income_2_q) ~ 
#                     age + 
#                     # factor(BMI_q) +
#                     factor(minority) +
#                     + factor(sex) + 
#                     # + factor(educ) +
#                     factor(income_1_q)  
#                     + factor(urban_1_q)
#                   , data = dat2, Hess = TRUE)
# pred_exp_2 <- predict(exp_mod_2, dat2, type = "class")
# 
# exp_mod_3 <- polr(factor(income_3_q) ~ 
#                     age + 
#                     factor(BMI_q) +
#                     factor(minority) 
#                     + factor(sex)
#                     + factor(income_2_q)
#                     + factor(urban_2_q)
#                   , data = dat2, Hess = TRUE)
# pred_exp_3 <- predict(exp_mod_3, dat2, type = "class")
# 
# exp_mod_4 <- polr(factor(income_4_q) ~ 
#                     age +
#                     # factor(BMI_q) +
#                     factor(minority)
#                     + factor(sex)
#                     + factor(educ)
#                     + factor(income_3_q)
#                     + factor(urban_3_q)
#                   , data = dat2, Hess = TRUE)
# pred_exp_4 <- predict(exp_mod_4, dat2, type = "class")
# 
# exp_mod_5 <- polr(factor(income_5_q) ~ 
#                     factor(BMI_q)
#                     + factor(minority) 
#                     + factor(sex)
#                     + factor(educ)
#                     + factor(income_4_q)
#                     + factor(urban_4_q)
#                   , data = dat2, Hess = TRUE)
# pred_exp_5 <- predict(exp_mod_5, dat2, type = "class")
# 
# # acc:
# length(which(pred_exp_1 == dat2$income_1_q))/N
# length(which(pred_exp_2 == dat2$income_2_q))/N
# length(which(pred_exp_3 == dat2$income_3_q))/N
# length(which(pred_exp_4 == dat2$income_4_q))/N
# length(which(pred_exp_5 == dat2$income_5_q))/N


## Mediator models:
med_mod_1 <- lm(PM_1 ~ 
                  age + factor(BMI_q) + 
                  factor(minority) + factor(sex) +
                  factor(educ) + factor(income_1_q) + factor(urban_1_q), data = dat2)
summary(med_mod_1)

med_mod_2 <- lm(PM_2 ~ 
                  age + factor(BMI_q) + 
                  factor(minority) + factor(sex) + 
                  factor(educ) + 
                  factor(income_2_q) + factor(urban_2_q) + PM_1, 
                data = dat2)
summary(med_mod_2)

med_mod_3 <- lm(PM_3 ~ 
                  age + factor(BMI_q) + 
                  factor(minority) + factor(sex) + 
                  factor(educ) + 
                  factor(income_3_q) + factor(urban_3_q) + PM_2, 
                data = dat2)
summary(med_mod_3)

med_mod_4 <- lm(PM_4 ~ 
                  age + factor(BMI_q) + 
                  factor(minority) + factor(sex) + 
                  factor(educ) + 
                  factor(income_4_q) + factor(urban_4_q) + PM_3, 
                data = dat2)
summary(med_mod_4)

med_mod_5 <- lm(PM_5 ~ 
                  age + 
                  factor(BMI_q) + 
                  factor(minority) + factor(sex) + 
                  factor(educ) + 
                  factor(income_5_q) + factor(urban_5_q) + PM_4, 
                data = dat2)
summary(med_mod_5)



## Time varying confounder models:
urban_mod_1 <- polr(factor(urban_1_q) ~ 
                     age 
                     + factor(BMI_q)
                     + factor(minority) 
                     + factor(sex)
                     # + factor(educ)
                     ,
                    data = dat2, Hess = TRUE) 
pred_urban_1 <- predict(urban_mod_1, dat2, type = "class")


urban_mod_2 <- polr(factor(urban_2_q) ~ 
                      age 
                      + factor(BMI_q)
                      + factor(minority) + factor(sex)
                      # + factor(educ)
                      # + factor(urban_1_q)
                      ,
                    data = dat2, Hess = TRUE) 
pred_urban_2 <- predict(urban_mod_2, dat2, type = "class")

urban_mod_3 <- polr(factor(urban_3_q) ~ 
                      age + 
                      factor(BMI_q)
                      + factor(minority) + factor(sex)
                      # + factor(educ)
                      + factor(urban_2_q)
                      ,
                    data = dat2, Hess = TRUE) 
pred_urban_3 <- predict(urban_mod_3, dat2, type = "class")

urban_mod_4 <- polr(factor(urban_4_q) ~ 
                      age + factor(BMI_q) + 
                      factor(minority) + factor(sex)
                      # + factor(educ)
                      + factor(urban_3_q)
                      ,
                    data = dat2, Hess = TRUE) 
pred_urban_4 <- predict(urban_mod_4, dat2, type = "class")

urban_mod_5 <- polr(factor(urban_5_q) ~ 
                      age 
                      + factor(BMI_q)
                      + factor(minority) + factor(sex)
                      # + factor(educ)
                      # + factor(urban_4_q)
                      ,
                    data = dat2, Hess = TRUE) 
pred_urban_5 <- predict(urban_mod_5, dat2, type = "class")


# acc:
length(which(pred_urban_1 == dat2$urban_1_q))/N
length(which(pred_urban_2 == dat2$urban_2_q))/N
length(which(pred_urban_3 == dat2$urban_3_q))/N
length(which(pred_urban_4 == dat2$urban_4_q))/N
length(which(pred_urban_5 == dat2$urban_5_q))/N




### 1(b). For all t, the probability of surviving by t given past covariate history and survival to t‐1
# Fit logistic model to predict survival probability for each time point
pooled.logistic <- glm(event_at_t ~ 
                         age + factor(sex) + factor(BMI_q) + factor(educ) + PM
                         + factor(income_q)
                         + factor(urban_q) +
                         + t
                         + factor(income_q)*t
                       ,  family = "binomial", 
                       data = indiv_pt) 
summary(pooled.logistic)






############################### STEP 2 ###############################
# Goal: Estimate the joint distribution of the terfactual mediators, M(1:T)*, 
# for survivors as identified by Equation 4.
# (2a.i) For t ≥ 1,generate time t confounders, exposure, and mediator based on the estimated model coefficients
# of (1a) and previously generated covariates under the time‐varying exposure intervention a(1:t‐1).
# (2a.ii) Assign time t exposure under the intervention a(1:T).
# ??? I use interventional exposure for the code block below. Not sure why observed exposure is needed.

# newdat_A0: a(1:T)
newdat_A0 = data.frame(dat2)[,c(1:8, 14)]
newdat_A0$urban_1_q = predict(urban_mod_1, newdat_A0, type = "class")
# newdat_A0$income_1_q = predict(exp_mod_1, newdat_A0, type = "class")
newdat_A0$income_1_q <- 1 # interventional level
newdat_A0$PM_1 = predict(med_mod_1, newdat_A0, type = "response")

newdat_A0$urban_2_q = predict(urban_mod_2, newdat_A0, type = "class")
# newdat_A0$income_2_q = predict(exp_mod_2, newdat_A0, type = "class")
newdat_A0$income_2_q <- 1 # interventional level
newdat_A0$PM_2 = predict(med_mod_2, newdat_A0, type = "response")

newdat_A0$urban_3_q = predict(urban_mod_3, newdat_A0, type = "class")
# newdat_A0$income_3_q = predict(exp_mod_3, newdat_A0, type = "class")
newdat_A0$income_3_q <- 1 # interventional level
newdat_A0$PM_3 = predict(med_mod_3, newdat_A0, type = "response")

newdat_A0$urban_4_q = predict(urban_mod_4, newdat_A0, type = "class")
# newdat_A0$income_4_q = predict(exp_mod_4, newdat_A0, type = "class")
newdat_A0$income_4_q <- 1 # interventional level
newdat_A0$PM_4 = predict(med_mod_4, newdat_A0, type = "response")

newdat_A0$urban_5_q = predict(urban_mod_5, newdat_A0, type = "class")
# newdat_A0$income_5_q = predict(exp_mod_5, newdat_A0, type = "class")
newdat_A0$income_5_q <- 1 # interventional level
newdat_A0$PM_5 = predict(med_mod_5, newdat_A0, type = "response")         


# ===== Convert to person-time format ==== #
indiv_A0_pt <- indiv_pt_fun(data = newdat_A0)[[1]]
indiv_A0_ls <- indiv_pt_fun(data = newdat_A0)[[2]]

# check % of event:
# not 'cvd': this is not the event at each time point
nrow(indiv_A0_pt[which(indiv_A0_pt$event_at_t == 1),])/N

# use the logistic model in 1(b) to predict survival probability for each time point
indiv_A0_pt$cond_surv_p <- 1-predict(pooled.logistic, 
                                     newdata = indiv_A0_pt, 
                                     type = "response")
for_surv_A0 <- indiv_A0_pt[, c("id", "t", "income_q", "PM", "cond_surv_p")]


cum.surv.A0.mat <- as.data.frame(matrix(NA, nrow = N, ncol = 5))
colnames(cum.surv.A0.mat) <- c("Cum_surv_1", "Cum_surv_2", "Cum_surv_3", "Cum_surv_4", "Cum_surv_5")
for(i in 1:N){
  if(nrow(indiv_A0_ls[[i]]) > 1){
    for(j in 2:nrow(indiv_A0_ls[[i]])){
      cum.surv.A0.mat[i,1] <- for_surv_A0[which(for_surv_A0$id == i),]$cond_surv_p[1]
      cum.surv.A0.mat[i,j] <- for_surv_A0[which(for_surv_A0$id == i),]$cond_surv_p[j] * cum.surv.A0.mat[i,j-1]
    }
  }
}





######################### Loop #########################
K <- 100
ITE <- IDE <- IIE <- list()
set.seed(3104)
for(x in 2:5){
  ITE[[x]] <- IDE[[x]] <- IIE[[x]] <- rep(NA, K)
  ITE[[1]] <- IDE[[1]] <- ITE[[1]] <- rep(0, K)
  
  # (2c). Repeat (2a) replacing intervention a(1:T) with a(1:T)*.
  # set A1 level: A1
  # newdat_A1: a(1:T)*
  newdat_A1 = data.frame(dat2)[,c(1:8, 14)]
  newdat_A1$urban_1_q = predict(urban_mod_1, newdat_A1, type = "class")
  # newdat_A1$income_1_q = predict(exp_mod_1, newdat_A1, type = "class")
  newdat_A1$income_1_q <- x # interventional level
  newdat_A1$PM_1 = predict(med_mod_1, newdat_A1, type = "response")
  
  newdat_A1$urban_2_q = predict(urban_mod_2, newdat_A1, type = "class")
  # newdat_A1$income_2_q = predict(exp_mod_2, newdat_A1, type = "class")
  newdat_A1$income_2_q <- x # interventional level
  newdat_A1$PM_2 = predict(med_mod_2, newdat_A1, type = "response")
  
  newdat_A1$urban_3_q = predict(urban_mod_3, newdat_A1, type = "class")
  # newdat_A1$income_3_q = predict(exp_mod_3, newdat_A1, type = "class")
  newdat_A1$income_3_q <- x # interventional level
  newdat_A1$PM_3 = predict(med_mod_3, newdat_A1, type = "response")
  
  newdat_A1$urban_4_q = predict(urban_mod_4, newdat_A1, type = "class")
  # newdat_A1$income_4_q = predict(exp_mod_4, newdat_A1, type = "class")
  newdat_A1$income_4_q <- x # interventional level
  newdat_A1$PM_4 = predict(med_mod_4, newdat_A1, type = "response")
  
  newdat_A1$urban_5_q = predict(urban_mod_5, newdat_A1, type = "class")
  # newdat_A1$income_5_q = predict(exp_mod_5, newdat_A1, type = "class")
  newdat_A1$income_5_q <- x # interventional level
  newdat_A1$PM_5 = predict(med_mod_5, newdat_A1, type = "response")          
  
  
  
  # ===== Convert to person-time format ==== #
  indiv_A1_pt <- indiv_pt_fun(data = newdat_A1)[[1]]
  indiv_A1_ls <- indiv_pt_fun(data = newdat_A1)[[2]]
  
  # check % of event:
  # not 'cvd': this is not the event at each time point
  nrow(indiv_A1_pt[which(indiv_A1_pt$event_at_t == 1),])/N
  
  # use the logistic model in 1(b) to predict survival probability for each time point
  indiv_A1_pt$cond_surv_p <- 1-predict(pooled.logistic, 
                                       newdata = indiv_A1_pt, 
                                       type = "response")
  for_surv_A1 <- indiv_A1_pt[, c("id", "t", "income_q", "PM", "cond_surv_p")]
  
  
  cum.surv.A1.mat <- as.data.frame(matrix(NA, nrow = N, ncol = 5))
  colnames(cum.surv.A1.mat) <- c("Cum_surv_1", "Cum_surv_2", "Cum_surv_3", "Cum_surv_4", "Cum_surv_5")
  for(i in 1:N){
    if(nrow(indiv_A1_ls[[i]]) > 1){
      for(j in 2:nrow(indiv_A1_ls[[i]])){
        cum.surv.A1.mat[i,1] <- for_surv_A1[which(for_surv_A1$id == i),]$cond_surv_p[1]
        cum.surv.A1.mat[i,j] <- for_surv_A1[which(for_surv_A1$id == i),]$cond_surv_p[j] * cum.surv.A1.mat[i,j-1]
      }
    }
  }
  
  
  for(k in 1:K){
    # (2b). For each t = 1, ..., T, randomly permute the n values of the joint mediators 
    # assigned under intervention a(1:T) in (2a).
    rand_perm_med_A0_time1 = sample(newdat_A0$PM_1, length(newdat_A0$PM_1), replace = F)
    rand_perm_med_A0_time2 = sample(newdat_A0$PM_2, length(newdat_A0$PM_2), replace = F)
    rand_perm_med_A0_time3 = sample(newdat_A0$PM_3, length(newdat_A0$PM_3), replace = F)
    rand_perm_med_A0_time4 = sample(newdat_A0$PM_4, length(newdat_A0$PM_4), replace = F)
    rand_perm_med_A0_time5 = sample(newdat_A0$PM_5, length(newdat_A0$PM_5), replace = F)
    
    
    # 2(d). Repeat (2b) replacing intervention a(1:T) with a(1:T)*: don't need this for now
    # rand_perm_med_A1_time1 = sample(newdat_A1$PM_1, length(newdat_A1$PM_1), replace = F)
    # rand_perm_med_A1_time2 = sample(newdat_A1$PM_2, length(newdat_A1$PM_2), replace = F)
    # rand_perm_med_A1_time3 = sample(newdat_A1$PM_3, length(newdat_A1$PM_3), replace = F)
    # rand_perm_med_A1_time4 = sample(newdat_A1$PM_4, length(newdat_A1$PM_4), replace = F)
    # rand_perm_med_A1_time5 = sample(newdat_A1$PM_5, length(newdat_A1$PM_5), replace = F)
    
    ############################### STEP 3 ###############################
    # Goal: Recursively for each time t = 1,..., T and each subject i = 1, ..., n:
    # (3a.i) Repeat (2a.i)
    # we estimate Q(a1, G_a0) below, aka Q(a(1:T), a(1:T)*) 
    newdat_3ai = data.frame(dat2)[,c(1:8, 14)]
    newdat_3ai$urban_1_q = predict(urban_mod_1, newdat_3ai, type = "class")
    newdat_3ai$income_1_q <- x # interventional level
    newdat_3ai$PM_1 = rand_perm_med_A0_time1 # random draw
    
    newdat_3ai$urban_2_q = predict(urban_mod_2, newdat_3ai, type = "class")
    newdat_3ai$income_2_q <- x # interventional level
    newdat_3ai$PM_2 = rand_perm_med_A0_time2 # random draw
    
    newdat_3ai$urban_3_q = predict(urban_mod_3, newdat_3ai, type = "class")
    newdat_3ai$income_3_q <- x # interventional level
    newdat_3ai$PM_3 = rand_perm_med_A0_time3 # random draw
    
    newdat_3ai$urban_4_q = predict(urban_mod_4, newdat_3ai, type = "class")
    newdat_3ai$income_4_q <- x # interventional level
    newdat_3ai$PM_4 = rand_perm_med_A0_time4 # random draw
    
    newdat_3ai$urban_5_q = predict(urban_mod_5, newdat_3ai, type = "class")
    newdat_3ai$income_5_q <- x # interventional level
    newdat_3ai$PM_5 = rand_perm_med_A0_time5 # random draw
    
    
    # (3a.ii) Assign the time t mediator as the ith component of the permuted vector for time t from (2b)
    # (3a.iii) Assign time t exposure under the intervention a(1:T)1.
    # ??? the above steps are already done in (3a.i)
    
    
    # (3a.iv) Estimate the probability of surviving by t given past covariate history and survival to t‐1 given the
    # generated history in (3.a.i) through (3.a.iii) based on the estimated regression coefficients from (1b)
    
    # ===== Convert to person-time format ==== #
    indiv_3b_pt <- indiv_pt_fun(data = newdat_3ai)[[1]]
    indiv_3b_ls <- indiv_pt_fun(data = newdat_3ai)[[2]]
    
    # check % of event:
    # not 'cvd': this is not the event at each time point
    nrow(indiv_3b_pt[which(indiv_3b_pt$event_at_t == 1),])/N
    
    # use the logistic model in 1(b) to predict survival probability for each time point
    indiv_3b_pt$cond_surv_p <- 1 - predict(pooled.logistic, newdata = indiv_3b_pt, type = "response")
    for_surv_3b <- indiv_3b_pt[, c("id", "t", "income_q", "PM", "cond_surv_p")]
    
    
    cum.surv.3b.mat <- as.data.frame(matrix(NA, nrow = N, ncol = 5))
    colnames(cum.surv.3b.mat) <- c("Cum_surv_1", "Cum_surv_2", "Cum_surv_3", "Cum_surv_4", "Cum_surv_5")
    for(i in 1:N){
      if(nrow(indiv_3b_ls[[i]]) > 1){
        for(j in 2:nrow(indiv_3b_ls[[i]])){
          cum.surv.3b.mat[i,1] <- for_surv_3b[which(for_surv_3b$id == i),]$cond_surv_p[1]
          cum.surv.3b.mat[i,j] <- for_surv_3b[which(for_surv_3b$id == i),]$cond_surv_p[j] * cum.surv.3b.mat[i,j-1]
        }
      }
    }
    
    ITE[[x]][k] <- mean(cum.surv.A1.mat$Cum_surv_5) - mean(cum.surv.A0.mat$Cum_surv_5)
    IDE[[x]][k] <- mean(cum.surv.3b.mat$Cum_surv_5) - mean(cum.surv.A0.mat$Cum_surv_5)
    IIE[[x]][k] <- mean(cum.surv.A1.mat$Cum_surv_5) - mean(cum.surv.3b.mat$Cum_surv_5)
    
  }

  
}







range(ITE)
hist(IDE)
range(IDE)
hist(IIE)
range(IIE)




### Results
# higher income --> lower PM -- > lower risk
# ITE:
mean(cum.surv.A1.mat$Cum_surv_5) - mean(cum.surv.A0.mat$Cum_surv_5)
# -0.1800697
# IDE:
mean(cum.surv.3b.mat$Cum_surv_5) - mean(cum.surv.A0.mat$Cum_surv_5)
# -0.5843618
# IIE:
mean(cum.surv.A1.mat$Cum_surv_5) - mean(cum.surv.3b.mat$Cum_surv_5)
# 0.4042921

mean(IIE[[2]])
sd(IIE[[2]])
mean(IIE[[3]])
sd(IIE[[3]])
mean(IIE[[4]])
sd(IIE[[4]])
mean(IIE[[5]])
sd(IIE[[5]])

mean(ITE[[2]])
mean(ITE[[3]])
mean(ITE[[4]])
mean(ITE[[5]])









