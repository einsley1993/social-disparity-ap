

library(locfit)
seed <- 3104
set.seed(seed)
N <- 5000

id <- c(1:N)
# Treatment: time-fixed
minority <- rbinom(N, 1, 0.1)

# Time-fixed baseline covariates
age <- floor(rnorm(N, 50, 10))
hist(age)
sex <- rbinom(N, 1, 0.5)
educ <- rbinom(N, 1, 0.7) # < high school = 0, >= high school = 1
BMI <- rnorm(N, 23, 2)

set.seed(3104)
# Time-varying covariate: unit = 10k
urban <- list()
urban[[1]] <- 0.3*minority + 0.01*age - 0.2*educ + rnorm(1, 0, 1)
hist(urban[[1]])


# Time-varying exposure: unit = 10k
income <- list()
income[[1]] <- 2 + 0.3*minority + 0.3*age + 0.1*sex + 0.05*BMI + 0.6*educ + 0.3*urban[[1]] + rnorm(1, 0, 1)
hist(income[[1]])

# Mediator: time-varying. mean ~ 8
PM <- list()
PM[[1]] <- 8 + 0.3*minority - 0.05*age + 0.1*educ - 0.1*income[[1]] - 0.2*urban[[1]] + rnorm(1, 0, 0.1)
hist(PM[[1]])
mean(PM[[1]])





# Generate subsequent PM and income
# assume T time points
T <- 5
for(i in 2:T){
  urban[[i]] <- 0.3*minority + 0.01*age - 0.2*educ + 0.2*urban[[i-1]]
  income[[i]] <- -0.1 + 1*income[[i-1]] - 0.1*PM[[i-1]] - 0.2*urban[[i]] + rnorm(1, 0, 2)
  PM[[i]] <- 5 - 0.1*income[[i]] + 0.6*PM[[i-1]] + rnorm(1, 0, 1)
}

urban_df <- do.call(cbind.data.frame, urban)
income_df <- do.call(cbind.data.frame, income)
PM_df <- do.call(cbind.data.frame, PM)

for(i in 1:T){
  colnames(urban_df)[i] <- paste0("urban_", i)
  colnames(income_df)[i] <- paste0("income_", i)
  colnames(PM_df)[i] <- paste0("PM_", i)
}

apply(urban_df, 2, range)
apply(income_df, 2, range)
apply(PM_df, 2, range)









# Generate survival time
set.seed(3104)
survt <- exp(2 - 0.4*rowMeans(PM_df, na.rm = T) + rnorm(1, 2, 0.5))
# survt <- exp(0 + 0.3*sqrt(rowMeans(income_df, na.rm = T)) 
             # - 0.3*rowMeans(PM_df, na.rm = T) + rnorm(1, 2, 0.5))
hist(survt)
range(survt)

# Event
cvd <- ifelse(survt < T, 1, 0)
mean(cvd)
# 0.24

# censor obs beyond survt
urban_df_copy <- urban_df
PM_df_copy <- PM_df
income_df_copy <- income_df


for(i in 1:N){
  for(j in 1:T){
    if(j >= ceiling(survt[i])){
      urban_df_copy[i, j:T] <- NA
      PM_df_copy[i, j:T] <- NA
      income_df_copy[i, j:T] <- NA
    }
  }
}

id <- 1:N

dat1 <- cbind.data.frame(id, 
                        minority, # treatment
                        cvd,
                        sex, age, educ, BMI, # time-fixed covariates
                        urban_df_copy, # time-varying covariayes
                        income_df_copy, # time-varying exposure
                        PM_df_copy, # time-varying mediator
                        survt # survival time
)

# censor at time point T
dat1$end <- ifelse(dat1$survt <= T, survt, T)
dat1$end.int <- floor(dat1$end)




# ============== Distribution plots ============== #
library(ggplot2)
# income vs. time
library(reshape)
pdf("prelim/PM_income.pdf", width = 5, height = 4)
dat1_in <- dat1[, 13:17]
dat1_in$id <- 1:N
dat1_in_melt <- melt(dat1_in, id.vars = "id", 
                     measure.vars = colnames(dat1_in)[1:T])

ggplot(data = dat1_in_melt, aes(x = variable, y = value)) +
  geom_boxplot() +
  scale_x_discrete(breaks = colnames(dat1_in)[1:T],
                   labels = 1:T) +
  xlab("Time") +
  ylab("Income (Unit: 10k)") +
  ggtitle("Income distribution in 5 years")



dat1_PM <- dat1[, 18:22]
dat1_PM$id <- 1:N
dat1_PM_melt <- melt(dat1_PM, id.vars = "id", 
                     measure.vars = colnames(dat1_PM)[1:T])

ggplot(data = dat1_PM_melt, aes(x = variable, y = value)) +
  geom_boxplot() +
  scale_x_discrete(breaks = colnames(dat1_PM)[1:T],
                   labels = 1:T) +
  xlab("Time") +
  ylab("PM 2.5") +
  ggtitle("PM 2.5 distribution in 5 years")
dev.off()



dat <- dat1

# ============== Crude analyses ============== #
### 1. Kaplan-Meier curve
library(survival)
survival <- Surv(dat$end, dat$cvd) # 0: censored
survival
# + sign = censored individual
survival1 <- survfit(survival ~ 1, data = dat)
survival1

library(survminer)
surv_treat <- survfit(survival ~ as.factor(minority), data = dat) 
ggsurvplot(surv_treat, 
           data = dat,
           # conf.int = TRUE,
           xlim = c(1, 5),
           break.time.by = 1,
           xlab = "Time") +
  ggtitle("Crude Kaplan-Meier survival estimates") 
ggsave("prelim/KM.pdf", height = 5, width = 5)













