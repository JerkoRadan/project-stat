


# ------------ Question A ------------------

#install.packages("SemiPar")

library(SemiPar)
library(MASS)

DataQA = read.table("covidDATA.txt")
modulo11=function(x) {x- floor(x/11)*11}
studentnumber=773435            
mycol=modulo11(studentnumber)
mydataA=DataQA[,c((1+mycol*4):(4+mycol*4))]


hist(mydataA$NEW_IN.3)
hist(log(mydataA$NEW_IN.3+1))

mydataA$DATE.3 <- as.Date(mydataA$DATE.3, '%Y-%m-%d')
mydataA['Days'] <- seq(from=1, to=56, by=1)
mydataA['log_new_in'] <- log(mydataA$NEW_IN.3+1)

attach(mydataA)

par(mfrow=c(1,1)) 
plot(DATE.3,log_new_in)

#Lm is not a good fit
model0 <- spm(log_new_in~Days, spar.method="ML")
summary(model0)

# Penalized regression splines with truncated polynomial basis of degree 1
model1 <- spm(log_new_in~ f(Days, basis="trunc.poly", degree=1))
summary(model1)
model1$fit$coef
model1$info$pen$knots

# Penalized regression splines with truncated polynomial basis of degree 2
model2 <- spm(log_new_in~ f(Days, basis="trunc.poly", degree=2))
summary(model2)
model2$fit$coef
model2$info$pen$knots

# Penalized regression splines with radial basis of degree 3
model3 <- spm(log_new_in~ f(Days))
summary(model3)
model3$fit$coef
model3$info$pen$knots


####### Part A.1 #######

par(mfrow=c(1,1))
plot(model1, se=F, ylim=c(min(log_new_in), max(log_new_in)), main= "Models estimators", ylab = "Number of new intakes",  xlab = "Days")
lines(model2, se=F, col="blue")
lines(model3, se=F, col="green")
legend("topright", legend=c("model 1", "model 2", "model 3"), col=c("black", "blue", "green"), lty = 1)
points(Days,log_new_in)

####### Part A.2 #######

predict1 <- predict(model1, newdata = mydataA, se=TRUE)
predict2 <- predict(model2, newdata = mydataA, se=TRUE)
predict3 <- predict(model3, newdata = mydataA, se=TRUE)

#Pointwise confidence interval using a t-distribution since the population variance is unknown
n=length(Days)
par(mfrow=c(1,1))
plot(model1, se=F, ylim=c(min(log_new_in), max(log_new_in)), main= "Model 1 estimator and confidence interval", ylab = "Number of new intakes",  xlab = "Days")
lines(unlist(Days),predict1$fit+qt(0.99,df=n-1)*predict1$se,col="blue", lty=2)
lines(unlist(Days),predict1$fit-qt(0.99,df=n-1)*predict1$se,col="blue", lty=2)
points(Days,log_new_in)

par(mfrow=c(1,1))
plot(model2, se=F, ylim=c(min(log_new_in), max(log_new_in)), main= "Model 2 estimator and confidence interval", ylab = "Number of new intakes",  xlab = "Days")
lines(unlist(Days),predict2$fit+qt(0.99,df=n-1)*predict2$se,col="blue", lty=2)
lines(unlist(Days),predict2$fit-qt(0.99,df=n-1)*predict2$se,col="blue", lty=2)
points(Days,log_new_in)

par(mfrow=c(1,1))
plot(model3, se=F, ylim=c(min(log_new_in), max(log_new_in)), main= "Model 3 estimator and confidence interval", ylab = "Number of new intakes",  xlab = "Days")
lines(unlist(Days),predict3$fit+qt(0.99,df=n-1)*predict3$se,col="blue", lty=2)
lines(unlist(Days),predict3$fit-qt(0.99,df=n-1)*predict3$se,col="blue", lty=2)
points(Days,log_new_in)

detach(mydataA)



# ------------ Question B ------------------
#install.packages("glmnet")

library(glmnet)


DataQB = read.table("BankDefaultData.txt",header=T)
studentnumber = 773435
set.seed(studentnumber)
rownumbers = sample(1:6436,size=1000)
mydataB = DataQB[rownumbers,]

attach(mydataB)

####### Part B.1 #######
mydataB$Term <- as.factor(Term)

plot(Income, Default)
plot(Employment, Default)
plot(Loan, Default)

# Full model 1
glm.full=glm(Default~.,data=mydataB, family=binomial)
summary(glm.full)

# Full model 2
glm.full2=glm(Default~.^2,data=mydataB, family=binomial)
summary(glm.full2)

# StepAIC Model 1
bm_aic1<- stepAIC(glm.full, k=2, direction='both', scope=list(upper=~.^2, lower=~1))
summary(bm_aic1)

# StepAIC Model 2
bm_aic2<- stepAIC(glm.full, k=2, direction='backward', scope=list(upper=~., lower=~1))
summary(bm_aic2)

# StepAIC Model 3
bm_aic3<- stepAIC(glm.full2, k=2, direction='backward', scope=list(upper=~.^2, lower=~1))
summary(bm_aic3)

all_models <- list(glm.full, glm.full2, bm_aic1, bm_aic2, bm_aic3)

tab = matrix(0, ncol=1, nrow=length(all_models))
for( i in 1:length(all_models)){
  tab[i,] = AIC(all_models[[i]], k = 2)
}

tab = cbind(round(tab, digits=2), c("glm.full","glm.full2", "Model 1", "Model 2", "Model 3"))
colnames(tab) = c("AIC", "model")
order_models <- tab[order(tab[,1]),]
order_models


####### Part B.2 #######
source("PostAICupdate.R")

# All subsets of glm.full including the intercept

# Postaic intervals
Postaic0.99 <- PostAIC(y=Default, mydataB[,-which(names(mydataB) == "Default")], model.set = "partsubsets", quant=0.99, family = binomial, common=c(), intercept = T, linearcomb = F )
Postaic0.99$`PostAIC intervals`

# Naive intervals (bm_aic2 is the same model as the one obtained by PostAIC, if it wasn't the case should fit the glm model with the corresponding features and levels)
summary(bm_aic2)
summary <- cbind(bm_aic2$coefficients, confint(bm_aic2),confint(bm_aic2)[,2]-confint(bm_aic2)[,1], Postaic0.99$'PostAIC intervals', Postaic0.99$'PostAIC intervals'[,2]-Postaic0.99$'PostAIC intervals'[,1])

colnames(summary) <- c("Estimate", "Naive 0.5%", "Naive 99.5%", "Naive CI width", "PostAIC 0.5%", "PostAIC 99.5%", "PostAIC CI width")
summary

detach(mydataB)

# ------------ Question C -----------------
library(MASS)
library(glmnet)

DataQC = read.table("bikestations.txt")
digitsum = function(x) sum(floor(x/10^(0:(nchar(x)-1)))%%10)
studentnumber=773435 
mysum = digitsum(studentnumber)
Y = DataQC[,mysum]
X = DataQC[,-mysum]
hist(Y)

####### Part C.1 #######

x<- as.matrix(X)
fit_glmnet <- glmnet(x = x, y = Y, family= "poisson")

#### (a) Ridge regression

alpha1 <- 0

cv_fit1 <- cv.glmnet(x, Y, alpha =alpha1, family= "poisson", nfold=7)
plot(cv_fit1)

cv_fit1$lambda.min #Lambda = 12.21
coef(fit_glmnet, s=cv_fit1$lambda.min)


#### (b) Lasso estimation

alpha2 <- 1

cv_fit2 <- cv.glmnet(x, Y, alpha =alpha2, family= "poisson", nfold=7)
plot(cv_fit2)

cv_fit2$lambda.min # Lambda = 0.4187
coef(fit_glmnet, s=cv_fit2$lambda.min)


#### (c) Elastic net with alpha = 0.5

alpha3 <- 0.5

cv_fit3 <- cv.glmnet(x, Y, alpha =alpha3, family= "poisson", nfold=7)
plot(cv_fit3)

cv_fit3$lambda.min # Lambda = 0.7284
coef(fit_glmnet, s=cv_fit3$lambda.min)

#### (d) Elastic net with alpha = 0.8

alpha4 <- 0.8

cv_fit4 <- cv.glmnet(x, Y, alpha =alpha4, family= "poisson", nfold=7)
plot(cv_fit4)

cv_fit4$lambda.min # Lambda = 0.5484
coef(fit_glmnet, s=cv_fit4$lambda.min)

coef_summary <- cbind(coef(fit_glmnet, s=cv_fit1$lambda.min), coef(fit_glmnet, s=cv_fit2$lambda.min), coef(fit_glmnet, s=cv_fit3$lambda.min), coef(fit_glmnet, s=cv_fit4$lambda.min))
colnames(coef_summary) <- c("Ridge (alpha=0)", "Lasso (alpha=1)", "Elastic net (alpha=0.5)", "Elastic net (alpha=0.8)")
coef_summary

glmnet:::cv.glmnet.raw
