## This is the sample codes from Zhu's paper
## use Average absolute correlation coeffect (AAC) as the stopping rule
## only incorportate with small dataset, use gbm package

F.aac.iter <- function(i, data, ps.model, ps.num, rep, criterion) { 
  # i: number of iterations (trees) 
  # data: dataset containing the treatment and the covariates 
  # ps.model: the boosting model to estimate p(T_iX_i) 
  # ps.num: the estimated p(T_i) 
  # rep: number of replications in bootstrap 
  # criterion: the correlation metric used as the stopping criterion 
  GBM.fitted <- predict(ps.model, newdata = data, n.trees = floor(i),
                        type = "response") 
  ps.den <- dnorm((data$T - GBM.fitted)/sd(data$T - GBM.fitted),0,1)
  wt <- ps.num/ps.den
  aac_iter <- rep(NA,rep) 
  for (i in 1:rep){ 
    bo <- sample(1:dim(data)[1], replace = TRUE, prob = wt) # change
    newsample <- data[bo,] 
    j.drop <- match(c("T"),names(data)) 
    j.drop <- j.drop[!is.na(j.drop)] 
    x <- newsample[,-j.drop] 
    if(criterion == "spearman" | criterion == "kendall"){
      ac <- apply(x, MARGIN = 2, FUN = cor, y = newsample$T, 
                  method = criterion)
    } else if (criterion == "distance"){
        ac <- apply(x, MARGIN = 2, FUN = dcor, y = newsample$T)
        } else if (criterion == "pearson"){ 
          ac <- matrix(NA,dim(x)[2],1) 
          for (j in 1:dim(x)[2]){ 
            ac[j] = ifelse (!is.factor(x[,j]), cor(newsample$T, x[,j], 
                                                   method = criterion),
                            polyserial(newsample$T, x[,j])) 
          }
        } else print("The criterion is not correctly specified") 
    aac_iter[i] <- mean(abs(1/2*log((1+ac)/(1-ac))), na.rm = TRUE)
  }
  aac <- mean(aac_iter)
  return(aac)
}
  
# Find the optimal number of trees using Pearson/polyserial correlation 
library(gbm) 
library(polycor) 
mydata <- data.frame(T = M2WTCON, X = x) 
model.num = lm(T~1, data = mydata) 
ps.num <- dnorm((mydata$T-model.num$fitted)/(summary(model.num))$sigma,0,1)
model.den <- gbm(T~.,data = mydata, shrinkage = 0.0005, 
                 interaction.depth = 4, distribution = "gaussian",
                 n.trees = 20000) 
opt <- optimize(F.aac.iter,interval = c(1,20000), data = mydata, 
                ps.model = model.den, ps.num = ps.num, rep = 50, 
                criterion = "pearson") 
best.aac.iter <- opt$minimum
best.aac <- opt$objective
          
          
          
          