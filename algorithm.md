# Grid Search Algorithm for Tree Methods in Generalized Propensity Score Model

## Major difference

- Incorporate big data with the H2O machine (multiple cores)
- Use the average of absolute correlations (AAC) as the stopping rule instead of mean squared error
- Grid search for (1) **Number of Trees** (`ntrees`), (2) **Maximum Tree Depth**(`max_depth`) and (3) **Column Sample Rate**(`col_sample_rate`)

## Basic idea

search for the best `ntrees` for each combination of `max_depth` and `col_sample_rate`

*Notes:*

- set the potential values of `max_depth` large due to the relatively large number of variables.
- set `ntrees` to a large enough value (e.g. `ntrees=25000`) so that we should not find the optimal at the boundary
- use the `h2o.staged_predict_proba` function to predict the GBM results for different numbers of trees that less than the number we set to `ntrees`
  - for each prediction result (e.g. when `ntrees=25000`, there will be 25000 prediction results each with different number of trees), do balancing checking and find the optimal `ntrees` for each combination of `max_depth` and `col_sample_rate`

## Steps

### Set default parameters

```R
min.rows <- 10 # fewest allowed (weighted) observations in a leaf
learn.rate <- 0.005 # learning rate 
rep.bootstrap <- 50 # number of bootstrap replicates
```

### Set the range of hyperparameters

```R
n.trees <- 25000
max.depth <- c(6, 8, 10, 12)
col.sample.rate <- c(0.8, 0.9, 1.0)
```

### Fit generalized propensity score (GPS) model, do prediction

```R
## take bc_30d exposure as an example
independent <- c("year","sex","married","mage","m_edu", "cigdpp","cigddp",
                   "clinega","kotck","pncgov", "rf_db_gest","rf_db_other",
                   "rf_hbp_chronic", "rf_hbp_pregn","rf_cervix","rf_prev_4kg",
                   "rf_prev_sga","firstborn","m_wg_cat",
                   "log_mhincome", "log_mhvalue", "percentPoverty",
                   "mrace_1", "mrace_2", "mrace_3", "mrace_4",
                   "bc_3090d", "bc_90280d")
birth.hex <- as.h2o(birth, destination_frame = "birth.hex")
gbm_30d <- h2o.gbm(y = "T",
                   x = independent,
                   training_frame = birth.hex,
                   ntrees = n.trees, 
                   max_depth = max.depth[1], # change
                   min_rows = min.rows,
                   learn_rate = learn.rate, 
                   col_sample_rate = rate, # change
                   distribution = "gaussian")
pred.gbm_30d <- h2o.staged_predict_proba(object = gbm_30d, newdata = birth.hex)
```

### Check balance and find the optimal `n.trees`

```R
## the function to calculate the AAC for different number of trees -
## - for each combination of max_depth and col_sample_rate
F.aac.iter <- function(i, data, data.hex, ps.num, ps.model.pred, rep) {
  # i: number of iterations (number of trees) 
  # data: dataset containing the treatment and the covariates not in h2o structure.
  # data.hex: dataset containing the treatment and the covariates in h2o env.
  # ps.model.pred: the staged prediction results of boosting model to estimate (p(T_iX_i)) 
  # ps.num: the estimated p(T_i) 
  # rep: number of replications in bootstrap 
  GBM.fitted <- as.vector(ps.model.pred[,floor(i)])
  ps.den <- dnorm((data$T - GBM.fitted)/sd(data$T - GBM.fitted),0,1)
  wtnt <- ps.num/ps.den
  wt <- wtnt
  wt <- fifelse(wt>quantile(wtnt, 0.99), quantile(wtnt, 0.99), wt)
  wt <- fifelse(wt<quantile(wtnt, 0.01), quantile(wtnt, 0.01), wt)
  aac_iter <- rep(NA,rep) 
  for (i in 1:rep){
    # set.seed(i)
    bo <- sample(1:dim(data)[1], size = floor((dim(data)[1])^0.7), 
                 replace = TRUE, prob = wt) # change
    newsample <- data[bo,]
    newsample <- Filter(function(x)(length(funique(x))>1), newsample)
    # j.drop <- match(c("T"), names(data))
    # j.drop <- j.drop[!is.na(j.drop)]
    x <- newsample %>% select(-T)
    ac <- matrix(NA,dim(x)[2],1)
    x <- as.data.frame(x)
    for (j in 1:dim(x)[2]){
      ac[j] = ifelse (!is.factor(x[,j]), stats::cor(newsample$T, x[,j],
                                                    method = "pearson"),
                      polyserial(newsample$T, x[,j]))
    }
    aac_iter[i] <- mean(abs(ac), na.rm = TRUE)
  }
  aac <- mean(aac_iter)
  return(aac)
}
```

- **calculate the stabilized inverse-probability weights (SIPW)** using the staged GPS model prediction results (saved in `pred.gbm_30d` 
  $$
  SIPW=\frac{f_X(X;\mu_1,\sigma^2_1)}{f_{X|C}(X|C=c;\mu_2,\sigma^2_2)}
  $$

- where $f_{\bullet}(\bullet)$ denotes the probability density function with mean $\mu$ and variance $\sigma^2$, $X$ the continuous exposure (black carbon exposure), and $C$ the set of confounders

- **truncate the SIPW** at 1\% and 99\%: those larger than 99 percentile are set to 99 percentile, less than 1 percentile are set to 1 percentile

- **calculate the balances, get AAC for each combination**
  
  - *do bootstrapping* to create data samples from the pseudo-population (bootstrap sampling with replacements with SIPW) for a designated number of replicates (`rep.bootstrap`)  
    - the sample size is set to `population^{0.7}` according to [references]
    - `rep.bootstrap = 50` in consideration of computation sources.
    - in each bootstrap sample, calculate the Pearson correlation coefficients for binary and continuous confounders, and calculate the polyserial correlation for the ordinal confounders
      - I am aware of that the most popular way to account for the balance results with the binary variable is to calculate the standardized mean difference, but according to <https://www.sciencedirect.com/topics/psychology/standardized-mean-difference>, the variance-accounted-for correlation is not really different statistically from it.
      - polyserial correlation is used by Zhu in his paper in 2015 *A Boosting Algorithm for Estimating Generalized Propensity Scores with Continuous Treatments* <https://pubmed.ncbi.nlm.nih.gov/26877909/>, and see <http://documentation.sas.com/?cdcId=pgmsascdc&cdcVersion=9.4_3.5&docsetId=procstat&docsetTarget=procstat_corr_details16.htm&locale=en> for method details
  - *average the absolute correlations over the bootstrapping replicates* as AAC
      - previously, I did the Fisher transformation for each correlation and then average it, and this put more weights on the extreme value. In this case, the extreme values are the correlation between exposures for different time windows.
      - while all the other balancing results are good according to previous results, we do not want to put all of efforts on the exposure correlation, which may result in overfitting/overadjustment for the other reasonable confounders
  
- **search the minimum of AACs**, and find its corresponding combination of `ntrees`, `max_depth` and `col_sample_rate`.

