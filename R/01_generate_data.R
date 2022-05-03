remove(list = ls())
library(MASS)

lambda = c(0.5, 0.75, 1.00)

n0        = 50
trt.prob  = 0.50
n         = c(50, 75, 100)
Nmax      = max(n)
Ndatasets = 5000

fmla     = cd4 ~ treatment + race + age
rhs.fmla = fmla[-2]

## read Data
setwd('~/Projects/Paper3/ACTG')
actg.hist = read.table('actg019.bios262', skip = 13)
actg.cur  = read.table('actg036.bios262', skip = 12)
colnames(actg.hist) = 
  colnames(actg.cur) = 
  c('age', 'treatment', 'race', 'cd4')

actg.cur$cd4 = round(actg.cur$cd4, 0)

fit.hist = glm.nb(fmla, actg.hist)
fit.cur  = glm.nb(fmla, actg.cur)
beta0    = coef(fit.cur)
beta     = beta0
beta['treatment'] = 0
trteff.true = c(beta0['treatment'], beta['treatment'])
m           = summary(fit.hist)$theta


## find data set that has significant effect for historical data
for ( i in 1:10000 ) {
  set.seed(i)
  ## generate covariates and treatment assignment
  age       = with(actg.hist, rlnorm(n0, meanlog = mean(log(age)), sdlog = sd(log(age))))
  race      = rbinom(n0, size = 1, prob = mean(actg.hist$race))
  treatment = rbinom(n0, size = 1, prob = trt.prob)
  hist      = data.frame(treatment, age, race)
  
  ## generate negative binomial response
  X0            = model.matrix(rhs.fmla, hist)
  # hist$cd4      = rpois(n0, lambda = exp(X0 %*% beta0))
  hist$cd4      = rnegbin(n0, mu = exp(X0 %*% beta0), theta = m)
  fit           = glm.nb(fmla, hist)
  pval          = summary(fit)$coefficients['treatment', 'Pr(>|z|)']
  betahat0       = coef(fit)['treatment']
  if ((betahat0 > 0) ) {
    norm   = sum( abs(betahat0 - trteff.true[1]) )
    if ( norm < 1e-3 ){
      print(i)
      break
    }
  }
}

sigmahat0 = summary(fit)$coefficients[, 'Std. Error']

## save historical data
res = list(
  'histdata'    = hist,
  'fmla'        = fmla,
  'rhs.fmla'    = rhs.fmla,
  'betahat0'    = coef(fit),
  'betahat0.se' = as.numeric(sigmahat0),
  'm'           = m,
  'mhat'        = fit$theta
)
setwd('~/Projects/Paper3/negbin/Data')
saveRDS(res, file = 'histdata.rds')



## construct historical data summary statistics
sigmahat0 = summary(fit)$coefficients[, 'Std. Error']
betahat0  = cbind(coef(fit), coef(fit))
betahat0['treatment', 2] = 0

dir = '~/Projects/Paper3/negbin/Data'
for ( j in 1:2) {
  beta = betahat0[, j]
  hyp  = c('h1', 'h0')[j]
  for ( i in 1:Ndatasets ) {
    repeat{
      ## Generate covariates
      data = data.frame(
        'treatment' = rbinom(Nmax, 1, 0.50),
        'race'      = rbinom(Nmax, 1, 0.90),
        'age'       = rlnorm(Nmax, meanlog = mean(log(actg.hist$age)), sdlog = sd(log(actg.hist$age)))
      )
      prop.race = mean(data$race[1:n[1]])
      if(prop.race <= 0.95)
        break
    }
    Xnew      = model.matrix(rhs.fmla, data)
    data$cd4  = rpois( n = nrow(Xnew), lambda = exp( Xnew %*% beta ) )
    etahat0   = Xnew %*% betahat0[, 1]   ## prior prediction always based on historical data
    mu0       = exp(etahat0)
    tau0      = mu0^2 * (Xnew^2 %*% sigmahat0^2)
    lambda0   = mu0 / tau0

    saveRDS(
      object = list(
        'data'     = data,
        'formula'  = fmla,
        'mu0'      = mu0,
        'lambda0'  = lambda0
      )
      , file = file.path(dir, hyp, paste0('data_', i, '.rds'))
    )
  }
}



