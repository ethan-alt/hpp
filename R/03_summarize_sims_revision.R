remove(list = ls())

library(tidyverse)

setwd('~/Projects/Paper3/extensive_sims/Data')
histdata = readRDS('histdata.rds')
beta0    = histdata$betahat0['treatment']
remove('histdata')

## HPP loop
setwd('~/Projects/Paper3/extensive_sims/Results/sims_raw/hpp')
nfiles = length(list.files(pattern = '.rds'))
for ( i in 1:nfiles ) {
  filename  = paste0('res_', i, '.rds')
  if(!file.exists(filename))
    next
  temp      = readRDS(filename)
  hyp       = temp$hyp
  method    = temp$method
  n         = temp$n
  lambda    = temp$lambda
  id        = i
  temp      = temp$post
  beta.true = ifelse(hyp == 'h0', 0.0, beta0)
  ndatasets = sum(complete.cases(temp))
  temp      = temp[complete.cases(temp), ]
  bias      = mean( temp[,'mean'] - beta.true ) 
  mse       = mean( (temp[,'mean'] - beta.true)^2 )
  post.prob = mean(temp[, 'post.prob'] >= 0.975)
  hpd.width = mean( temp[,'p97.5'] - temp[, 'p2.5'] )
  hpd.cov   = mean( (beta.true > temp[, 'p2.5']) & (beta.true < temp[, 'p97.5']) )
  temp      = data.frame(
    id = id, ndatasets = ndatasets, method = method, n = n, lambda = lambda, 
    hyp = hyp, beta.true = beta.true, bias = bias, mse = mse, post.prob = post.prob,
    hpd.width = hpd.width, hpd.cov = hpd.cov
  )
  if ( i == 1 ) {
    res.hpp = temp
    next
  }
  res.hpp = rbind(res.hpp, temp)
}

res.hpp = res.hpp %>%
  group_by(method, n, lambda, hyp, beta.true) %>%
  summarize(
    bias      = weighted.mean(bias,      ndatasets),
    mse       = weighted.mean(mse,       ndatasets),
    post.prob = weighted.mean(post.prob, ndatasets),
    hpd.width = weighted.mean(hpd.width, ndatasets),
    hpd.cov   = weighted.mean(hpd.cov,   ndatasets),
    ndatasets = sum(ndatasets)
  )


res.napp <- readRDS('~/Projects/Paper3/extensive_sims/napp/res_napp.rds')
res.comm <- readRDS('~/Projects/Paper3/extensive_sims/commensurate/res_commensurate.rds')
res.hier <- rbind(res.napp, res.comm)
res.hier <- res.hier %>%
  rename(hpd.cov = ci.cov)

res <- rbind(res.hpp, res.hier)


res$compatibility = ifelse(res$hyp == 'h0', 'Incompatible', 'Compatible')


## compute relative MSE
res = res %>% 
  group_by(n, hyp) %>%
  mutate(mse.hpp.5 = mse[method == 'hpp' & lambda == 0.5]) %>%
  mutate(mse.hpp.5 = mse / mse.hpp.5)

res$beta.true <- ifelse(res$hyp == 'h0', 0, 0.04780295)
saveRDS(res, '~/Projects/Paper3/extensive_sims/Results/extensive_sims_revision.rds')


