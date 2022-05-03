

library(dplyr)

histdata <- readRDS('~/Projects/Paper3/extensive_sims/Data/histdata.rds')
truth.h1 <- histdata$betahat0['treatment']

save.dir <- '~/Projects/Paper3/extensive_sims/commensurate'
res.dir  <- file.path(save.dir, 'results')
nfiles   <- length(list.files(res.dir, '.rds'))

for ( i in 1:nfiles ) {
  res_i <- readRDS( file.path( res.dir, paste0('res_', i, '.rds') ) )
  df_i <- data.frame(
    'method' = res_i$method, 'n' = res_i$n, 'hyp' = res_i$hyp, 
    'data.indx' = res_i$data.indx, res_i$post
  )
  if ( i == 1 ) {
    res = df_i
  } else {
    res = rbind(res, df_i)
  }
}

res.summ <- res %>%
  group_by(method, n, hyp) %>%
  summarise(
      'bias'   = mean( mean - truth )
    , 'mse'    = mean( (mean - truth)^2 )
    , 'pow'    = mean( (post.prob > 0.95) )
    , 'ci.cov' = mean( if_else( (truth >= p2.5) & (truth <= p97.5), 1, 0 ) )
  )

saveRDS(res.summ, file.path(save.dir, 'res_commensurate.rds'))


