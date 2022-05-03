# remove(list = ls())
library(cmdstanr)
library(posterior)

cmdstanr::set_cmdstan_path('~/.cmdstanr/cmdstan-2.27.0')


## source comm prior R function
stanfile.path <- '~/Projects/Paper3/extensive_sims/commensurate'
setwd(stanfile.path)
source('glm_commensurate_R.R')


## MCMC sampling arguments
chains = 1
warmup = 1000
thin   = 1
nsmpl  = 20000

## sim / cluster arguments
n         = c(50, 75, 100)
hyp       = c('h0', 'h1')
method    = 'commensurate'
ndatasets  = 5000
each.cl    = 100
# each.cl    = 2   ## for testing
output.dir = '/pine/scr/e/t/ethanalt/Projects/Paper3'


# ## FOR TESTING--COMMENT OUT
# ndatasets  = 1
# each.cl    = 1
# nsmpl      = 15000

save.every = min(50, each.cl)

grid       = expand.grid('n' = n, 'hyp' = hyp, 'method' = method, 'lambda' = NA, stringsAsFactors = F)
grid       = grid[rep(seq_len(nrow(grid)), each = ndatasets / each.cl), ]
grid$end   = seq(from = each.cl, to = ndatasets, by = each.cl)
grid$start = grid$end - (each.cl - 1)

id = as.numeric( Sys.getenv('SLURM_ARRAY_TASK_ID') )
if(is.na(id)) {
  output.dir = NULL
  id = 1
}


n.id        = grid[id, 'n']
hyp.id      = grid[id, 'hyp']
method.id   = grid[id, 'method']
lambda.id   = grid[id, 'lambda']
start.id    = grid[id, 'start']
end.id      = grid[id, 'end']
indx.id     = start.id:end.id
remove(grid)


## get save directory
save.dir = file.path('~/Projects/Paper3/extensive_sims', method, 'results')

## read historical data
histdata = readRDS('~/Projects/Paper3/extensive_sims/Data/histdata.rds')
betahat0    = histdata$betahat0
betahat0.se = histdata$betahat0.se
formula     = histdata$fmla
histdata    = histdata$histdata
y0          = histdata$cd4
X0          = model.matrix(formula, histdata)
n0          = nrow(histdata)

tau         = rep(1, ncol(X0))

res = matrix(nrow = each.cl, ncol = 5)
colnames(res) = c('mean', 'sd', 'post.prob', 'p2.5', 'p97.5')
rownames(res) = indx.id
for ( i in 1:each.cl ) {
  ## get current data; subset to sample size
  data    = readRDS(file.path('~/Projects/Paper3/extensive_sims/Data', hyp.id, paste0('data_', indx.id[i], '.rds')))
  data    = data$data[1:n.id, ]
  
  ## Conduct sampling
  smpl = tryCatch(
    glm_commensurate(
      formula, poisson('log'), data = data, histdata = histdata, tau = tau, chains = chains, iter_warmup = warmup, iter_sampling = nsmpl, thin = thin,
      refresh = 0, output_dir = output.dir, output_basename = paste0("smpl_", id, '_', indx.id[i] )
    )
    , error = function(e) NA
  )
  
  
  if ("CmdStanMCMC" %in% class(smpl)) {
    smpl = smpl$summary('beta[2]', c(mean = mean, sd = sd, pp = ~mean(.x > 0), q=~quantile(.x, probs = c(0.025, 0.975))))
    smpl = as.numeric(smpl[,-1])
  } else {
    smpl = rep(NA, times = ncol(res))
  }
  res[i, ] = smpl
  if ( i %% save.every == 0 ) {
    saveRDS(
      list('post' = res, 'id' = i, 'data.indx' = indx.id, 'lambda' = lambda.id, 'tau' = tau, 'n' = n.id, 'method' = method.id, 'hyp' = hyp.id),
      file.path(save.dir, paste0('res_', id, '.rds'))
    )
  }
}

saveRDS(
  list('post' = res, 'id' = i, 'data.indx' = indx.id, 'lambda' = lambda.id, 'n' = n.id, 'method' = method.id, 'hyp' = hyp.id),
  file.path(save.dir, paste0('res_', id, '.rds'))
)

