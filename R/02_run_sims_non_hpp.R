remove(list = ls())

library(cmdstanr)
library(posterior)
# rstan_options(auto_write = TRUE)
# options(mc.cores = 1)

## MCMC sampling arguments
chains = 1
warmup = 1000
thin   = 1
nsmpl  = 15000
# iter   = warmup + thin * nsmpl

## sim / cluster arguments
n         = c(50, 75, 100)
hyp       = c('h0', 'h1')
method    = c('gpp', 'pp', 'ci')
# method    = 'hpp'
lambda     = c(0.5, 0.75, 1.0)
ndatasets  = 5000
each.cl    = 500
output.dir = '/pine/scr/e/t/ethanalt/Projects/Paper3'


# ## FOR TESTING--COMMENT OUT
# ndatasets  = 5
# each.cl    = 5 
# nsmpl      = 15000

save.every = min(25, each.cl)

grid       = expand.grid('n' = n, 'hyp' = hyp, 'method' = method, 'lambda' = lambda, stringsAsFactors = F)
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
filename.id = c('hierarchical2', 'gpp', 'powerprior', 'conjprior')[method.id == method]
filename.id = paste0(filename.id, '.stan')
remove(grid)

filename.id = NA
if(method.id == 'hpp') filename.id = 'hierarchical2.stan'
if(method.id == 'gpp') filename.id = 'gpp.stan'
if(method.id == 'pp')  filename.id = 'powerprior.stan'
if(method.id == 'ci')  filename.id = 'conjprior.stan'


## get save directory
save.dir = file.path('~/Projects/Paper3/extensive_sims/Results/sims_raw', ifelse(method.id == 'hpp', 'hpp', 'non_hpp'))


## read historical data
histdata = readRDS('~/Projects/Paper3/extensive_sims/Data/histdata.rds')
betahat0    = histdata$betahat0
betahat0.se = histdata$betahat0.se
formula     = histdata$fmla
histdata    = histdata$histdata
y0          = histdata$cd4
X0          = model.matrix(formula, histdata)
n0          = nrow(histdata)
remove(histdata)

## compile proper stan file
stanmod = cmdstan_model( file.path('~/Projects/Paper3/extensive_sims/Stan', filename.id) )


sample_stan = function(formula, data, lambda, ...) {
  freq.fit = glm(formula, poisson(), data)
  X        = model.matrix(freq.fit)
  standat = list(
    'y'           = data[, all.vars(formula)[1]],
    'X'           = X,
    'n'           = nrow(X),
    'p'           = ncol(X),
    'lambda'      = lambda,
    'post_sample' = 1
  )
  remove(X)
  if(method.id == 'hpp') {
    standat = c( standat,
      list(
        'mu0'         = mu0,
        'lambda0'     = lambda0,
        'start'       = as.vector(coef(freq.fit)),
        'maxit'       = 10,
        'tol'         = 1e-6,
        'post_sample' = 1
      )
    )
  } else if (method.id == 'gpp') {
    standat = c(standat,
      list(
        'betahat' = betahat0,
        'cov_betahat' = diag(betahat0.se^2)
      )
    )
  } else if (method.id == 'pp') {
    standat = c(standat,
      list(
       'n0' = n0,
       'y0' = y0,
       'X0' = X0
      )
    )
  } else if (method.id == 'ci') {
    standat = c(standat,
                list('mu0' = mu0)
    )
  } else { stop("method not recognized") }
  
  chains = as.numeric( list(...)$chains )
  if(is.null(chains))
    chains = 1
  
  if(method.id == 'hpp') {
    init = lapply(1:chains, function(x) list('beta' = as.vector(coef(freq.fit)), 'm' = as.vector(mu0)))
  } else
    init = lapply(1:chains, function(x) list('beta' = coef(freq.fit)) )
  
  ## return sampling
  # rstan::sampling(stanmod, standat, init = init, ...)
  stanmod$sample(data = standat, init = init, ...)
}

res = matrix(nrow = each.cl, ncol = 5)
colnames(res) = c('mean', 'sd', 'post.prob', 'p2.5', 'p97.5')
rownames(res) = indx.id
for ( i in 1:each.cl ) {
  ## get current data
  data    = readRDS(file.path('~/Projects/Paper3/extensive_sims/Data', hyp.id, paste0('data_', indx.id[i], '.rds')))
  mu0     = as.vector( data$mu0[1:n.id] )
  lambda0 = as.vector( data$lambda0[1:n.id] )
  data    = data$data[1:n.id, ]
  
  ## get samples from posterior (wrap in tryCatch to continue loop if any errors)
  ## get samples from posterior
  smpl = tryCatch(
      sample_stan(
        formula, data, lambda = lambda.id, chains = chains, iter_warmup = warmup, iter_sampling = nsmpl, thin = thin,
        refresh = 0, output_dir = output.dir, output_basename = paste0("smpl_", id, '_', indx.id[i] )
    )
    , error = function(e) NA
  )
  ## Compute posterior statistics
  if ('CmdStanMCMC' %in% class(smpl)) {
    smpl = tryCatch(
      smpl$summary('beta[2]', c(mean = mean, sd = sd, pp = ~mean(.x > 0), q=~quantile(.x, probs = c(0.025, 0.975))))[, -1],
      error = function(e) rep(NA, ncol(res))
    )
  } else {
    smpl = rep(NA, times = ncol(res))
  }
  res[i, ] = as.numeric(smpl)
  if ( i %% save.every == 0 ) {
    saveRDS(
      list('post' = res, 'id' = i, 'data.indx' = indx.id, 'lambda' = lambda.id, 'n' = n.id, 'method' = method.id, 'hyp' = hyp.id),
      file.path(save.dir, paste0('res_', id, '.rds'))
    )
  }
}



