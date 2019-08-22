source("data_prep.R")
source("evaluate_cross_validation.R")

model.dump.string <- getwd()

#### Evaluate RMSPE ####

#model <- "lasso"

compare_prior_rmspe <- function(model) {
  
  prior.mods <- list.files(model.dump.string)
  prior.mods <- prior.mods[grepl(paste0(model, "_lambda_k_sims_"), prior.mods)]
  prior.mods <- prior.mods[grepl("basic", prior.mods) == FALSE]
  
  ## Loop through priors and get relevant information
  
  mod.rmspe <- list()

  for (i in 1:length(prior.mods)) {
    
    par.k.mods <- readRDS(paste0(model.dump.string, prior.mods[i])) 

    # Get prior settings
    
    lambda.r <- par.k.mods[[1]]$r
    lambda.d <- par.k.mods[[1]]$d
        
    # Get RMSPE
    
    pp.rmspe <- calculate_pp_k_rmspe(par.k.mods)
    
    map.rmspe <- calculate_map_k_rmspe(par.k.mods)
    
    d <- data.table(r = lambda.r, 
                    d = lambda.d,
                    pp_rmspe = pp.rmspe$rmse,
                    map_rmspe = map.rmspe$rmse)
    
    mod.rmspe[[i]] <- d
    
    }
    
   mod.rmspe <- rbindlist(mod.rmspe)  
  
   # Plot grid of RMSPE
   
   l.mod.rmspe <- data.table(gather(mod.rmspe, index, value, pp_rmspe:map_rmspe))
   l.mod.rmspe[, type := "Posterior predictive"]
   l.mod.rmspe[index == "map_rmspe", type := "MAP predictive"]
   l.mod.rmspe[, r := as.character(format(r, scientific = F))]
   l.mod.rmspe[, d := as.character(format(d, scientific = F))]
   
   p.map <- ggplot(l.mod.rmspe[index == "map_rmspe",], aes(x = r, y = d, fill = value)) + 
              geom_tile() + facet_wrap(~ type) +
                geom_text(aes(label = round(value, 4))) + 
                  theme_bw(base_size = 12) +
                    theme(strip.background = element_rect(fill = "white")) +
                      scale_fill_continuous(name = "RMSPE", high = "grey", low = "white") +
                        xlab("r") + ylab("delta")
                    
   p.pp <- ggplot(l.mod.rmspe[index == "pp_rmspe",], aes(x = r, y = d, fill = value)) + 
              geom_tile() + facet_wrap(~ type) +
                geom_text(aes(label = round(value, 4))) + 
                  theme_bw(base_size = 12) +
                    theme(strip.background = element_rect(fill = "white")) +
                      scale_fill_continuous(name = "RMSPE", high = "grey", low = "white") +
                         xlab("r") + ylab("delta")                
   
   # Plot histograms
   
   p.hist <- ggplot(l.mod.rmspe, aes(x = value)) + 
              geom_histogram(bins = 10, colour = "black", fill = "white") + 
                facet_wrap(~type, scales = "free") +
                    theme_bw(base_size = 12) +
                      theme(strip.background = element_rect(fill = "white")) +
                        xlab("RMSPE") + ylab("Count")
     
   r.list <- list(mod.rmspe = mod.rmspe,
                  p.pp = p.pp,
                  p.map = p.map,
                  p.hist)
   
   return(r.list)
   
  }

#### Evaluate MAP estimates ####

#model <- "lasso"
   
compare_prior_maps <- function(model) {
  
  prior.mods <- list.files(model.dump.string)
  prior.mods <- prior.mods[grepl(paste0(model, "_lambda_sims_"), prior.mods)]
  prior.mods <- prior.mods[grepl("basic", prior.mods) == FALSE]
  
  ## Loop through priors and get relevant information
  
  par.mod.breakdowns <- list()
  par.mod.densities <- list()
  
  for (i in 1:length(prior.mods)) {
    
    par.mods <- readRDS(paste0(model.dump.string, prior.mods[i])) 
    
    # Get prior settings
    
    lambda.r <- par.mods$r
    lambda.d <- par.mods$d

    # Load frequentist model (for convenience)
    
    freq.mod <- readRDS(paste0(model.dump.string, "freq_",model,"_model_basic.RDS"))
                            
    # Get parameter estimates
    
    par.par <- estimate_parameters(par.mods, lambda = TRUE, sigma_sq = TRUE, 
                                   freq.model = freq.mod)

    par.breakdown <- data.table(par.par$par.summary)
    
    par.densities <- data.table(par.par$par.density)
    
    par.breakdown[, r := lambda.r]
    par.breakdown[, d := lambda.d]
    
    par.densities[, r := lambda.r]
    par.densities[, d := lambda.d] 
    
    par.mod.breakdowns[[i]] <- par.breakdown
    par.mod.densities[[i]] <- par.densities
    
    }

    par.mod.breakdowns <- rbindlist(par.mod.breakdowns)
    
    par.mod.densities <- rbindlist(par.mod.densities)     
    
    p.points <- ggplot(par.mod.breakdowns, aes(x = var, y = map)) + 
                  geom_violin() +
                    geom_dotplot(binaxis = "y", stackdir = "center") + 
                      facet_wrap(~var, scales = "free")

    r.list <- list(pars = par.mod.breakdowns,
                   dens = par.mod.densities,
                   p.points = p.points)
    
    return(r.list)
    
    }

#### Ridge

ridge.prior.pp.rmspe <- compare_prior_rmspe("ridge")

ridge.prior.map.rmspe <- compare_prior_maps("ridge")
 
#### LASSO

lasso.prior.pp.rmspe <- compare_prior_rmspe("lasso")

lasso.prior.map.rmspe <- compare_prior_maps("lasso")

#### Save results

prior.tests <- list(ridge.prior.pp.rmspe = ridge.prior.pp.rmspe,
                          ridge.prior.map.rmspe = ridge.prior.map.rmspe,
                          lasso.prior.pp.rmspe = lasso.prior.pp.rmspe,
                          lasso.prior.map.rmspe = lasso.prior.map.rmspe)

saveRDS(prior.tests, "prior_test_results.RDS")



   
    