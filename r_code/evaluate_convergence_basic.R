source("data_prep_basic.R")

model.dump.string <- getwd()

#### Define trace plot functions

build_traceplots <- function(sims, lambda = TRUE, sigma_sq = TRUE, tau_sq = TRUE) {
  
  # Extract traceplot information
  
  par.list <- paste0("beta[", seq(1, sims$p, 1), "]")
  
  if(sigma_sq == TRUE) {par.list <- c(par.list, "sigma_sq")}
  
  if(lambda == TRUE) {par.list <- c(par.list, "lambda")}
  
  if(tau_sq == TRUE) {par.list <- c(par.list, paste0("tau_sq[", seq(1, sims$p, 1), "]"))}
  
  t.plot <- traceplot(sims$model.object, pars = par.list)
  
  t.data <- data.table(t.plot$data)
  
  # Getting variable names
  
  beta.names <- data.table(parameter = paste0("beta[", seq(1, sims$p, 1), "]"),
                           var_name = colnames(sims$X))
  
  t.data[, parameter := as.character(parameter)]
  t.data <- data.table(left_join(t.data, beta.names, by = "parameter"))
  t.data[is.na(var_name), var_name := parameter] 
  
  var.order <- colnames(sims$X)
  var.order <- c(var.order, "lambda", "sigma_sq")
  
  if(tau_sq == TRUE) {var.order <- c(var.order, paste0("tau_sq[", seq(1, sims$p, 1), "]"))}
  
  t.data[, var_name := factor(var_name, levels = var.order)]
  
  # Thin data for nicer plotting
  
  t.thin <- t.data[seq(0, nrow(t.data), by = 15),]
  
  # Trace plot main parameters
  
  tp <- ggplot(t.thin[grepl("tau", var_name) == FALSE,], aes(x = iteration, y = value, colour = as.character(chain))) +
    geom_line(size = 0.3) +
    facet_wrap(~ var_name, scales = "free_y") + 
    theme_bw(base_size = 12) +
    theme(strip.background = element_rect(fill = "white")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(legend.position = "bottom") +
    scale_colour_grey(name = "Chain") +
    xlab("Iteration") +
    ylab("Value")    
  
  # Trace plot tau sq parameters
  
  if(tau_sq == TRUE) {
    
    tp.tau <- ggplot(t.thin[grepl("tau", var_name) == TRUE,], aes(x = iteration, y = value, colour = as.character(chain))) +
      geom_line(size = 0.3) +
      facet_wrap(~ var_name, scales = "free_y") + 
      theme_bw(base_size = 12) +
      theme(strip.background = element_rect(fill = "white")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      theme(legend.position = "bottom") +
      scale_colour_grey(name = "Chain") +
      xlab("Iteration") +
      ylab("Value")
    
  }
  
  # Extract Rhat
  
  r.sum <- data.table(summary(sims$model.object)$summary, keep.rownames = TRUE)
  r.sum <- r.sum[, c("rn", "Rhat")]
  
  setnames(r.sum, "rn", "parameter")
  
  r.sum <- data.table(left_join(r.sum, beta.names, by = "parameter"))
  r.sum[is.na(var_name), var_name := parameter] 
  
  r.sum <- r.sum[grepl("y_tilde", parameter) == FALSE,]
  r.sum <- r.sum[parameter != "lp__",]
  
  r.sum[, var_name := factor(var_name, levels = var.order)]          
  r.sum <- r.sum[order(var_name),]
  
  # Replace Rhat with split-Rhat - hack in
  
  split.r.list <- list()
  
  for (i in 1:length(par.list)) {
    
    t.par <- t.data[parameter == par.list[i],]
    t.par <- data.table(spread(t.par, chain, value))
    
    split.r <- Rhat(as.matrix(t.par[, 4:7]))
  
    split.r.list[[i]] <- data.table(parameter = par.list[i], Rhat = split.r)
    
    }
    
  split.r.list <- rbindlist(split.r.list)
  
  r.sum[, Rhat := NULL]
  
  r.sum <- data.table(inner_join(r.sum, split.r.list, by = "parameter"))
  
  # Plot violin of R hats
  
  vp <- ggplot(r.sum, aes(x = "parameter", y = Rhat)) + 
    geom_violin() +
    geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5) +
    geom_hline(aes(yintercept = 1), colour = "black", linetype = "dashed") + 
    theme_bw(base_size = 12) +
    xlab("") + ylab("R hat")
  
  r.list <- list(tp = tp,
                 vp = vp,
                 r.sum = r.sum)
  
  if(tau_sq == TRUE) {r.list$tp.tau <- tp.tau}
  
  return(r.list)
  
}

#### Load models 

blm.sims <- readRDS(paste0(model.dump.string, "blm_sims_basic.RDS"))

ridge.sims <- readRDS(paste0(model.dump.string, "ridge_sims_basic.RDS"))

lasso.sims <- readRDS(paste0(model.dump.string, "lasso_sims_basic.RDS"))


# summary(lasso.sims$model.object)$summary
# 
# Rhat(lasso.sims$model.object, ""
# 
#      traceplot(lasso.sims$model.object, pars = "beta[1]")
#   
#      t.plot <- traceplot(lasso.sims$model.object, pars = par.list)   
#      
# test <- rstan::extract(lasso.sims$model.object, "beta[1]")[[2]]

# test <- as.matrix(spread(traceplot(lasso.sims$model.object, pars = "lambda")$data, chain, value)[, 3:6])
# 
# 
# Rhat(test)
#### Look at trace-plots

### BLM

blm.conv <- build_traceplots(blm.sims, lambda = FALSE, sigma_sq = TRUE, tau_sq = FALSE)

# blm.conv$tp
# 
# blm.conv$r.sum

#max(blm.conv$r.sum$Rhat)

### Ridge

ridge.conv <- build_traceplots(ridge.sims, lambda = TRUE, sigma_sq = TRUE, tau_sq = FALSE)

# ridge.conv$tp
# 
# ridge.conv$r.sum

#max(ridge.conv$r.sum$Rhat)

### LASSO

lasso.conv <- build_traceplots(lasso.sims, lambda = TRUE, sigma_sq = TRUE, tau_sq = TRUE)
# 
# lasso.conv$tp
#  
# lasso.conv$r.sum
#
#lasso.conv$tp.tau

#max(lasso.conv$r.sum$Rhat)
