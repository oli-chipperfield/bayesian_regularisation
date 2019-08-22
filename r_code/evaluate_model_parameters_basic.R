source("data_prep_basic.R")

model.dump.string <- getwd()

#### Using full parameter set ####

#### Assess parameter estimates  ####

## Define function that'll summarise parameters

# sims <- blm.sims
# lambda <- FALSE
# sigma_sq <- TRUE
# freq.mod <- ols.mod

estimate_parameters <- function(sims, lambda = TRUE, sigma_sq = TRUE, freq.model = NULL) {
  
  b.mat <- sims$beta.matrix
  colnames(b.mat) <- colnames(sims$X)
  
  b.long <- data.table(melt(b.mat))
  names(b.long) <- c("sim", "var", "value")
  
  b.long[, var := as.character(var)]
  
  if(lambda == TRUE) {
    
    l.long <- data.table(sim = seq(1, length(sims$lambda.draws), 1), 
                         var = "lambda", 
                         value = sims$lambda.draws)
    
    b.long <- rbind(b.long, l.long)
    
  }
  
  if(sigma_sq == TRUE) {
    
    s.long <- data.table(sim = seq(1, length(sims$sigma.sq.draws), 1), 
                         var = "sigma_sq", 
                         value = sims$sigma.sq.draws)
    
    b.long <- rbind(b.long, s.long)
    
  }
  
  var.order <- unique(b.long$var)
  var.order <- var.order[order(var.order)]
  var.order <- var.order[!(var.order %in% c("sigma_sq", "lambda"))]
  
  if(sigma_sq == TRUE) {var.order <- c(var.order, "sigma_sq")}
  
  if(lambda == TRUE) {var.order <- c(var.order, "lambda")}
  
  b.long[, var := factor(var, levels = var.order)]
  
  #### Table summary
  
  b.summary <- b.long[, list(lower = quantile(value, 0.025),
                             average = mean(value),
                             upper = quantile(value, 0.975)),
                      by = "var"]
  
  #### Derive density graphs
  
  b.density <- list()
  
  for (i in 1:length(var.order)) {
    
    b.density[[i]] <- with(density(b.long[var == var.order[i],]$value), data.table(x, y))
    b.density[[i]]$var <- var.order[i]
    
  }  
  
  b.density <- rbindlist(b.density)
  b.density[, var := factor(var, levels = var.order)]
  
  b.density <- data.table(inner_join(b.density, b.summary, by = "var"))
  
  #### Estimate MAP
  
  maps <- inner_join(b.density, b.density[, list(y = max(y)), by = "var"], by = c("y", "var"))[, c("var", "x")]
  setnames(maps, "x", "map")
  
  b.summary <- inner_join(b.summary, maps, by = "var")
  
  ##### Quick hack #####
  
  b.summary$map <- b.summary$average
  
  ##### 
  
  # Add frequentist estimates
  
  b.summary$freq_point <- if(is.null(freq.model$coefficients)) {
    
    c(freq.model$beta[, 1], rep(NA, NROW(b.summary) - length(freq.model$beta[,1])))
    
    } else {
    
    c(freq.model$coefficients, rep(NA, NROW(b.summary) - length(freq.model$coefficients)))
      
    }
  
  #### Plot posterior densities
  
  b.density[, area_y := ifelse(x > lower & x < upper, y, 0)]
  
  b.plot <- ggplot(b.density, aes(x = x)) +  
    facet_wrap(~ var, scales = "free") +
    geom_line(aes(y = y)) +
    geom_area(aes(y = area_y), alpha = 0.3) +
    geom_vline(aes(xintercept = 0), colour = "red", linetype = "dashed") +
    geom_vline(data = b.summary, aes(xintercept = map), colour = "black", linetype = "dashed") +
    theme_bw(base_size = 12) +
    theme(strip.background = element_rect(fill = "white")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("Parameter value") +
    ylab("Density")
  
  # Add frequentist point estimates
  
  freq.points <- na.omit(b.summary[, c("var", "freq_point")])
  
  freq.points <- data.table(inner_join(freq.points, b.density, by = "var"))
  freq.points[, dist_point := abs(freq_point - x)]
  freq.points[, rank_dist := rank(dist_point), by = "var"]
  
  freq.points <- freq.points[rank_dist == 1,]

  b.plot <- b.plot + geom_point(data = freq.points, aes(x = freq_point, y = y), shape = 18, size = 3)
  
  # Add MAP estimates
  
  r.list <- list(par.summary = b.summary,
                 par.plot = b.plot,
                 par.long = b.long,
                 par.density = b.density)
  
}

## Define function that'll summarise parameter covariance

parameter_covar <- function(sims, lambda = TRUE, sigma_sq = TRUE) {
  
  b.mat <- sims$beta.matrix
  colnames(b.mat) <- colnames(sims$X)
  
  b.mat <- data.table(b.mat)
  setcolorder(b.mat, names(b.mat)[order(names(b.mat))])
  
  if(sigma_sq == TRUE) {
    
    b.mat$sigma_sq <- sims$sigma.sq.draws
    
  }
  
  if(lambda == TRUE) {
    
    b.mat$lambda <- sims$lambda.draws
    
  }
  
  p.plot <- ggpairs(b.mat, 
                    upper = list(continuous = wrap("density", colour = "black", bins = 5)),
                    lower = list(continuous = "blank"),
                    diag = list(continuous = "blankDiag"))
  
  p.plot <- p.plot + theme_bw(base_size = 10)
  p.plot <- p.plot + theme(strip.background = element_rect(fill = "white"))
  
  ### Plot individually if you need a close up
  
  par.list <- colnames(b.mat)
  
  par.combs <- data.table(expand.grid(par.list, par.list))
  par.combs <- par.combs[Var1 != Var2,]
  
  
  p.plots <- vector("list",nrow(par.combs))
  
  names(p.plots) <- paste0(par.combs$Var1,"_", par.combs$Var2)
  
  for (i in 1:nrow(par.combs)) {
    
    d <- cbind(b.mat[, par.combs[i, 1]$Var1, with = FALSE],
               b.mat[, par.combs[i, 2]$Var2, with = FALSE])
    
    setnames(d, c("x", "y"))
    
    d.plot <- ggplot(d, aes(x = x, y = y)) + geom_density_2d(bins = 15, colour = "black") +
      theme_bw(base_size = 12) + 
      xlab(par.combs[i, 1]$Var1) + ylab(par.combs[i, 2]$Var2)
    
    p.plots[[i]] <- d.plot
    
  }
  
  return(list(p.plot = p.plot,
              p.plots = p.plots))
  
}

#### Load models

ols.mod <- readRDS(paste0(model.dump.string, "ols_model_basic.RDS"))

freq.ridge.mod <- readRDS(paste0(model.dump.string, "freq_ridge_model_basic.RDS"))

freq.lasso.mod <- readRDS(paste0(model.dump.string, "freq_lasso_model_basic.RDS"))

blm.sims <- readRDS(paste0(model.dump.string, "blm_sims_basic.RDS"))

ridge.sims <- readRDS(paste0(model.dump.string, "ridge_sims_basic.RDS"))

lasso.sims <- readRDS(paste0(model.dump.string, "lasso_sims_basic.RDS"))

## BLM ##

blm.parameters <- estimate_parameters(blm.sims, lambda = FALSE, sigma_sq = TRUE, freq.model = ols.mod)

# blm.parameters$par.summary
#
#blm.parameters$par.plot

### Look at the covariates

blm.covar <- parameter_covar(blm.sims, lambda = FALSE, sigma_sq = TRUE)

#blm.covar$p.plot

## Ridge ##

ridge.parameters <- estimate_parameters(ridge.sims, lambda = TRUE, sigma_sq = TRUE, freq.model = freq.ridge.mod)

#ridge.parameters$par.summary
# 
#ridge.parameters$par.plot

### Look at the covariates

ridge.covar <- parameter_covar(ridge.sims, lambda = TRUE, sigma_sq = TRUE)

#ridge.covar$p.plot

#### Lasso

lasso.parameters <- estimate_parameters(lasso.sims, lambda = TRUE, sigma_sq = TRUE, freq.model = freq.lasso.mod)

#lasso.parameters$par.summary
 
#lasso.parameters$par.plot

### Look at the covariates

lasso.covar <- parameter_covar(lasso.sims, lambda = TRUE, sigma_sq = TRUE)

lasso.covar$p.plot

#### Compare the model parameters

## Compare parameters

compare_parameters <- function(blm.parameters,
                               ridge.parameters,
                               lasso.parameters,
                               ols.mod,
                               freq.ridge.mod,
                               freq.lasso.mod) {
  
  blm.parameters$par.long[, model := "Linear"]
  
  ridge.parameters$par.long[, model := "Ridge"]
  
  lasso.parameters$par.long[, model := "LASSO"]
  
  par.summaries <- rbind(blm.parameters$par.long,
                         ridge.parameters$par.long,
                         lasso.parameters$par.long)
  
  par.summaries[, model := factor(model, levels = c("Linear", "Ridge", "LASSO"))]
  
  ### Get frequentist point estimates
  
  ols.point <- data.table(var = unique(par.summaries$var))
  ols.point <- ols.point[var != "sigma_sq" & var != "lambda",]
  ols.point$point_est <- ols.mod$coefficients
  ols.point[, model := "Linear"]
  
  ridge.point <- data.table(var = unique(par.summaries$var))
  ridge.point <- ridge.point[var != "sigma_sq" & var != "lambda",]
  ridge.point$point_est <- freq.ridge.mod$beta[,1]
  ridge.point[, model := "Ridge"]
  
  lasso.point <- data.table(var = unique(par.summaries$var))
  lasso.point <- lasso.point[var != "sigma_sq" & var != "lambda",]
  lasso.point$point_est <- freq.lasso.mod$beta[,1]
  lasso.point[, model := "Lasso"]
  
  point.tab <- rbind(ols.point, ridge.point, lasso.point)
  point.tab[, model := factor(model, levels = c("Linear", "Ridge", "LASSO"))]
  
  ### Box plots
  
  pb <- ggplot(par.summaries, aes(x = model, y = value)) + geom_boxplot(outlier.size = 0.1) +
    facet_wrap(~ var, scales = "free_y") +
    geom_hline(aes(yintercept = 0), colour = "red", linetype = "dashed") +
    theme_bw(base_size = 12) +
    theme(strip.background = element_rect(fill = "white")) +
    xlab("Model") +
    ylab("Posterior density")
  
  # Add point ests
  
  pb <- pb + geom_point(data = point.tab, aes(x = model, y = point_est), shape = 4, size = 3)
  
  ### Densities
  
  # b.density <- list()
  # 
  # for (i in 1:length(var.order)) {
  #   
  #   b.density[[i]] <- with(density(b.long[var == var.order[i],]$value), data.table(x, y))
  #   b.density[[i]]$var <- var.order[i]
  #   
  # }  
  # 
  # b.density <- rbindlist(b.density)
  # b.density[, var := factor(var, levels = var.order)]
  # 
  # b.density <- data.table(inner_join(b.density, b.summary, by = "var"))
  # 
  b.density <- list()
  
  var.mod <- unique(par.summaries[, c("var", "model")])
  
  for (i in 1:nrow(var.mod)) {
    
    b.data <- par.summaries[var == var.mod$var[i] & model == var.mod$model[i],]
    
    b.density[[i]] <- with(density(b.data$value), data.table(x, y))
    b.density[[i]]$var <- var.mod$var[i]
    b.density[[i]]$model <- var.mod$model[i]
    
  }
  
  b.density <- rbindlist(b.density)
  
  p <- ggplot(b.density, aes(x = x, y = y, linetype = model)) +
    geom_line() +
    facet_wrap(~ var, scales = "free", nrow = 3) +
    geom_vline(aes(xintercept = 0), colour = "red", linetype = "dashed") +
    scale_linetype_manual(name = "Model", values=c("solid", "dashed", "dotted")) +
    theme_bw(base_size = 12) +
    theme(strip.background = element_rect(fill = "white")) +
    theme(legend.position = "bottom") +
    xlab("Value") +
    ylab("Posterior density")    
  
  # Add freq point estimates
  
  point.den <- data.table(inner_join(point.tab, b.density, by = c("var", "model")))
  
  point.den[, dist_point := abs(point_est - x)]
  point.den[, rank_dist := rank(dist_point), by = c("var", "model")]
  
  point.den <- point.den[rank_dist == 1,]
  
  
  
  # p <- p + geom_point(data = point.den, aes(x = point_est, y = y, shape = model)) +
  #         scale_shape_manual(name = "Model", values = c(3, 4, 8))
  # 
  # Add point ests
  
  #p <- p + geom_vline(data = point.tab, aes(xintercept = point_est, linetype = model))
  
  ### MAP estimates
  
  blm.parameters$par.summary$model <- "Linear"
  
  ridge.parameters$par.summary$model <- "Ridge"
  
  lasso.parameters$par.summary$model <- "Lasso"
  
  par.maps <- rbind(blm.parameters$par.summary,
                    ridge.parameters$par.summary,
                    lasso.parameters$par.summary)
  
  par.maps$model <- factor(par.maps$model, levels = c("Linear", "Ridge", "Lasso"))
  
  p.map <- ggplot(par.maps, aes(x = model, y = map)) + 
    geom_bar(stat = "identity", alpha = 0, colour = "black") + 
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.5) +
    facet_wrap(~ var, scales = "free") +
    theme_bw(base_size = 12) +
    theme(strip.background = element_rect(fill = "white")) +
    geom_hline(aes(yintercept = 0), colour = "red", linetype = "dashed") +
    xlab("Model") + ylab("MAP")
  
  #### Which overlap with zero
  
  par.maps$includes_zero <- data.table::between(0, par.maps$lower, par.maps$upper)
  
  inc.zero <- spread(par.maps[, c("var", "model", "includes_zero")], model,  includes_zero)
  
  r.list <- list(box.p = pb,
                 den.p = p,
                 map.p = p.map,
                 inc.zero = inc.zero)
  
  return(r.list) 
  
}

comp.par <- compare_parameters(blm.parameters,
                               ridge.parameters,
                               lasso.parameters,
                               ols.mod,
                               freq.ridge.mod,
                               freq.lasso.mod)

# comp.par$box.p
# 
# comp.par$den.p

#comp.par$inc.zero






