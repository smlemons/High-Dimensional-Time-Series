suppressPackageStartupMessages({
  library(BigVAR)
  library(vars)
  library(rrpack)
})

test <- data[,1:50]

## VAR (OLS)

VARselect(test, lag.max=6, type = "const")
#                   1            2            3            4            5
# AIC(n) 1.664884e+02 1.652258e+02 1.649224e+02 1.645570e+02 1.619101e+02
# HQ(n)  1.727930e+02 1.777115e+02 1.835890e+02 1.894047e+02 1.929388e+02
# SC(n)  1.828130e+02 1.975549e+02 2.132560e+02 2.288951e+02 2.422527e+02
# FPE(n) 2.043016e+72 6.282690e+71 5.829240e+71 6.389346e+71 9.926225e+70

#                  6
# AIC(n) 1.591686e+02
# HQ(n)  1.963783e+02
# SC(n)  2.555156e+02
# FPE(n) 2.204654e+70

# Hannan-Quinn and BIC (here called Schwarz Criterion) suggest p = 1

mod_var <- VAR(test, p = 1)
plot(fitted(mod_var), residuals(mod_var), pch = 16, col = rgb(0,0,0,0.1),
     xlab = "fitted values", ylab = "residuals", 
     main = "Residuals of VAR OLS Model with p = 1")

## BigVAR (LASSO)

T1 <- floor(nrow(test) * 0.6)#Learning sample {1,...,T_1}
T2 <- floor(nrow(test) * 0.8)#Test sample {T_1+1,...,T_2}

res <- NULL
lambda <- NULL

for(p in 1:5){
  Model <- constructModel(
    Y = test, p = p,
    struct = "Basic",# means LASSO penalty
    gran = c(50, 10),# grid search for lambda,50 points+10 points for refinement
    cv = "Rolling",
    T1 = T1, T2 = T2,
    verbose = FALSE,#no message printed
    model.controls = list(intercept = TRUE)
  )
  fit <- cv.BigVAR(Model)
  res[p] <- min(fit@OOSMSFE)
  lambda[p] <- fit@lambda_evolve_path[which.min(fit@OOSMSFE)]
}

print(res)
# [1] 1259.806 1305.848 1394.446 1418.317 1390.853
# we should take p=1

Model <- constructModel(
  Y = test, p = 1,
  struct = "Basic",# means LASSO penalty
  gran = c(50, 10),# grid search for lambda,50 points+10 points for refinement
  cv = "Rolling",
  T1 = T1, T2 = T2,
  verbose = FALSE,#no message printed
  model.controls = list(intercept = TRUE)
)
fit <- cv.BigVAR(Model)

plot(fit@fitted, fit@resids, pch = 16, col = rgb(0,0,0,0.1),
     xlab = "fitted values", ylab = "residuals", 
     main = "Residuals of BigVAR LASSO Model with p = 1")

beta <- fit@betaPred
rownames(beta) <- colnames(X)


beta_to_Alist <- function(beta, p) {
  K <- nrow(beta)
  stopifnot(ncol(beta) == 1 + K*p)
  A <- vector("list", p)
  for (lag in 1:p) {
    j0 <- 2 + (lag - 1) * K
    j1 <- 1 + lag * K
    A[[lag]] <- beta[, j0:j1, drop = FALSE]
    rownames(A[[lag]]) <- rownames(beta)   # target
    colnames(A[[lag]]) <- rownames(beta)   # source
  }
  A
}

Ahat <- beta_to_Alist(beta, p = 1)
#Important: these are the matrices corresponding to scaled data!!! 
#For non-scaled data, A should be replaced by D*A*D^{-1}, with D diagonal matrix
#with standard deviations of each time series.



# heatmap binary for zero vs non-zero coefficients


# Reset graphics 
graphics.off()

plot_sparsity_states_simple <- function(M, main = "", tol = 1e-10, cex_lab = 0.5) {
  B <- (abs(M) > tol) * 1
  K <- nrow(B)
  
  image(
    t(apply(B, 2, rev)),
    axes = FALSE,
    col = c("grey90", "red3"),      # 0 = gris, 1 = rouge
    breaks = c(-0.5, 0.5, 1.5),
    main = main
  )
  
  axis(1,
       at = seq(0, 1, length.out = K),
       labels = colnames(B),
       las = 2, cex.axis = cex_lab)
  
  axis(2,
       at = seq(0, 1, length.out = K),
       labels = rev(rownames(B)),
       las = 2, cex.axis = cex_lab)
  
  box()
}

# --------- PLOTS ----------
p <- length(Ahat)

par(
  mfrow = c(1, p),
  mar = c(8, 8, 4, 2)   # marges fixes, stables
)

for (lag in 1:p) {
  plot_sparsity_states_simple(
    Ahat[[lag]],
    main = paste0("A", lag, "  (grey = 0, red = non zero)"),
    cex_lab = 0.45
  )
}


## rrr (low rank)

make_var_design <- function(Y, p, intercept = TRUE) {
  TT <- nrow(Y); K <- ncol(Y)
  stopifnot(TT > p + 1)
  Yresp <- Y[(p+1):TT, , drop = FALSE]                # (TT-p) x K
  Xlag  <- matrix(NA_real_, nrow = TT - p, ncol = K*p)
  for (lag in 1:p) {
    Xlag[, ((lag-1)*K+1):(lag*K)] <- Y[(p+1-lag):(TT-lag), , drop = FALSE]
  }
  Xbig <- if (intercept) cbind(1, Xlag) else Xlag
  list(Y = Yresp, X = Xbig)
}


lags <- 1:5      # Search lags from 1 to 8

# Store results
results_grid <- expand.grid(lag = lags)
results_grid$BIC <- NA
results_grid$AIC <- NA
results_grid$HQ <- NA

for(p in 1:5){
  cat("Testing lag ", p)
  
  des <- make_var_design(test, p = p, intercept = TRUE)
  Yresp <- des$Y
  Chi   <- des$X
  
  # reduced-rank regression with adaptive nuclear norm (ann) + GCV
  fit <- rrr(Yresp, Chi, penaltySVD = "ann", ic.type = "GCV")
  
  # Fit reduced-rank regression
    r <- fit$r
    # Calculate information criteria
    n_obs <- nrow(Yresp)
    n_params <- K * r + r * K * p  # C matrix + D matrices
    
    # Residual variance
    resids <- Yresp -(Chi %*% fit$coef)
    resid_cov <- cov(resids)
    log_det_sigma <- log(det(resid_cov))
    
    # Information criteria
    results_grid$AIC[p] <- n_obs * K * log_det_sigma + 2 * n_params
    results_grid$BIC[p] <- n_obs * K * log_det_sigma + log(n_obs) * n_params
    results_grid$HQ[p]  <- n_obs * K * log_det_sigma + 2 * log(log(n_obs)) * n_params
}

results_grid # suggests p=5
# lag     BIC     AIC      HQ
# 1   1 5885099 5872281 5877230
# 2   2 5805036 5787188 5794079
# 3   3 5712304 5688513 5697699
# 4   4 5628008 5598279 5609759
# 5   5 5514366 5478702 5492475


