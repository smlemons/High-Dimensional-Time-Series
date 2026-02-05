

# compute CV for every p 

Model <- constructModel(
  Y = data, p = 1,
  struct = "Basic",# means LASSO penalty
  gran = c(300, 10),# grid search for lambda,50 points+10 points for refinement
  cv = "Rolling",
  T1 = T1, T2 = T2,
  verbose = FALSE,#no message printed
  model.controls = list(intercept = TRUE)
)
fit_test1 <- cv.BigVAR(Model)

Model <- constructModel(
  Y = data, p = 2,
  struct = "Basic",# means LASSO penalty
  gran = c(300, 10),# grid search for lambda,50 points+10 points for refinement
  cv = "Rolling",
  T1 = T1, T2 = T2,
  verbose = FALSE,#no message printed
  model.controls = list(intercept = TRUE)
)
fit_test2 <- cv.BigVAR(Model)


Model <- constructModel(
  Y = data, p = 3,
  struct = "Basic",# means LASSO penalty
  gran = c(300, 10),# grid search for lambda,50 points+10 points for refinement
  cv = "Rolling",
  T1 = T1, T2 = T2,
  verbose = FALSE,#no message printed
  model.controls = list(intercept = TRUE)
)
fit_test3 <- cv.BigVAR(Model)


Model <- constructModel(
  Y = data, p = 4,
  struct = "Basic",# means LASSO penalty
  gran = c(300, 10),# grid search for lambda,50 points+10 points for refinement
  cv = "Rolling",
  T1 = T1, T2 = T2,
  verbose = FALSE,#no message printed
  model.controls = list(intercept = TRUE)
)
fit_test4 <- cv.BigVAR(Model)

Model <- constructModel(
  Y = data, p = 5,
  struct = "Basic",# means LASSO penalty
  gran = c(300, 10),# grid search for lambda,50 points+10 points for refinement
  cv = "Rolling",
  T1 = T1, T2 = T2,
  verbose = FALSE,#no message printed
  model.controls = list(intercept = TRUE)
)
fit_test5 <- cv.BigVAR(Model)


## select p based on proportinal improvement

improv <- matrix(NA, ncol = 149, nrow = 5)
improv[1, ] <- (fit_test1@MeanMSFE - fit_test1@OOSMSFE) / fit_test1@MeanMSFE
improv[2, ] <- (fit_test2@MeanMSFE - fit_test2@OOSMSFE) / fit_test2@MeanMSFE
improv[3, ] <- (fit_test3@MeanMSFE - fit_test3@OOSMSFE) / fit_test3@MeanMSFE
improv[4, ] <- (fit_test4@MeanMSFE - fit_test4@OOSMSFE) / fit_test4@MeanMSFE
improv[5, ] <- (fit_test5@MeanMSFE - fit_test5@OOSMSFE) / fit_test5@MeanMSFE
improv_res <- apply(improv, 1, mean)
print(improv_res)

p <- which.max(improv_res)


#pdf("reds_lasso.pdf", height = 5, width = 8)
#plot(fit_test1@fitted, fit_test1@resids, pch = 16, col = rgb(0,0,0,0.1),
#     xlab = "fitted values", ylab = "residuals", 
#     main = "Residuals of BigVAR LASSO Model with p = 1")
#dev.off()

#pdf("reds_lasso2.pdf", height = 5, width = 8)
#plot_average_acf_pacf(fit_test1@resids)
#dev.off()

#beta <- fit@betaPred
#rownames(beta) <- colnames(X)

# beta_to_Alist <- function(beta, p) {
#   K <- nrow(beta)
#   stopifnot(ncol(beta) == 1 + K*p)
#   A <- vector("list", p)
#   for (lag in 1:p) {
#     j0 <- 2 + (lag - 1) * K
#     j1 <- 1 + lag * K
#     A[[lag]] <- beta[, j0:j1, drop = FALSE]
#     rownames(A[[lag]]) <- rownames(beta)   # target
#     colnames(A[[lag]]) <- rownames(beta)   # source
#   }
#   A
# }
# 
# Ahat <- beta_to_Alist(beta, p = 1)
# #Important: these are the matrices corresponding to scaled data!!!
# #For non-scaled data, A should be replaced by D*A*D^{-1}, with D diagonal matrix
# #with standard deviations of each time series.
# 
# 
# 
# # heatmap binary for zero vs non-zero coefficients
# 
# 
# # Reset graphics
# graphics.off()
# 
# plot_sparsity_states_simple <- function(M, main = "", tol = 1e-10, cex_lab = 0.5) {
#   B <- (abs(M) > tol) * 1
#   K <- nrow(B)
# 
#   image(
#     t(apply(B, 2, rev)),
#     axes = FALSE,
#     col = c("grey90", "red3"),      # 0 = gris, 1 = rouge
#     breaks = c(-0.5, 0.5, 1.5),
#     main = main
#   )
# 
#   axis(1,
#        at = seq(0, 1, length.out = K),
#        labels = colnames(B),
#        las = 2, cex.axis = cex_lab)
# 
#   axis(2,
#        at = seq(0, 1, length.out = K),
#        labels = rev(rownames(B)),
#        las = 2, cex.axis = cex_lab)
# 
#   box()
# }
# 
# # --------- PLOTS ----------
# p <- length(Ahat)
# 
# par(
#   mfrow = c(1, p),
#   mar = c(8, 8, 4, 2)   # marges fixes, stables
# )
# 
# for (lag in 1:p) {
#   pdf("Heatmap_Coefficients_Lasso.pdf", height = 5, width = 8)
#   plot_sparsity_states_simple(
#     Ahat[[lag]],
#     main = paste0("A", lag, "  (grey = 0, red = non zero)"),
#     cex_lab = 0.45
#   )
#   dev.off()
# }

suppressPackageStartupMessages({
  library(BigVAR)
  library(vars)
  library(rrpack)
})

#air <- read.csv("TaiwanAirBox032017.csv", col.names = c("time", paste0("col",2:517)))
#data <- diff(ts(air[, -1]), lag = 24)
#data <- data[, -c(29, 70)]
data <- X
# -----------------------------
# Helpers
# -----------------------------
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

betaPred_to_Alist <- function(betaPred, k, p, has_intercept = TRUE) {
  A <- vector("list", p)
  off <- if (has_intercept) 1 else 0
  for (lag in 1:p) {
    j0 <- off + (lag - 1) * k + 1
    j1 <- off + lag * k
    A[[lag]] <- betaPred[, j0:j1, drop = FALSE]
  }
  A
}

destandardize_Alist <- function(A_std, sdv) {
  D <- diag(sdv); Dinv <- diag(1/sdv)
  lapply(A_std, function(M) D %*% M %*% Dinv)
}

compute_c_from_muA <- function(mu, A_list) {
  A_sum <- Reduce("+", A_list)
  as.numeric(mu - A_sum %*% mu)
}

forecast_1step_VAR <- function(Ytrain, A_list, c_vec, p) {
  # yhat_{t+1} = c + sum_{l=1..p} A_l y_{t+1-l}
  yhat <- c_vec
  for (lag in 1:p) yhat <- yhat + A_list[[lag]] %*% Ytrain[nrow(Ytrain) - (lag-1), ]
  as.numeric(yhat)
}

mse_global <- function(err) mean(err^2)
mse_by_series <- function(err) colMeans(err^2)

# -----------------------------
# 1) BigVAR-LASSO fit -> A(orig), c(orig), forecast (1-step)
# -----------------------------
fit_bigvar_forecast1 <- function(Ytrain, p, gran = c(50,10), cv = "Rolling") {
  mu  <- colMeans(Ytrain)
  sdv <- apply(Ytrain, 2, sd); sdv[sdv == 0] <- 1
  
  TT <- nrow(Ytrain)
  # T1/T2 only used for CV splitting inside BigVAR; keep them feasible
  T1 <- floor(0.60 * TT)
  T2 <- floor(0.80 * TT)
  
  Model <- constructModel(
    Y = Ytrain, p = p, struct = "Basic",
    gran = gran, cv = cv, T1 = T1, T2 = T2,
    verbose = FALSE,
    model.controls = list(intercept = TRUE)
  )
  fit <- cv.BigVAR(Model)
  
  beta_std <- fit@betaPred                  # k x (1+k*p)
  A_std <- betaPred_to_Alist(beta_std, k = ncol(Ytrain), p = p, has_intercept = TRUE)
  A_org <- destandardize_Alist(A_std, sdv)
  c_org <- compute_c_from_muA(mu, A_org)
  
  yhat <- forecast_1step_VAR(Ytrain, A_org, c_org, p)
  list(yhat = yhat, A1 = A_org[[1]])
}

# -----------------------------
# 2) Low-rank VAR via rrpack::rrr (global low-rank on stacked lag matrices)
#    Y = Chi * B + e,  B is (1+k*p) x k, low-rank encouraged via penaltySVD="ann"
# -----------------------------
fit_rrr_forecast1 <- function(Ytrain, p) {
  mu <- colMeans(Ytrain)
  
  des <- make_var_design(Ytrain, p = p, intercept = TRUE)
  Yresp <- des$Y
  Chi   <- des$X
  
  # reduced-rank regression with adaptive nuclear norm (ann) + GCV
  LRE <- rrr(Yresp, Chi, penaltySVD = "ann", ic.type = "GCV")
  B <- LRE$coef             # (1+k*p) x k
  
  # Extract A_l (k x k) and c
  c_hat <- as.numeric(B[1, ])
  A_list <- vector("list", p)
  for (lag in 1:p) {
    j0 <- 2 + (lag - 1)*ncol(Ytrain)
    j1 <- 1 + lag*ncol(Ytrain)
    A_list[[lag]] <- t(B[j0:j1, , drop = FALSE])  # k x k
  }
  
  # 1-step forecast using estimated (c, A_l)
  yhat <- forecast_1step_VAR(Ytrain, A_list, c_hat, p)
  list(yhat = yhat, A1 = A_list[[1]], rank = LRE$rank)
}

# -----------------------------
# 3) VAR OLS via vars::VAR
# -----------------------------
fit_varols_forecast1 <- function(Ytrain, p) {
  fit <- vars::VAR(as.data.frame(Ytrain), p = p, type = "const")
  # matrices A_l
  Alist <- Acoef(fit)            # list length p, each k x k
  # intercept vector (const per equation)
  c_hat <- sapply(coef(fit), function(m) m["const", "Estimate"])
  yhat <- forecast_1step_VAR(Ytrain, Alist, c_hat, p)
  list(yhat = yhat, A1 = Alist[[1]])
}





# ============================================================
# Rolling evaluation (1-step) — p is a parameter
# ============================================================
# test <- data[,1:50]

k <- ncol(data)

t_start <- floor(0.7 * nrow(data))        # <-- align with your old code style (ensure > p)
t_end   <- nrow(data) -100

pred_big  <- matrix(NA_real_, nrow = t_end - t_start + 1, ncol = k)
pred_rrr  <- pred_big
pred_ols  <- pred_big
y_true    <- pred_big
rank_rrr  <- integer(nrow(pred_big))

colnames(pred_big) <- colnames(pred_rrr) <- colnames(pred_ols) <- colnames(y_true) <- colnames(data) 

idx <- 1
for (t in t_start:t_end) {
  Ytrain <- data[1:t, , drop = FALSE]
  y_next <- data[t+1, ]
  y_true[idx, ] <- y_next
  
  # BigVAR-LASSO (attention: internal standardization -> we de-standardize A and forecast ourselves)
  fb <- fit_bigvar_forecast1(Ytrain, p = p)
  pred_big[idx, ] <- fb$yhat
  
  # Low-rank (global on stacked lags)
  fr <- fit_rrr_forecast1(Ytrain, p = 3)
  pred_rrr[idx,] <- fr$yhat
  rank_rrr[idx] <- fr$rank
  
  # VAR OLS
  fo <- fit_varols_forecast1(Ytrain, p = 1)
  pred_ols[idx, ] <- fo$yhat
  
  if (idx %% 10 == 0) cat("t =", t, "done\n")
  idx <- idx + 1
}

# ============================================================
# Performance: global and per series
# ============================================================

err_big <- y_true - pred_big
err_rrr <- y_true - pred_rrr
err_ols <- y_true - pred_ols

cat("\n=== Global MSE (1-step) ===\n")
print(c(
  BigVAR_LASSO = mse_global(err_big),
  LowRank_RRR  = mse_global(err_rrr),
  VAR_OLS      = mse_global(err_ols)
))
res_comp <- c(BigVAR_LASSO = mse_global(err_big),LowRank_RRR  = mse_global(err_rrr),
  VAR_OLS = mse_global(err_ols))

suppressPackageStartupMessages({
  library(BigVAR)
  library(glasso)
  library(igraph)
  library(maps)
  library(fields)
})

# -----------------------------
# 1) Data
# -----------------------------

#test <- data
#air <- read.csv("TaiwanAirBox032017.csv", col.names = c("time", paste0("col",2:517)))
#data <- diff(ts(air[, -1]), lag = 24)
#data <- data[, -c(29, 70)]
data <- X

# ----------------------------------
harmonic <- function(x) x
# 2) VAR(2) BigVAR
# ----------------------------------

TT <- nrow(data)
T1 <- floor(0.60 * TT)
T2 <- floor(0.80 * TT)

Model <- constructModel(
  Y = data, p = p,
  struct = "Basic",
  gran = c(300, 10),
  cv = "Rolling",
  T1 = T1, T2 = T2,
  verbose = FALSE,
  model.controls = list(intercept = TRUE)
)

fit <- cv.BigVAR(Model)
beta <- fit@betaPred

# -----------------------------
# 3) Residuals extraction + original scale
# -----------------------------
make_var_design <- function(Y, p, intercept = TRUE) {
  T <- nrow(Y); k <- ncol(Y)
  Yresp <- Y[(p+1):T, , drop = FALSE]
  Xlag  <- matrix(NA_real_, nrow = T - p, ncol = k*p)
  for (lag in 1:p) {
    Xlag[, ((lag-1)*k+1):(lag*k)] <- Y[(p+1-lag):(T-lag), , drop = FALSE]
  }
  Xbig <- if (intercept) cbind(1, Xlag) else Xlag
  list(Y = Yresp, X = Xbig)
}

mu  <- colMeans(data)
sdv <- apply(data, 2, sd); sdv[sdv == 0] <- 1
Z   <- sweep(sweep(data, 2, mu, "-"), 2, sdv, "/")

des   <- make_var_design(Z, p = p, intercept = TRUE)
Zresp <- des$Y
Xbig  <- des$X

Zhat  <- Xbig %*% t(beta)
U_std <- Zresp - Zhat
U_orig <- sweep(U_std, 2, sdv, "*")
colnames(U_orig) <- colnames(data)
cat("Residuals dim =", dim(U_orig), "\n")

# -----------------------------
# 4) Graphical LASSO standardized residuals
# -----------------------------
U_sdv  <- apply(U_orig, 2, sd); U_sdv[U_sdv == 0] <- 1
U_std2 <- sweep(U_orig, 2, U_sdv, "/")

log_likelihood <- function(precision, emp_cov) {
  p <- nrow(precision)
  logdet <- determinant(precision, logarithm = TRUE)$modulus
  0.5 * (as.numeric(logdet) - sum(emp_cov * precision) - p * log(2*pi))
}

glasso_cv_time <- function(ts, rholist = seq(0.01, 1.5, 0.01)) {
  TT <- nrow(ts)
  cut <- floor(2*TT/3)
  S_train <- cov(ts[1:cut, , drop=FALSE])
  S_test  <- cov(ts[(cut+1):TT, , drop=FALSE])
  
  ll <- numeric(length(rholist))
  for (i in seq_along(rholist)) {
    G <- glasso(S_train, rho = rholist[i], trace = FALSE)
    ll[i] <- log_likelihood(G$wi, S_test)
  }
  rho_best <- rholist[which.max(ll)]
  Gfull <- glasso(cov(ts), rho = rho_best, trace = FALSE)
  list(Theta = Gfull$wi, rho = rho_best)
}

gl <- glasso_cv_time(U_std2)
Theta_std <- gl$Theta
colnames(Theta_std) <- rownames(Theta_std) <- colnames(data)
cat("Selected rho (glasso) =", gl$rho, "\n")

# -----------------------------
# 5) Convert to partial correlations
# -----------------------------
partial_corr <- function(Theta) {
  d <- sqrt(diag(Theta))
  P <- -Theta / (d %*% t(d))
  diag(P) <- 0
  P
}

P <- partial_corr(Theta_std)
colnames(P) <- rownames(P) <- colnames(data)




#---------------------------------------------------------------
# Heatmap for Census categories of States
#---------------------------------------------------------------
plot_heatmap <- function(
    P,
    main = "Partial correlations (Census category)"
) {
  
  pal <- colorRampPalette(
    c("#2166ac", "#f7f7f7", "#b2182b")
  )(101)
  
  a <- max(abs(P), na.rm = TRUE)
  
  op <- par(no.readonly = TRUE)
  on.exit(par(op), add = TRUE)
  
  # extra right margin for legend
  par(mar = c(6, 6, 4, 9))
  
  image(
    1:nrow(P), 1:ncol(P),
    t(P[nrow(P):1, ]),
    col = pal,
    zlim = c(-a, a),
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = main
  )
  
  axis(
    1,
    at = 1:nrow(P),
    labels = paste0("TS ", 1:nrow(P)),
    las = 2,
    cex.axis = .6
  )
  
  axis(
    2,
    at = 1:ncol(P),
    labels = rev(paste0("TS ", 1:nrow(P))),
    las = 2,
    cex.axis = .6
  )
  
  box()
  
  # legend
  fields::image.plot(
    zlim = c(-a, a),
    legend.only = TRUE,
    col = pal,
    legend.lab = "\n\n Partial correlation",
    add = TRUE,
    legend.line=3
  )
}

pdf("Heatmap_Partial_Correlation.pdf", height = 5, width = 8)
plot_heatmap(P, main = "Partial Correlations (Graphical LASSO)")
dev.off()

