suppressPackageStartupMessages({
  library(BigVAR)
  library(vars)
  library(rrpack)
})


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
test <- data[,1:50]

p <- 1                 # <-- change number of lags here
k <- ncol(test)

t_start <- floor(0.7 * nrow(test))        # <-- align with your old code style (ensure > p)
t_end   <- nrow(test) -100

pred_big  <- matrix(NA_real_, nrow = t_end - t_start + 1, ncol = k)
pred_rrr  <- pred_big
pred_ols  <- pred_big
y_true    <- pred_big
rank_rrr  <- integer(nrow(pred_big))

colnames(pred_big) <- colnames(pred_rrr) <- colnames(pred_ols) <- colnames(y_true) <- colnames(test)

idx <- 1
for (t in t_start:t_end) {
  Ytrain <- test[1:t, , drop = FALSE]
  y_next <- test[t+1, ]
  y_true[idx, ] <- y_next
  
  # BigVAR-LASSO (attention: internal standardization -> we de-standardize A and forecast ourselves)
  fb <- fit_bigvar_forecast1(Ytrain, p = p)
  pred_big[idx, ] <- fb$yhat
  
  # Low-rank (global on stacked lags)
  fr <- fit_rrr_forecast1(Ytrain, p = p)
  pred_rrr[idx, ] <- fr$yhat
  rank_rrr[idx] <- fr$rank
  
  # VAR OLS
  fo <- fit_varols_forecast1(Ytrain, p = p)
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

# BigVAR_LASSO  LowRank_RRR      VAR_OLS 
# 55.02024     53.09655     52.73878 