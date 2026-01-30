suppressPackageStartupMessages({
  library(BigVAR)
  library(Matrix)
  library(igraph)
})

# -----------------------------
# 1) Reading the file
# -----------------------------

X0 <- read.csv("PElectricity1344.csv", col.names = paste0("col",1:1344))
X0 <- as.matrix(X0)
X0 <- X0[1:100,1:50]

# -----------------------------
# 2) differences of log
# -----------------------------
# Security for very small values
eps <- 1e-6
Xlog <- log(pmax(X0, eps))
X <- diff(Xlog)  # Δlog
colnames(X) <- colnames(X0)
dimension=dim(X)
T=dimension[1]
k=dimension[2]

# -----------------------------
# 3) VAR(p) LASSO wit CV "Rolling" (temporal ordering)
# -----------------------------
p <- 2   #can be changed (eg: p=1,2,3)
T1 <- floor(T * 0.6)#Learning sample {1,...,T_1}
T2 <- floor(T * 0.8)#Test sample {T_1+1,...,T_2}

Model <- constructModel(
  Y = X, p = p,
  struct = "Basic",#means LASSO penalty
  gran = c(50, 10),#grid search for lambda,50 points+10 points for refinement
  cv = "Rolling",
  T1 = T1, T2 = T2,
  verbose = FALSE,#no message printed
  model.controls = list(intercept = TRUE)
)

fit <- cv.BigVAR(Model)
lambda_star <- fit@OptimalLambda
lambda_star
# betaPred: k x (1 + k*p) (Intercept + lags parameters)
beta <- fit@betaPred
rownames(beta) <- colnames(X)

cat("beta dim:", dim(beta), "\n")

# -----------------------------
# 4) Extraction of matrices A_j (k x k) for each lag
# -----------------------------
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

Ahat <- beta_to_Alist(beta, p = p)
#Important: these are the matrices corresponding to scaled data!!! 
#For non-scaled data, A should be replaced by D*A*D^{-1}, with D diagonal matrix
#with standard deviations of each time series.



# -----------------------------
# 5) Visualization 1: heatmap binary for zero vs non-zero coefficients
# -----------------------------
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

#------------------------------
#6) Extraction of residuals
errors<-fit@resids
acf(errors[,1])
acf(X[,1])

# -----------------------------
# 7) Visualization of dependencies (graph)
#--------------------------------

library(igraph)

# -----------------------------
# 1) full names of the 50 states
# -----------------------------
map_full_state_names <- function(abb) {
  m <- setNames(state.name, state.abb)
  m["DC"] <- "District of Columbia"
  ifelse(abb %in% names(m), unname(m[abb]), abb)
}

# -----------------------------
# 2) Extraction of edges+aggregation (max |coef|)
# -----------------------------
edges_from_A <- function(A_list, tol = 1e-10) {
  E <- data.frame()
  for (lag in seq_along(A_list)) {
    M <- A_list[[lag]]
    idx <- which(abs(M) > tol, arr.ind = TRUE)
    if (nrow(idx) == 0) next
    tmp <- data.frame(
      from = colnames(M)[idx[,2]],
      to   = rownames(M)[idx[,1]],
      lag  = lag,
      w    = M[idx]
    )
    E <- rbind(E, tmp)
  }
  E
}

build_state_graph <- function(Ahat,
                              tol = 1e-10,
                              top_m_cross = 220,
                              drop_loops = TRUE) {
  
  abb <- colnames(Ahat[[1]])
  Eall <- edges_from_A(Ahat, tol = tol)
  
  # Final weigths: max |coef| for all lags
  Eagg <- aggregate(abs(w) ~ from + to, data = Eall, FUN = max)
  colnames(Eagg)[3] <- "absw"
  
  if (drop_loops) Eagg <- subset(Eagg, from != to)
  
  # Only keep the strongest links (to avoid dense graphs)
  Eagg <- Eagg[order(-Eagg$absw), ]
  if (nrow(Eagg) > top_m_cross) Eagg <- Eagg[1:top_m_cross, ]
  
  g <- graph_from_data_frame(Eagg, directed = TRUE, vertices = data.frame(name = abb))
  V(g)$label <- map_full_state_names(V(g)$name)
  E(g)$absw  <- Eagg$absw
  g
}

layout_ultra_spread <- function(g,
                                niter = 7000,
                                charge = 0.10,
                                spring.length = 1100,
                                jitter_sd = 0.02,
                                seed = 123) {
  set.seed(seed)
  L <- layout_with_graphopt(g, niter = niter, charge = charge, spring.length = spring.length)
  L + matrix(rnorm(length(L), sd = jitter_sd), ncol = 2)
}

# -----------------------------
# 4) Export PNG
# -----------------------------
export_state_graph_png <- function(g, L,
                                   file = "state_dependency_graph.png",
                                   width = 3200, height = 2400, res = 220,
                                   label_cex = 0.9, label_dist = 1.35,
                                   curved = 0.28) {
  
  # Nodes
  V(g)$size <- 6 + 1.5 * degree(g, mode = "out")
  
  #Edges
  edge_width <- 0.6 + 3.0 * (E(g)$absw / max(E(g)$absw))
  
  #Arrows
  arrow_size  <- 0.30
  arrow_width <- 1.10
  
  png(filename = file, width = width, height = height, res = res)
  par(mar = c(1, 1, 3, 1))
  
  
  plot(
    g,
    layout = L,
    rescale = FALSE,
    xlim = range(L[,1]) * 1.18,
    ylim = range(L[,2]) * 1.18,
    vertex.label = V(g)$label,
    vertex.label.cex = label_cex,
    vertex.label.dist = label_dist,
    vertex.label.color = "red3",
    vertex.frame.color = "grey30",
    vertex.color = "white",
    edge.width = edge_width,
    edge.color = "grey45",
    edge.arrow.size = 0.5,  
    edge.curved = curved,
    # main = "Graph interactions for  VAR(2)-LASSO"
  )
  
  dev.off()
  message("PNG écrit dans : ", normalizePath(file, winslash = "/"))
}

# ============================================================
# EXECUTION
# ============================================================

g <- build_state_graph(Ahat, tol = 1e-10, top_m_cross = 220, drop_loops = TRUE)
L <- layout_ultra_spread(g, niter = 7000, charge = 0.10, spring.length = 1100, jitter_sd = 0.02)

export_state_graph_png(
  g, L,
  file = "state_dependency_graph.png",
  width = 3200, height = 2400, res = 220,
  label_cex = 0.9, label_dist = 1.35, curved = 0.28
)


# -----------------------------
# 8) Summary of "leaders" / "followers"
# -----------------------------
outdeg <- sort(degree(g, mode = "out"), decreasing = TRUE)
indeg  <- sort(degree(g, mode = "in"),  decreasing = TRUE)

cat("\nTop 10 out-degree (most influential series):\n")
print(head(outdeg, 10))

cat("\nTop 10 in-degree (most inluenced series):\n")
print(head(indeg, 10))


#---------------------------------
# 9) Comparison of 1-step predictions (LASSO-postLASSO-StandardVAR)
#---------------------------------

suppressPackageStartupMessages({
  library(vars)
})

# ---------- (de-)standardization of time series ----------
std_fit <- function(Ytrain) {
  mu  <- colMeans(Ytrain)
  sdv <- apply(Ytrain, 2, sd)
  sdv[sdv == 0] <- 1
  list(mu = mu, sd = sdv)
}
std_apply <- function(Y, mu, sdv) sweep(sweep(Y, 2, mu, "-"), 2, sdv, "/")
std_invert_point <- function(z, mu, sdv) as.numeric(mu + sdv * z)

# ---------- design VAR(p) for standardized data ----------
make_var_design <- function(Y, p, intercept = TRUE) {
  T <- nrow(Y); K <- ncol(Y)
  Yresp <- Y[(p+1):T, , drop = FALSE]
  Xlag  <- matrix(NA_real_, nrow = T - p, ncol = K*p)
  for (lag in 1:p) {
    Xlag[, ((lag-1)*K+1):(lag*K)] <- Y[(p+1-lag):(T-lag), , drop = FALSE]
  }
  Xbig <- if (intercept) cbind(Intercept = 1, Xlag) else Xlag
  list(Y = Yresp, X = Xbig)
}

# ---------- post-lasso refit OLS on the support given by BigVAR (standardized) ----------
post_lasso_refit_models <- function(betaPred, Ytrain_std, p, tol = 1e-12) {
  des <- make_var_design(Ytrain_std, p, intercept = TRUE)
  Yresp <- des$Y
  Xbig  <- des$X
  K <- ncol(Ytrain_std)
  
  models <- vector("list", K)
  for (i in 1:K) {
    sel <- which(abs(betaPred[i, ]) > tol)
    if (!(1 %in% sel)) sel <- c(1, sel) # keep intercept
    sel <- sort(unique(sel))
    fit <- lm.fit(x = Xbig[, sel, drop=FALSE], y = Yresp[, i])
    models[[i]] <- list(sel = sel, coef = fit$coefficients)
  }
  models
}

# ---------- reconstruction of betaPred_post (K x (1+Kp)) ----------
post_models_to_betaPred <- function(models, K, p) {
  P <- 1 + K*p
  B <- matrix(0, nrow = K, ncol = P)
  for (i in 1:K) B[i, models[[i]]$sel] <- models[[i]]$coef
  B
}

# ---------- 1-step forecast post-lasso (for standardized data) ----------
forecast_postlasso_1step_std <- function(beta_post, Ytrain_std, p) {
  K <- ncol(Ytrain_std)
  x <- 1
  for (lag in 1:p) x <- c(x, Ytrain_std[nrow(Ytrain_std) - (lag-1), ])
  yhat <- as.numeric(beta_post %*% x)  # (KxP) %*% (Px1)
  yhat
}

# ============================================================
# Rolling-origin compare
# ============================================================
rolling_compare_1step_original <- function(Y, p = 2, n_test = 36,
                                           bigvar_struct = "Basic",
                                           bigvar_gran = c(50, 10),
                                           cv = "Rolling",
                                           verbose = FALSE) {
  stopifnot(is.matrix(Y))
  T <- nrow(Y); K <- ncol(Y)
  stopifnot(T > n_test + p + 20)
  
  test_idx <- (T - n_test + 1):T
  
  pred_ols   <- matrix(NA_real_, n_test, K, dimnames = list(NULL, colnames(Y)))
  pred_lasso <- matrix(NA_real_, n_test, K, dimnames = list(NULL, colnames(Y)))
  pred_post  <- matrix(NA_real_, n_test, K, dimnames = list(NULL, colnames(Y)))
  y_true     <- matrix(NA_real_, n_test, K, dimnames = list(NULL, colnames(Y)))
  
  for (h in seq_along(test_idx)) {
    t0 <- test_idx[h]
    Ytrain <- Y[1:(t0-1), , drop=FALSE]
    y_next <- Y[t0, ]
    y_true[h, ] <- y_next
    
    # -------- VAR OLS (original scale, no standardization) --------
    fit_ols <- vars::VAR(as.data.frame(Ytrain), p = p, type = "const")
    fc_ols  <- predict(fit_ols, n.ahead = 1)
    pred_ols[h, ] <- sapply(colnames(Y), function(v) fc_ols$fcst[[v]][1, "fcst"])
    
    # -------- BigVAR LASSO --------
    st <- std_fit(Ytrain)
    
    TT <- nrow(Ytrain)
    T1 <- floor(0.60 * TT)
    T2 <- floor(0.80 * TT)
    
    Model <- BigVAR::constructModel(
      Y = Ytrain, p = p,
      struct = bigvar_struct,
      gran = bigvar_gran,
      cv = cv, T1 = T1, T2 = T2,
      verbose = verbose
    )
    fit_lasso <- BigVAR::cv.BigVAR(Model)
    
    # 1 step ahead prediction with BigVAR (data are scaled)
    fc_std <- BigVAR::predict(fit_lasso, n.ahead = 1)  # K x 1
    zhat <- as.numeric(fc_std[, 1])
    
    # de-standardize the prediction
    pred_lasso[h, ] <- std_invert_point(zhat, st$mu, st$sd)
    
    # -------- Post-LASSO --------
    # We need standardized data on train
    Ytrain_std <- std_apply(Ytrain, st$mu, st$sd)
    
    beta_std <- fit_lasso@betaPred
    models_post <- post_lasso_refit_models(beta_std, Ytrain_std, p = p)
    beta_post_std <- post_models_to_betaPred(models_post, K = K, p = p)
    
    zhat_post <- forecast_postlasso_1step_std(beta_post_std, Ytrain_std, p = p)
    pred_post[h, ] <- std_invert_point(zhat_post, st$mu, st$sd)
    
    if (h %% 5 == 0) cat("Rolling step", h, "/", n_test, "\n")
  }
  
  # -------- Performances --------
  mse_global <- function(pred) mean((y_true - pred)^2)
  mse_by_var <- function(pred) colMeans((y_true - pred)^2)
  
  out <- list(
    mse_global = c(
      VAR_OLS   = mse_global(pred_ols),
      VAR_LASSO = mse_global(pred_lasso),
      PostLASSO = mse_global(pred_post)
    ),
    mse_by_series = data.frame(
      series = colnames(Y),
      VAR_OLS   = mse_by_var(pred_ols),
      VAR_LASSO = mse_by_var(pred_lasso),
      PostLASSO = mse_by_var(pred_post)
    ),
    preds = list(ols = pred_ols, lasso = pred_lasso, post = pred_post),
    y_true = y_true,
    test_idx = test_idx
  )
  
  # Best improvements vs OLS
  out$mse_by_series$Delta_LASSO_minus_OLS <- out$mse_by_series$VAR_LASSO - out$mse_by_series$VAR_OLS
  out$mse_by_series$Delta_Post_minus_OLS  <- out$mse_by_series$PostLASSO - out$mse_by_series$VAR_OLS
  out$mse_by_series <- out$mse_by_series[order(out$mse_by_series$Delta_LASSO_minus_OLS), ]
  
  out
}


res <- rolling_compare_1step_original(X, p = 2, n_test = 36)
print(res$mse_global)
head(res$mse_by_series, 50)
