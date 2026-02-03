# ============================================================
# VAR(2) with Graphical LASSO on residuals
# ============================================================

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

test <- data

# ----------------------------------
harmonic <- function(x) x
# 2) VAR(2) BigVAR
# ----------------------------------
p <- 1
T1 <- floor(0.60 * T)
T2 <- floor(0.80 * T)

Model <- constructModel(
  Y = test, p = p,
  struct = "Basic",
  gran = c(50, 10),
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

mu  <- colMeans(test)
sdv <- apply(test, 2, sd); sdv[sdv == 0] <- 1
Z   <- sweep(sweep(test, 2, mu, "-"), 2, sdv, "/")

des   <- make_var_design(Z, p = p, intercept = TRUE)
Zresp <- des$Y
Xbig  <- des$X

Zhat  <- Xbig %*% t(beta)
U_std <- Zresp - Zhat
U_orig <- sweep(U_std, 2, sdv, "*")
colnames(U_orig) <- colnames(test)
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
colnames(Theta_std) <- rownames(Theta_std) <- colnames(test)
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
colnames(P) <- rownames(P) <- colnames(test)




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

plot_heatmap(P, main = "Partial Correlations (Graphical LASSO)")

