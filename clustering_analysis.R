############################################################
# Taiwan AirBox — Clustering benchmark (K-means + HC)
# Methods:
# 1) K-means on ACF coefficients
# 2) K-means on time-series vectors (Euclidean)
# 3) HC on ACF FEATURE distance (Euclidean on ACF coeffs)
# 4) HC on GCC distance
# 5) HC on CC distance
# 6) HC on AR distance (AR.PIC)
# 7) HC on ACF DISSIMILARITY (TSclust diss METHOD="ACF")
#
# Model selection:
# - For all methods: choose K by maximum average silhouette over K=2..Kmax
# - For K-means methods (1) and (2): also compute CH index as an option
#
# Requirements:
# - Your saved file: airbox_stl_deseasonalized.rds (contains X_adj, datetime, sensor_table)
# - diss() is available in the global env (as in your setup; do NOT use TSclust::)
# - GCCmatrix() is available in the global env (as in your setup)
############################################################



download_and_untar <- function(url, name) {
  tf <- tempfile(fileext = ".tar.gz")
  download.file(url, tf, mode = "wb")
  out <- file.path(tempdir(), paste0(name, "_src"))
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  untar(tf, exdir = out)
  
  pkgdir <- file.path(out, list.dirs(out, full.names = FALSE, recursive = FALSE)[1])
  pkgdir
}


source_functions_from_pkg <- function(pkgdir, fun_names) {
  rdir <- file.path(pkgdir, "R")
  rfiles <- list.files(rdir, pattern = "\\.R$", full.names = TRUE)
  
  needed <- character()
  for (f in rfiles) {
    txt <- paste(readLines(f, warn = FALSE), collapse = "\n")
    for (fn in fun_names) {
      # cherche "fn <- function" ou "fn=function"
      if (grepl(paste0("\\b", fn, "\\b\\s*<-\\s*function"), txt) ||
          grepl(paste0("\\b", fn, "\\b\\s*=\\s*function"), txt)) {
        needed <- c(needed, f)
      }
    }
  }
  needed <- unique(needed)
  if (length(needed) == 0) stop("Aucun fichier R ne contient les fonctions demandées.")
  
  for (f in needed) {
    source(f, local = .GlobalEnv)
  }
  invisible(needed)
}

## TSclust & SLBDD from CRAN
tsclust_dir <- download_and_untar(
  "https://cran.r-project.org/src/contrib/Archive/TSclust/TSclust_1.3.2.tar.gz",
  "TSclust"
)

slbdd_dir <- download_and_untar(
  "https://cran.r-project.org/src/contrib/Archive/SLBDD/SLBDD_0.0.4.tar.gz",
  "SLBDD"
)


## - TSclust : diss() + méthodes "ACF" and "AR.PIC"
## - SLBDD : GCCmatrix() + silh.clus()
source_functions_from_pkg(tsclust_dir, c("diss", "diss.ACF", "diss.AR.PIC"))
source_functions_from_pkg(slbdd_dir,  c("GCCmatrix", "silh.clus"))


stopifnot(exists("diss"), exists("GCCmatrix"), exists("silh.clus"))

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(maps)
  library(gridExtra)
  library(cluster)
  library(stats)
  library(dtw)     # for DTW distances; dtwclust used if installed (optional)
})

# ==========================================================
# 0) Load data
# ==========================================================
obj <- readRDS("airbox_stl_deseasonalized.rds")
X_adj <- obj$X_adj                 # T x K
datetime <- obj$datetime
sensor_table <- obj$sensor_table   # sensor_id, longitude, latitude, ...

# Standardize per sensor (recommended)
X_std <- scale(X_adj, center = TRUE, scale = TRUE)  # T x K

# Sensor ids from column names: expects "S123" or "123"
cn <- colnames(X_std)
sensor_id <- suppressWarnings(as.integer(sub("^S", "", cn)))

# ==========================================================
# 1) Helper functions: distances, indices, selection, plotting
# ==========================================================

# ---- ACF feature matrix (K x h_acf) ----
acf_features <- function(X_TK, h_acf = 5, show_progress = TRUE) {
  K <- ncol(X_TK)
  F <- matrix(NA_real_, nrow = K, ncol = h_acf)
  rownames(F) <- colnames(X_TK)
  colnames(F) <- paste0("lag", 1:h_acf)
  
  if (show_progress) {
    pb <- txtProgressBar(min = 0, max = K, style = 3)
    on.exit(close(pb), add = TRUE)
  }
  
  for (j in 1:K) {
    a <- acf(X_TK[, j], lag.max = h_acf, plot = FALSE)$acf
    F[j, ] <- as.numeric(a[2:(h_acf + 1)])
    if (show_progress) setTxtProgressBar(pb, j)
  }
  
  F
}

# ---- CC distance matrix (k x k) ----
cc_distance <- function(X_TK, h_cc = 5, show_progress = TRUE) {
  K <- ncol(X_TK)
  MCC <- matrix(0, K, K)
  
  # progress bar (one tick per i)
  if (show_progress) {
    pb <- txtProgressBar(min = 0, max = K, style = 3)
    on.exit(close(pb), add = TRUE)
  }
  
  for (i in 1:K) {
    # precompute ACF for i once
    aa <- acf(X_TK[, i], lag.max = h_cc, plot = FALSE)$acf
    aa2 <- c(aa, aa); aa2 <- aa2[-1]
    ma <- sqrt(sum(aa2 * aa2))
    
    for (j in 1:K) {
      ab <- ccf(X_TK[, i], X_TK[, j], lag.max = h_cc, plot = FALSE)$acf
      
      # precompute ACF for j (could be cached globally, but keeping simple)
      bb <- acf(X_TK[, j], lag.max = h_cc, plot = FALSE)$acf
      bb2 <- c(bb, bb); bb2 <- bb2[-1]
      mb <- sqrt(sum(bb2 * bb2))
      
      mab <- sum(ab * ab)
      MCC[i, j] <- 1 - mab / (ma * mb)
    }
    
    if (show_progress) setTxtProgressBar(pb, i)
  }
  
  MCC
}

# ---- Silhouette (average) given a distance matrix and membership ----
avg_silhouette <- function(D, memb) {
  ss <- cluster::silhouette(memb, as.dist(D))
  mean(ss[, "sil_width"])
}

# ---- CH index for k-means object ----
ch_index <- function(km, n) {
  K <- km$centers %>% nrow()
  (km$betweenss / (K - 1)) / (km$tot.withinss / (n - K))
}

# ---- Select K by silhouette (and optionally CH for kmeans) ----
select_K_kmeans <- function(X_km, Kmax = 10, criterion = c("silhouette", "CH"), seed = 1, nstart = 20) {
  criterion <- match.arg(criterion)
  n <- nrow(X_km)
  
  sil <- rep(NA_real_, Kmax)
  ch  <- rep(NA_real_, Kmax)
  
  D_euc <- as.matrix(dist(X_km))
  
  set.seed(seed)
  for (K in 2:Kmax) {
    km <- kmeans(X_km, centers = K, nstart = nstart)
    sil[K] <- avg_silhouette(D_euc, km$cluster)
    ch[K]  <- ch_index(km, n)
  }
  
  if (criterion == "silhouette") {
    K_opt <- which.max(sil)
  } else {
    K_opt <- which.max(ch)
  }
  
  set.seed(seed)
  km_final <- kmeans(X_km, centers = K_opt, nstart = max(nstart, 50))
  
  list(K_opt = K_opt, km = km_final, sil = sil, ch = ch, D = D_euc)
}

select_K_hclust <- function(D, Kmax = 10, hc_method = "complete") {
  hc <- hclust(as.dist(D), method = hc_method)
  sil <- rep(NA_real_, Kmax)
  
  for (K in 2:Kmax) {
    memb <- cutree(hc, k = K)
    sil[K] <- avg_silhouette(D, memb)
  }
  K_opt <- which.max(sil)
  memb_opt <- cutree(hc, k = K_opt)
  
  list(K_opt = K_opt, hc = hc, memb = memb_opt, sil = sil)
}

# ---- Taiwan map plot colored by membership ----
plot_taiwan_clusters <- function(sensor_table, memb, title = "Clusters", zoom = TRUE) {
  df <- sensor_table %>%
    mutate(cluster = factor(memb)) %>%
    filter(!is.na(longitude), !is.na(latitude), !is.na(cluster))
  
  world_map  <- map_data("world")
  taiwan_map <- world_map %>% filter(region == "Taiwan")
  
  p <- ggplot() +
    geom_polygon(data = taiwan_map, aes(long, lat, group = group),
                 fill = "grey95", color = "grey40", linewidth = 0.3) +
    geom_point(data = df, aes(longitude, latitude, color = cluster),
               size = 1.8, alpha = 0.85) +
    labs(title = title, x = "Longitude", y = "Latitude", color = "Cluster") +
    theme_minimal()
  
  if (zoom) {
    p <- p + coord_quickmap(xlim = c(119.5, 122), ylim = c(21.5, 25.5))
  } else {
    lon_range <- range(df$longitude, na.rm = TRUE)
    lat_range <- range(df$latitude,  na.rm = TRUE)
    lon_pad <- diff(lon_range) * 0.05
    lat_pad <- diff(lat_range) * 0.05
    p <- p + coord_quickmap(
      xlim = c(lon_range[1] - lon_pad, lon_range[2] + lon_pad),
      ylim = c(lat_range[1] - lat_pad, lat_range[2] + lat_pad)
    )
  }
  p
}

# ==========================================================
# 2) Build feature spaces and distance matrices needed
# ==========================================================
Kmax <- 10
h_acf <- 5
h_cc  <- 5

# ACF features (K x h)
F_acf <- acf_features(X_std, h_acf = h_acf)

# Time-series vectors for Euclidean k-means: each sensor is a row (K x T)
X_vec <- t(X_std)   # K x T

# Distances for HC methods:
# (3) HC on ACF feature Euclidean distance
D_acf_feat <- as.matrix(dist(F_acf))

# (4) HC on GCC distance (orientation: GCCmatrix typically expects T x K)
res_gcc <- GCCmatrix(X_std, 5)
D_gcc <- as.matrix(res_gcc$DM)

# (5) HC on CC distance
D_cc <- cc_distance(X_std, h_cc = h_cc)

# (6) HC on AR distance via diss METHOD="AR.PIC"
D_ar <- as.matrix(diss(t(X_std), METHOD = "AR.PIC"))

# (7) HC on ACF dissimilarity via diss METHOD="ACF"
# diss() expects series in rows: K x T
D_acf_diss <- as.matrix(diss(t(X_std), METHOD = "ACF", lag.max = h_acf))

# ==========================================================
# 3) METHOD 1: K-means on ACF coefficients
# ==========================================================
km1_sil <- select_K_kmeans(F_acf, Kmax = Kmax, criterion = "silhouette", seed = 1)
km1_ch  <- select_K_kmeans(F_acf, Kmax = Kmax, criterion = "CH", seed = 1)

memb_km_acf_sil <- km1_sil$km$cluster
memb_km_acf_ch  <- km1_ch$km$cluster

p_km_acf_sil <- plot_taiwan_clusters(sensor_table, memb_km_acf_sil,
                                     title = paste0("K-means on ACF coeffs | K* (sil) = ", km1_sil$K_opt))
p_km_acf_ch  <- plot_taiwan_clusters(sensor_table, memb_km_acf_ch,
                                     title = paste0("K-means on ACF coeffs | K* (CH) = ", km1_ch$K_opt))

# ==========================================================
# 4) METHOD 2: K-means on time-series vectors (Euclidean)
# ==========================================================
km2_sil <- select_K_kmeans(X_vec, Kmax = Kmax, criterion = "silhouette", seed = 1)
km2_ch  <- select_K_kmeans(X_vec, Kmax = Kmax, criterion = "CH", seed = 1)

memb_km_euc_sil <- km2_sil$km$cluster
memb_km_euc_ch  <- km2_ch$km$cluster

p_km_euc_sil <- plot_taiwan_clusters(sensor_table, memb_km_euc_sil,
                                     title = paste0("K-means on series (Eucl.) | K* (sil) = ", km2_sil$K_opt))
p_km_euc_ch  <- plot_taiwan_clusters(sensor_table, memb_km_euc_ch,
                                     title = paste0("K-means on series (Eucl.) | K* (CH) = ", km2_ch$K_opt))

# ==========================================================
# 5) METHOD 3: HC on ACF FEATURE distance
# ==========================================================
hc3 <- select_K_hclust(D_acf_feat, Kmax = Kmax, hc_method = "complete")
p_hc_acf_feat <- plot_taiwan_clusters(sensor_table, hc3$memb,
                                      title = paste0("HC on ACF features | K* (sil) = ", hc3$K_opt))

# ==========================================================
# 6) METHOD 4: HC on GCC distance
# ==========================================================
hc4 <- select_K_hclust(D_gcc, Kmax = Kmax, hc_method = "complete")
p_hc_gcc <- plot_taiwan_clusters(sensor_table, hc4$memb,
                                 title = paste0("HC on GCC distance | K* (sil) = ", hc4$K_opt))

# ==========================================================
# 7) METHOD 5: HC on CC distance
# ==========================================================
hc5 <- select_K_hclust(D_cc, Kmax = Kmax, hc_method = "complete")
p_hc_cc <- plot_taiwan_clusters(sensor_table, hc5$memb,
                                title = paste0("HC on CC distance | K* (sil) = ", hc5$K_opt))

# ==========================================================
# 8) METHOD 6: HC on AR distance
# ==========================================================
hc6 <- select_K_hclust(D_ar, Kmax = Kmax, hc_method = "complete")
p_hc_ar <- plot_taiwan_clusters(sensor_table, hc6$memb,
                                title = paste0("HC on AR.PIC distance | K* (sil) = ", hc6$K_opt))

# ==========================================================
# 9) METHOD 7: HC on ACF dissimilarity (TSclust diss METHOD="ACF")
# ==========================================================
hc7 <- select_K_hclust(D_acf_diss, Kmax = Kmax, hc_method = "complete")
p_hc_acf_diss <- plot_taiwan_clusters(sensor_table, hc7$memb,
                                      title = paste0("HC on ACF dissimilarity | K* (sil) = ", hc7$K_opt))

# ==========================================================
# 10) Summary table + plots
# ==========================================================
summary_tbl <- data.frame(
  method = c(
    "K-means ACF (sil)", "K-means ACF (CH)",
    "K-means series Eucl (sil)", "K-means series Eucl (CH)",
    "HC ACF features (sil)",
    "HC GCC (sil)",
    "HC CC (sil)",
    "HC AR.PIC (sil)",
    "HC ACF dissimilarity (sil)"
  ),
  K_opt = c(
    km1_sil$K_opt, km1_ch$K_opt,
    km2_sil$K_opt, km2_ch$K_opt,
    hc3$K_opt,
    hc4$K_opt,
    hc5$K_opt,
    hc6$K_opt,
    hc7$K_opt
  ),
  score = c(
    max(km1_sil$sil, na.rm = TRUE), max(km1_ch$ch, na.rm = TRUE),
    max(km2_sil$sil, na.rm = TRUE), max(km2_ch$ch, na.rm = TRUE),
    max(hc3$sil, na.rm = TRUE),
    max(hc4$sil, na.rm = TRUE),
    max(hc5$sil, na.rm = TRUE),
    max(hc6$sil, na.rm = TRUE),
    max(hc7$sil, na.rm = TRUE)
  )
)

print(summary_tbl)

print(p_km_acf_sil)
print(p_km_euc_sil)
print(p_hc_acf_feat)
print(p_hc_gcc)
print(p_hc_cc)
print(p_hc_ar)
print(p_hc_acf_diss)

saveRDS(
  list(
    km1_sil = km1_sil,
    km1_ch = km1_ch,
    km2_sil = km2_sil,
    km2_ch = km2_ch,
    hc3 = hc3,
    hc4 = hc4,
    hc5 = hc5,
    hc6 = hc6,
    hc7 = hc7
  ),
  file = "airbox_clusters.rds"
)
