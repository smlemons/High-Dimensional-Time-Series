############################################################
# Taiwan AirBox — Seasonality removal using STL (per sensor)
# Goal: produce a matrix suitable for VAR by removing the STL seasonal component
# Steps:
# 1) Load data (as before), remove V29 & V70
# 2) Plot ACF/PACF BEFORE (subset + average)
# 3) Apply STL per sensor with frequency = 24 (daily seasonality)
# 4) Subtract the estimated seasonal component (optionally also remove trend)
# 5) Plot ACF/PACF AFTER
# 6) Save the deseasonalized matrix for VAR
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(ggplot2)
  library(forecast)   # ggAcf, ggPacf, na.interp
  library(gridExtra)
  library(grid)
})

# -----------------------------
# 0) Theme
# -----------------------------
theme_project <- function(base_size = 12) {
  theme_bw() +
    theme(
      plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size, margin = margin(b = 6)),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = element_blank()
    )
}
theme_set(theme_project())

# -----------------------------
# 1) Load data (as before)
# -----------------------------
air_raw <- fread("TaiwanAirBox032017.csv")
loc <- fread("locations032017.csv")
setnames(loc, old = "V1", new = "sensor_id")

sensor_cols <- setdiff(names(air_raw), "time")

start_time <- ymd_hms("2017-03-01 00:00:00", tz = "Asia/Taipei")

air <- air_raw %>%
  mutate(
    time_index = as.integer(time),
    datetime   = start_time + hours(time_index - 1)
  ) %>%
  select(datetime, all_of(sensor_cols))

# Remove V29 & V70
remove_cols <- paste0("V", c(29, 70))
sensor_cols <- setdiff(sensor_cols, remove_cols)
air <- air %>% select(datetime, all_of(sensor_cols))

sensor_id_from_name <- as.integer(sub("^V", "", sensor_cols))

X <- as.matrix(air[, ..sensor_cols])
colnames(X) <- paste0("S", sensor_id_from_name)

cat("Loaded: T =", nrow(X), "hours; K =", ncol(X), "sensors\n")

# -----------------------------
# 2) ACF/PACF plotting helpers
# -----------------------------
plot_acf_pacf_subset <- function(Xmat, sensor_ids, lag_max = 120, title_tag = "") {
  ids_available <- as.integer(sub("^S", "", colnames(Xmat)))
  sensor_ids <- sensor_ids[sensor_ids %in% ids_available]
  
  acf_plots <- lapply(sensor_ids, function(sid) {
    j <- which(ids_available == sid)
    ggAcf(Xmat[, j], lag.max = lag_max) +
      labs(title = paste0("Sensor ", sid, " ACF ", title_tag)) +
      theme_project(9)
  })
  
  pacf_plots <- lapply(sensor_ids, function(sid) {
    j <- which(ids_available == sid)
    ggPacf(Xmat[, j], lag.max = lag_max) +
      labs(title = paste0("Sensor ", sid, " PACF ", title_tag)) +
      theme_project(9)
  })
  
  grid.arrange(grobs = acf_plots, ncol = 3,
               top = textGrob(paste0("ACF (subset) ", title_tag),
                              gp = gpar(fontsize = 14, fontface = "bold")))
  grid.arrange(grobs = pacf_plots, ncol = 3,
               top = textGrob(paste0("PACF (subset) ", title_tag),
                              gp = gpar(fontsize = 14, fontface = "bold")))
}

plot_average_acf_pacf <- function(Xmat, lag_max = 120, title_tag = "") {
  acf_mat <- sapply(seq_len(ncol(Xmat)), function(j) {
    stats::acf(Xmat[, j], lag.max = lag_max, plot = FALSE, na.action = na.pass)$acf[, 1, 1]
  })
  pacf_mat <- sapply(seq_len(ncol(Xmat)), function(j) {
    stats::pacf(Xmat[, j], lag.max = lag_max, plot = FALSE, na.action = na.pass)$acf
  })
  
  df_acf <- data.frame(
    lag = 0:lag_max,
    mean = rowMeans(acf_mat, na.rm = TRUE),
    sd   = apply(acf_mat, 1, sd, na.rm = TRUE)
  ) %>% mutate(lo = mean - sd, hi = mean + sd)
  
  df_pacf <- data.frame(
    lag = 1:lag_max,
    mean = rowMeans(pacf_mat, na.rm = TRUE),
    sd   = apply(pacf_mat, 1, sd, na.rm = TRUE)
  ) %>% mutate(lo = mean - sd, hi = mean + sd)
  
  p1 <- ggplot(df_acf, aes(lag, mean)) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_vline(xintercept = seq(24, lag_max, by = 24),
               linetype = "dashed", color = "grey60", linewidth = 0.6) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15) +
    geom_line(linewidth = 0.9) +
    labs(title = paste0("Average ACF across sensors ", title_tag),
         subtitle = "Dashed lines at 24 multiples",
         x = "Lag", y = "Mean ACF")
  
  p2 <- ggplot(df_pacf, aes(lag, mean)) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_vline(xintercept = seq(24, lag_max, by = 24),
               linetype = "dashed", color = "grey60", linewidth = 0.6) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15) +
    geom_line(linewidth = 0.9) +
    labs(title = paste0("Average PACF across sensors ", title_tag),
         subtitle = "Dashed lines at 24 multiples",
         x = "Lag", y = "Mean PACF")
  
  grid.arrange(p1, p2, ncol = 2)
}

# -----------------------------
# 3) BEFORE: ACF/PACF
# -----------------------------
subset_ids <- 1:6
lag_max <- 120
plot_acf_pacf_subset(X, subset_ids, lag_max = lag_max, title_tag = "(raw)")
plot_average_acf_pacf(X, lag_max = lag_max, title_tag = "(raw)")

# -----------------------------
# 4) STL deseasonalization (per sensor)
# -----------------------------
# STL decomposition assumes:
#   x_t = trend_t + seasonal_t + remainder_t
# For hourly data, we set frequency = 24 to capture DAILY seasonality.
# We remove seasonal_t by:
#   x_deseas_t = x_t - seasonal_t
#
# Important:
# - STL captures a potentially time-varying seasonal pattern (locally smoothed),
#   unlike simple hour-of-day averaging which forces a fixed pattern.
# - "s.window = 'periodic'" forces the seasonal component to be strictly periodic
#   (stable across time). This is often appropriate for only 31 days of data.
#
# If you want a seasonality that can slowly change over the month, set s.window
# to a number (e.g., s.window = 49 or 73). Larger -> smoother seasonal evolution.

stl_deseason <- function(x, freq = 24, s_window = "periodic", remove_trend = FALSE) {
  x <- as.numeric(x)
  if (any(!is.finite(x))) x <- forecast::na.interp(x)
  
  tsx <- ts(x, frequency = freq)
  fit <- stl(tsx, s.window = s_window, robust = TRUE)
  
  seas <- fit$time.series[, "seasonal"]
  trend <- fit$time.series[, "trend"]
  rem <- fit$time.series[, "remainder"]
  
  if (remove_trend) {
    # keep only remainder (removes both trend + seasonality)
    out <- rem
  } else {
    # remove seasonality only (keeps trend + remainder)
    out <- x - seas
  }
  
  list(
    y = as.numeric(out),
    seasonal = as.numeric(seas),
    trend = as.numeric(trend),
    remainder = as.numeric(rem)
  )
}

# Apply STL to all sensors
freq_season <- 24
s_window <- "periodic"
remove_trend <- FALSE   # set TRUE if you also want to remove slow trend

X_stl_adj <- matrix(NA_real_, nrow = nrow(X), ncol = ncol(X))
colnames(X_stl_adj) <- colnames(X)

# (Optional) store seasonal components for one sensor to illustrate
sid_demo <- subset_ids[subset_ids %in% as.integer(sub("^S", "", colnames(X)))][1]
demo_out <- NULL

for (j in seq_len(ncol(X))) {
  res <- stl_deseason(X[, j], freq = freq_season, s_window = s_window, remove_trend = remove_trend)
  X_stl_adj[, j] <- res$y
  
  # store decomposition for one demo sensor
  if (!is.na(sid_demo) && colnames(X)[j] == paste0("S", sid_demo)) {
    demo_out <- res
  }
  if (j %% 50 == 0) cat("STL processed", j, "/", ncol(X), "\n")
}

# Center each series (recommended before VAR with intercept)
X_stl_adj <- scale(X_stl_adj, center = TRUE, scale = FALSE)

# -----------------------------
# 5) AFTER: ACF/PACF
# -----------------------------
plot_acf_pacf_subset(X_stl_adj, subset_ids, lag_max = lag_max, title_tag = "(STL deseasonalized)")
plot_average_acf_pacf(X_stl_adj, lag_max = lag_max, title_tag = "(STL deseasonalized)")

# -----------------------------
# 6) Explain STL visually for one sensor (optional plot)
# -----------------------------
if (!is.null(demo_out)) {
  df_demo <- data.frame(
    datetime = air$datetime,
    raw = as.numeric(X[, which(colnames(X) == paste0("S", sid_demo))]),
    seasonal = demo_out$seasonal,
    trend = demo_out$trend,
    remainder = demo_out$remainder,
    deseasonalized = as.numeric(X_stl_adj[, which(colnames(X) == paste0("S", sid_demo))])
  ) %>%
    pivot_longer(-datetime, names_to = "component", values_to = "value")
  
  p_demo <- ggplot(df_demo, aes(x = datetime, y = value)) +
    geom_line(linewidth = 0.4) +
    facet_wrap(~ component, ncol = 1, scales = "free_y") +
    labs(
      title = paste0("STL decomposition components (Sensor ", sid_demo, ")"),
      subtitle = "raw = trend + seasonal + remainder; deseasonalized = raw - seasonal (then centered)",
      x = NULL, y = NULL
    )
  print(p_demo)
}

# -----------------------------
# 7) Save for VAR
# -----------------------------
saveRDS(
  list(
    X_adj = X_stl_adj,
    datetime = air$datetime,
    sensor_ids = as.integer(sub("^S", "", colnames(X_stl_adj))),
    stl = list(freq = freq_season, s_window = s_window, robust = TRUE, remove_trend = remove_trend),
    note = "X_adj is STL-deseasonalized (seasonal removed) and centered per sensor."
  ),
  file = "airbox_stl_deseasonalized_for_VAR.rds"
)

cat("\nSaved: airbox_stl_deseasonalized_for_VAR.rds\n")
