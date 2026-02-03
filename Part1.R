############################################################
# Part 1 — Taiwan AirBox EDA + stationarity/seasonality checks
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(forecast)   # Acf, Pacf, nsdiffs
  library(corrplot)
  library(maps)
  library(gridExtra)
  library(urca)       # ur.df (ADF)
  library(grid)
})

# -----------------------------
# 0) Global ggplot theme
# -----------------------------
theme_project <- function(base_size = 13) {
  theme_bw() +
    theme(
      plot.title = element_text(
        size = base_size + 3,
        face = "bold",
        hjust = 0
      ),
      plot.subtitle = element_text(
        size = base_size,
        margin = margin(b = 8)
      ),
      axis.title = element_text(
        size = base_size + 1
      ),
      axis.text = element_text(
        size = base_size - 1
      ),
      legend.title = element_text(
        size = base_size
      ),
      legend.text = element_text(
        size = base_size - 1
      ),
      panel.grid.major = element_line(
        color = "grey85",
        linewidth = 0.3
      ),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# Set as global default
theme_set(theme_project())


# -----------------------------
# 1) Load the data
# -----------------------------
air <- fread("TaiwanAirBox032017.csv")      # columns: time, V1..V516
loc <- fread("locations032017.csv")        # columns: V1, latitude, longitude
setnames(loc, old = "V1", new = "sensor_id")

sensor_cols <- setdiff(names(air), "time")
sensor_id_from_name <- as.integer(sub("^V", "", sensor_cols))

# -----------------------------
# 2) Build datetime index
# -----------------------------
start_time <- ymd_hms("2017-03-01 00:00:00", tz = "Asia/Taipei")

air <- air %>%
  mutate(
    time_index = as.integer(time),
    datetime   = start_time + hours(time_index - 1)
  ) %>%
  select(datetime, all_of(sensor_cols))

# -----------------------------
# 3) Data matrix (T x K)
# -----------------------------
X <- as.matrix(air[, ..sensor_cols])
colnames(X) <- paste0("S", sensor_id_from_name)

Tn <- nrow(X)
K  <- ncol(X)
cat("T =", Tn, "hours; K =", K, "sensors\n")

# -----------------------------
# 4) Locations table
# -----------------------------
sensor_table <- data.frame(
  sensor_id   = sensor_id_from_name,
  sensor_name = paste0("S", sensor_id_from_name)
) %>%
  left_join(loc, by = "sensor_id")

# -----------------------------
# 5) Descriptive summaries (raw scale) + average correlation
# -----------------------------
raw_mean <- colMeans(X)
raw_sd   <- apply(X, 2, sd)

desc <- sensor_table %>%
  mutate(
    mean_raw = raw_mean,
    sd_raw   = raw_sd
  )

print(desc %>% summarise(
  mean_of_means = mean(mean_raw),
  median_sd     = median(sd_raw),
  min_sd        = min(sd_raw),
  max_sd        = max(sd_raw)
))

# Correlation matrix on standardized series
X_scaled <- scale(X, center = TRUE, scale = TRUE)
C0_full <- cor(X_scaled)

# Global average correlation (off-diagonal)
off_diag <- C0_full[upper.tri(C0_full)]
avg_corr <- mean(off_diag)
cat("Average lag-0 correlation (off-diagonal):", avg_corr, "\n")

# Per-sensor average correlation with all other sensors (signed, not absolute)
avg_corr_per_sensor <- (rowSums(C0_full) - 1) / (K - 1)

desc <- desc %>%
  mutate(avg_corr = avg_corr_per_sensor)

# Identify the 3 sensors with the lowest mean
outlier_ids <- desc %>%
  arrange(mean_raw) %>%
  slice(1:3) %>%
  pull(sensor_id)

# Add indicator for the three identified series
desc_plot <- desc %>%
  mutate(is_outlier = sensor_id %in% outlier_ids)

# Scatter plot with different shape for the three series
ggplot(desc_plot, aes(x = mean_raw, y = sd_raw)) +
  geom_point(
    aes(shape = is_outlier),
    size = 3,
    alpha = 0.6
  )  +
  geom_text(
    data = desc %>% filter(sensor_id %in% outlier_ids),
    aes(label = sensor_id),
    hjust = -0.4,
    vjust = 0,
    size = 5
  ) +
  scale_shape_manual(
    values = c(`FALSE` = 19, `TRUE` = 24),   # circle vs triangle
    labels = c("Normal series", "Low-mean series")
  ) +
  labs(
    title = "Mean vs standard deviation across series",
    x = "Mean",
    y = "Standard deviation",
    shape = NULL
  )+ theme(
    legend.position = "inside",
    legend.position.inside = c(0.87, 0.10)
  )
# -----------------------------
# 7) Plot of the three selected outlier sensors (ordered)
# -----------------------------

outlier_cols <- paste0("V", outlier_ids)

outlier_df <- air %>%
  select(datetime, all_of(outlier_cols)) %>%
  pivot_longer(-datetime, names_to = "sensor", values_to = "value") %>%
  mutate(
    sensor_id = as.integer(sub("^V", "", sensor)),
    sensor_label = paste("Series", sensor_id)
  )

# Enforce numeric ordering of facets
sensor_levels <- paste("Series", sort(unique(outlier_df$sensor_id)))

outlier_df <- outlier_df %>%
  mutate(
    sensor_label = factor(sensor_label, levels = sensor_levels)
  )

ggplot(outlier_df, aes(x = datetime, y = value)) +
  geom_line() +
  facet_wrap(
    ~ sensor_label,
    scales = "free_y",
    ncol = 1
  ) +
  labs(
    title = "Outlier time series",
    x = NULL,
    y = "Value"
  )


# -----------------------------
# 7) Plot of all the time series + average series
#     Full period + date-defined zoom
# -----------------------------

# Average series
avg_series <- rowMeans(X)

# Long format once (full period)
air_long_full <- air %>%
  select(datetime, all_of(sensor_cols)) %>%
  pivot_longer(-datetime, names_to = "sensor", values_to = "value")

df_avg_full <- data.frame(
  datetime = air$datetime,
  avg = avg_series
)

p_full <- ggplot() +
  geom_line(
    data = air_long_full,
    aes(x = datetime, y = value, group = sensor),
    alpha = 0.03
  ) +
  geom_line(
    data = df_avg_full,
    aes(x = datetime, y = avg, color = "Average"),
    linewidth = 0.9
  ) +
  scale_color_manual(values = c("Average" = "red")) +
  labs(
    title = "All sensor time series (31 days)",
    x = NULL,
    y = "Value",
    color = NULL
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.05, 0.95),
    legend.justification = c(0, 1),
    legend.background = element_rect(
      fill = "white",
      color = "grey70",
      linewidth = 0.3
    ),
    legend.key = element_blank()
  )

# Date-defined zoom window
zoom_start <- ymd_hms("2017-03-20 00:00:00", tz = "Asia/Taipei")
zoom_end   <- ymd_hms("2017-03-24 00:00:00", tz = "Asia/Taipei")

air_zoom <- air %>%
  filter(datetime >= zoom_start, datetime <= zoom_end)

air_long_zoom <- air_zoom %>%
  select(datetime, all_of(sensor_cols)) %>%
  pivot_longer(-datetime, names_to = "sensor", values_to = "value")

df_avg_zoom <- data.frame(
  datetime = air_zoom$datetime,
  avg = rowMeans(as.matrix(air_zoom[, ..sensor_cols]))
)

# Midnight breaks for clean x-axis
breaks_midnight <- seq(
  from = floor_date(zoom_start, unit = "day"),
  to   = ceiling_date(zoom_end, unit = "day"),
  by   = "1 day"
)

p_zoom <- ggplot() +
  geom_line(
    data = air_long_zoom,
    aes(x = datetime, y = value, group = sensor),
    alpha = 0.05
  ) +
  geom_line(
    data = df_avg_zoom,
    aes(x = datetime, y = avg, color = "Average"),
    linewidth = 0.9
  ) +
  scale_color_manual(values = c("Average" = "red")) +
  scale_x_datetime(
    breaks = breaks_midnight,
    date_labels = "%d %b\n%H:%M"
  ) +
  labs(
    title = "All sensor time series (zoomed window)",
    x = NULL,
    y = "Value",
    color = NULL
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.05, 0.95),
    legend.justification = c(0, 1),
    legend.background = element_rect(
      fill = "white",
      color = "grey70",
      linewidth = 0.3
    ),
    legend.key = element_blank()
  )

grid.arrange(p_full, p_zoom, ncol = 2)

# -----------------------------
# 8) Full correlation heatmap (lag 0, all sensors)
# -----------------------------
# Labeling 500+ series is not readable; this renders the full matrix without labels.
corrplot(
  C0_full,
  method = "color",
  tl.pos = "n",
  order = "hclust",
  mar = c(0, 0, 2, 0),
  title = "Lag-0 correlation matrix (all sensors, clustered)"
)

# -----------------------------
# 9) Plot of sensors over the world and Taiwan
# -----------------------------
world_map  <- map_data("world")
taiwan_map <- world_map %>% filter(region == "Taiwan")

lon_range <- range(desc$longitude, na.rm = TRUE)
lat_range <- range(desc$latitude,  na.rm = TRUE)
lon_pad <- diff(lon_range) * 0.05
lat_pad <- diff(lat_range) * 0.05
lon_lim <- c(lon_range[1] - lon_pad, lon_range[2] + lon_pad)
lat_lim <- c(lat_range[1] - lat_pad, lat_range[2] + lat_pad)

p_global <- ggplot() +
  geom_polygon(
    data = world_map,
    aes(x = long, y = lat, group = group),
    fill = "grey95",
    color = "grey60",
    linewidth = 0.2
  ) +
  geom_point(
    data = desc,
    aes(x = longitude, y = latitude, color = mean_raw),
    size = 1.5,
    alpha = 0.6
  ) +
  coord_quickmap(xlim = lon_lim, ylim = lat_lim) +
  scale_color_viridis_c(option = "C", name = "Mean value") +
  labs(title = "Series coordinates (global)", x = "Longitude", y = "Latitude") +
  theme_minimal()

p_taiwan <- ggplot() +
  geom_polygon(
    data = taiwan_map,
    aes(x = long, y = lat, group = group),
    fill = "grey95",
    color = "grey40",
    linewidth = 0.3
  ) +
  geom_point(
    data = desc,
    aes(x = longitude, y = latitude, color = mean_raw),
    size = 1.5,
    alpha = 0.6
  ) +
  coord_quickmap(xlim = c(119.5, 122), ylim = c(21.5, 25.5)) +
  scale_color_viridis_c(option = "C", name = "Mean value") +
  labs(title = "Series coordinates (zoomed)", x = "Longitude", y = "Latitude") +
  theme_minimal()

grid.arrange(p_global, p_taiwan, ncol = 2)

# -----------------------------
# 10) ACF and PACF analysis (subset of sensors)
# Clean 2x3 layout
# -----------------------------

acf_pacf_ids <- 1:6
lag_max <- 80

# ACF plots
acf_plots <- lapply(acf_pacf_ids, function(sid) {
  j <- which(sensor_id_from_name == sid)
  ggAcf(
    X[, j],
    lag.max = lag_max
  ) +
    labs(title = paste("Sensor", sid)) +
    theme_project(8)
})

# PACF plots
pacf_plots <- lapply(acf_pacf_ids, function(sid) {
  j <- which(sensor_id_from_name == sid)
  ggPacf(
    X[, j],
    lag.max = lag_max
  ) +
    labs(title = paste("Sensor", sid)) +
    theme_project(8)
})

# Arrange in 2x3 grids
grid.arrange(
  grobs = acf_plots,
  ncol = 3,
  top = textGrob(
    "Autocorrelation functions (ACF)",
    gp = gpar(fontsize = 16, fontface = "bold")
  )
)

grid.arrange(
  grobs = pacf_plots,
  ncol = 3,
  top = textGrob(
    "Partial autocorrelation functions (PACF)",
    gp = gpar(fontsize = 16, fontface = "bold")
  )
)

# Compute ACF matrix: rows = lags (0..lag_max), cols = sensors
acf_mat <- sapply(seq_len(ncol(X)), function(j) {
  stats::acf(X[, j], lag.max = lag_max, plot = FALSE)$acf[, 1, 1]
})

# Compute PACF matrix: rows = lags (1..lag_max), cols = sensors
pacf_mat <- sapply(seq_len(ncol(X)), function(j) {
  stats::pacf(X[, j], lag.max = lag_max, plot = FALSE)$acf
})

# Mean and variability (sd across sensors) at each lag
acf_lags <- 0:lag_max
acf_mean <- rowMeans(acf_mat, na.rm = TRUE)
acf_sd   <- apply(acf_mat, 1, sd, na.rm = TRUE)

pacf_lags <- 1:lag_max
pacf_mean <- rowMeans(pacf_mat, na.rm = TRUE)
pacf_sd   <- apply(pacf_mat, 1, sd, na.rm = TRUE)

# Data frames for plotting
df_acf <- data.frame(
  lag = acf_lags,
  mean = acf_mean,
  lo = acf_mean - acf_sd,
  hi = acf_mean + acf_sd
)

df_pacf <- data.frame(
  lag = pacf_lags,
  mean = pacf_mean,
  lo = pacf_mean - pacf_sd,
  hi = pacf_mean + pacf_sd
)

# Plots (mean curve + +/- 1 sd band across sensors)
p_acf_avg <- ggplot(df_acf, aes(x = lag, y = mean)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_vline(
    xintercept = seq(24, max(df_acf$lag), by = 24),
    linetype = "dashed",
    color = "grey60",
    linewidth = 1
  ) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15) +
  geom_line(linewidth = 0.9) +
  labs(
    title = "Average ACF across all sensors",
    x = "Lag",
    y = "Mean ACF"
  )

p_acf_avg

p_pacf_avg <- ggplot(df_pacf, aes(x = lag, y = mean)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_vline(
    xintercept = seq(24, max(df_pacf$lag), by = 24),
    linetype = "dashed",
    color = "grey60",
    linewidth = 1
  ) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15) +
  geom_line(linewidth = 0.9) +
  labs(
    title = "Average PACF across all sensors",
    x = "Lag",
    y = "Mean PACF"
  )

p_pacf_avg


# -----------------------------
# 11) ADF and seasonality check (subset of sensors)
# -----------------------------
adf_type <- "drift"     # "none" / "drift" / "trend"
adf_lags <- 12          # ur.df selects lags by AIC within this cap
season_period <- 24

adf_results <- lapply(acf_pacf_ids, function(sid) {
  j <- which(sensor_id_from_name == sid)
  x <- X[, j]
  
  adf_fit <- ur.df(y = x, type = adf_type, lags = adf_lags, selectlags = "AIC")
  seas_diffs <- nsdiffs(x, m = season_period)
  
  list(
    sensor_id = sid,
    adf_summary = summary(adf_fit),
    nsdiffs = seas_diffs
  )
})

for (r in adf_results) {
  cat("\nADF (ur.df) for sensor", r$sensor_id, "\n")
  print(r$adf_summary)
  cat("nsdiffs (seasonal differences suggested, m=24):", r$nsdiffs, "\n")
}

# -----------------------------
# 12) Apply transformations and recheck ACF/PACF (subset)
# -----------------------------
transform_type <- "none"  # "none" / "diff" / "seas24" / "diff+seas24"

apply_transform <- function(x, type = c("none","diff","seas24","diff+seas24")) {
  type <- match.arg(type)
  if (type == "none") return(x)
  if (type == "diff") return(diff(x, lag = 1))
  if (type == "seas24") return(diff(x, lag = 24))
  y <- diff(x, lag = 24)
  y <- diff(y, lag = 1)
  y
}

par(mfrow = c(3, 4), mar = c(3, 3, 2, 1))
for (sid in acf_pacf_ids) {
  j <- which(sensor_id_from_name == sid)
  y <- apply_transform(X[, j], transform_type)
  Acf(y, lag.max = 50, main = paste("ACF", sid, transform_type))
}
par(mfrow = c(3, 4), mar = c(3, 3, 2, 1))
for (sid in acf_pacf_ids) {
  j <- which(sensor_id_from_name == sid)
  y <- apply_transform(X[, j], transform_type)
  Pacf(y, lag.max = 50, main = paste("PACF", sid, transform_type))
}
par(mfrow = c(1, 1))

# -----------------------------
# 13) AR order selection for each sensor (fast method)
# -----------------------------
p_max_ar <- 48
ar_method <- "burg"   # "burg" or "yw"

p_hat <- integer(K)

for (i in seq_len(K)) {
  x <- apply_transform(X[, i], transform_type)
  
  if (length(x) <= (p_max_ar + 10)) {
    p_hat[i] <- NA_integer_
    next
  }
  
  fit_aic <- ar(x, order.max = p_max_ar, aic = TRUE, method = ar_method)
  
  k_vec <- 0:p_max_ar
  n_eff <- length(x)
  bic_vec <- fit_aic$aic + (log(n_eff) - 2) * k_vec
  
  p_hat[i] <- which.min(bic_vec) - 1
  
  if (i %% 50 == 0) cat("AR order selection:", i, "/", K, "\n")
}

ar_dist <- data.frame(sensor = colnames(X), p_hat = p_hat)

print(summary(ar_dist$p_hat))

hist(
  ar_dist$p_hat,
  breaks = 0:p_max_ar,
  main = paste("Selected AR order distribution (BIC approx) — transform:", transform_type),
  xlab = "Selected p",
  col = "grey80",
  border = "white"
)

# -----------------------------
# Save objects
# -----------------------------
saveRDS(list(
  air = air,
  loc = loc,
  X = X,
  X_scaled = X_scaled,
  sensor_table = sensor_table,
  desc = desc,
  C0_full = C0_full,
  ar_dist = ar_dist,
  transform_type = transform_type
), file = "part1_objects.rds")

cat("Part 1 complete: saved to part1_objects.rds\n")
