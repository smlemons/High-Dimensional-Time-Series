############################################################
# Part 1 — Taiwan AirBox EDA + STL deseasonalization (VAR input)
#
# Changes vs previous Part1:
# - DO NOT use 24-hour seasonal differencing anymore
# - Replace it with STL deseasonalization (frequency = 24), per sensor
# - Keep: key EDA + ACF/PACF before/after + AR order selection (local)
# - Keep: plot of the 24-hour seasonal pattern (hour-of-day)
# - Skip: trend plot + average seasonal component over time
# - Save ONLY the cleaned + deseasonalized matrix + datetime + sensor_table
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(forecast)   # ggAcf, ggPacf
  library(corrplot)
  library(maps)
  library(gridExtra)
  library(grid)
  library(tseries)    # adf.test
})

# -----------------------------
# 0) Theme + plot saving
# -----------------------------
theme_set(
  theme_bw(base_size = 13) +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 13, margin = margin(b = 6)),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 12),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    )
)

plot_dir <- "plots_part1_stl"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

save_plot <- function(p, name, width = 10, height = 6) {
  ggsave(file.path(plot_dir, paste0(name, ".png")), p, width = width, height = height, dpi = 300)
}

save_grid_png <- function(expr, name, width = 14, height = 8) {
  png(file.path(plot_dir, paste0(name, ".png")), width = width, height = height, units = "in", res = 300)
  force(expr)
  dev.off()
}

# -----------------------------
# 1) Load data + datetime
# -----------------------------
air_raw <- fread("TaiwanAirBox032017.csv")  # time, V1..V516
loc     <- fread("locations032017.csv")     # V1, latitude, longitude
setnames(loc, old = "V1", new = "sensor_id")

start_time <- ymd_hms("2017-03-01 00:00:00", tz = "Asia/Taipei")

sensor_cols <- setdiff(names(air_raw), "time")

air <- air_raw %>%
  mutate(
    time_index = as.integer(time),
    datetime   = start_time + hours(time_index - 1)
  ) %>%
  select(datetime, all_of(sensor_cols))

# -----------------------------
# 2) Basic matrix + metadata
# -----------------------------
sensor_id_from_name <- as.integer(sub("^V", "", sensor_cols))

X <- as.matrix(air[, ..sensor_cols])
colnames(X) <- paste0("S", sensor_id_from_name)

cat("Initial: T =", nrow(X), "hours; K =", ncol(X), "sensors\n")

sensor_table <- data.frame(
  sensor_id   = sensor_id_from_name,
  sensor_name = paste0("S", sensor_id_from_name)
) %>%
  left_join(loc, by = "sensor_id")

# -----------------------------
# 3) Descriptives + outlier inspection (remove V29 & V70)
# -----------------------------
raw_mean <- colMeans(X, na.rm = TRUE)
raw_sd   <- apply(X, 2, sd, na.rm = TRUE)

desc <- sensor_table %>% mutate(mean_raw = raw_mean, sd_raw = raw_sd)

# Plot mean against standard deviation of the series
outlier_ids <- desc %>% arrange(mean_raw) %>% slice(1:3) %>% pull(sensor_id)
desc_plot <- desc %>% mutate(is_outlier = sensor_id %in% outlier_ids)

p_outlier_scatter <- ggplot(desc_plot, aes(mean_raw, sd_raw)) +
  geom_point(aes(shape = is_outlier), size = 3, alpha = 0.6) +
  geom_text(
    data = desc %>% filter(sensor_id %in% outlier_ids),
    aes(label = sensor_id),
    hjust = -0.4, vjust = 0, size = 4
  ) +
  scale_shape_manual(values = c(`FALSE` = 19, `TRUE` = 24)) +
  labs(title = "Mean vs SD across sensors", x = "Mean", y = "SD", shape = NULL)

save_plot(p_outlier_scatter, "outlier_mean_vs_sd", width = 9, height = 6)
print(p_outlier_scatter)

# Plot the outlier series identified above (for inspection)
outlier_cols <- paste0("V", outlier_ids)

outlier_df <- air %>%
  select(datetime, all_of(outlier_cols)) %>%
  pivot_longer(-datetime, names_to = "sensor", values_to = "value") %>%
  mutate(
    sensor_id = as.integer(sub("^V", "", sensor)),
    sensor_label = factor(paste("Series", sensor_id),
                          levels = paste("Series", sort(unique(sensor_id))))
  )

p_outlier_ts <- ggplot(outlier_df, aes(x = datetime, y = value)) +
  geom_line() +
  facet_wrap(~ sensor_label, scales = "free_y", ncol = 1) +
  labs(title = "Outlier time series", x = NULL, y = "Value")

save_plot(p_outlier_ts, "outlier_time_series", width = 9, height = 8)
print(p_outlier_ts)

# Remove V29 and V70 (fixed decision)
remove_cols <- paste0("V", c(29, 70))
sensor_cols <- setdiff(sensor_cols, remove_cols)

air <- air %>% select(datetime, all_of(sensor_cols))
sensor_id_from_name <- as.integer(sub("^V", "", sensor_cols))

X <- as.matrix(air[, ..sensor_cols])
colnames(X) <- paste0("S", sensor_id_from_name)

cat("After removing V29 & V70: T =", nrow(X), "hours; K =", ncol(X), "sensors\n")

sensor_table <- data.frame(
  sensor_id   = sensor_id_from_name,
  sensor_name = paste0("S", sensor_id_from_name)
) %>%
  left_join(loc, by = "sensor_id")

# Update desc (used for maps)
raw_mean <- colMeans(X, na.rm = TRUE)
raw_sd   <- apply(X, 2, sd, na.rm = TRUE)
desc <- sensor_table %>% mutate(mean_raw = raw_mean, sd_raw = raw_sd)

# -----------------------------
# 4) Overlay time series (full + zoom)
# -----------------------------
avg_series <- rowMeans(X, na.rm = TRUE)

air_long_full <- air %>%
  select(datetime, all_of(sensor_cols)) %>%
  pivot_longer(-datetime, names_to = "sensor", values_to = "value")

df_avg_full <- data.frame(datetime = air$datetime, avg = avg_series)

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
  labs(title = "All sensor time series (31 days)", x = NULL, y = "Value", color = NULL) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.05, 0.95),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", color = "grey70", linewidth = 0.3),
    legend.key = element_blank()
  )

zoom_start <- ymd_hms("2017-03-20 00:00:00", tz = "Asia/Taipei")
zoom_end   <- ymd_hms("2017-03-24 00:00:00", tz = "Asia/Taipei")

air_zoom <- air %>% filter(datetime >= zoom_start, datetime <= zoom_end)

air_long_zoom <- air_zoom %>%
  select(datetime, all_of(sensor_cols)) %>%
  pivot_longer(-datetime, names_to = "sensor", values_to = "value")

df_avg_zoom <- data.frame(
  datetime = air_zoom$datetime,
  avg = rowMeans(as.matrix(air_zoom[, ..sensor_cols]), na.rm = TRUE)
)

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
  scale_x_datetime(breaks = breaks_midnight, date_labels = "%d %b\n%H:%M") +
  labs(title = "All sensor time series (zoomed window)", x = NULL, y = "Value", color = NULL) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.05, 0.95),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = "white", color = "grey70", linewidth = 0.3),
    legend.key = element_blank()
  )

save_grid_png(
  grid.arrange(p_full, p_zoom, ncol = 2),
  "timeseries_full_and_zoom",
  width = 14, height = 6
)

grid.arrange(p_full, p_zoom, ncol = 2)

# -----------------------------
# 5) Correlation heatmap + maps
# -----------------------------
X_scaled <- scale(X, center = TRUE, scale = TRUE)
C0_full  <- cor(X_scaled, use = "pairwise.complete.obs")

png(file.path(plot_dir, "correlation_matrix_lag0_clustered.png"), width = 2000, height = 2000, res = 300)
corrplot(C0_full, method = "color", tl.pos = "n", order = "hclust",
         mar = c(0, 0, 2, 0), title = "Lag-0 correlation matrix (clustered)")
dev.off()

world_map  <- map_data("world")
taiwan_map <- world_map %>% filter(region == "Taiwan")

lon_range <- range(desc$longitude, na.rm = TRUE)
lat_range <- range(desc$latitude,  na.rm = TRUE)
lon_pad <- diff(lon_range) * 0.05
lat_pad <- diff(lat_range) * 0.05
lon_lim <- c(lon_range[1] - lon_pad, lon_range[2] + lon_pad)
lat_lim <- c(lat_range[1] - lat_pad, lat_range[2] + lat_pad)

p_global <- ggplot() +
  geom_polygon(data = world_map, aes(long, lat, group = group),
               fill = "grey95", color = "grey60", linewidth = 0.2) +
  geom_point(data = desc, aes(longitude, latitude, color = mean_raw),
             size = 1.5, alpha = 0.6) +
  coord_quickmap(xlim = lon_lim, ylim = lat_lim) +
  scale_color_viridis_c(option = "C", name = "Mean") +
  labs(title = "Sensor coordinates (global)", x = "Longitude", y = "Latitude") +
  theme_minimal()

p_taiwan <- ggplot() +
  geom_polygon(data = taiwan_map, aes(long, lat, group = group),
               fill = "grey95", color = "grey40", linewidth = 0.3) +
  geom_point(data = desc, aes(longitude, latitude, color = mean_raw),
             size = 1.5, alpha = 0.6) +
  coord_quickmap(xlim = c(119.5, 122), ylim = c(21.5, 25.5)) +
  scale_color_viridis_c(option = "C", name = "Mean") +
  labs(title = "Sensor coordinates (Taiwan zoom)", x = "Longitude", y = "Latitude") +
  theme_minimal()

save_grid_png(grid.arrange(p_global, p_taiwan, ncol = 2), "sensor_maps_global_and_taiwan", width = 12, height = 6)
grid.arrange(p_global, p_taiwan, ncol = 2)

# -----------------------------
# 6) ACF/PACF BEFORE STL (subset + averages)
# -----------------------------
acf_pacf_ids <- 1:6
lag_max <- 100

acf_plots_raw <- lapply(acf_pacf_ids, function(sid) {
  j <- which(sensor_id_from_name == sid)
  ggAcf(X[, j], lag.max = lag_max) + labs(title = paste("Sensor", sid))
})
pacf_plots_raw <- lapply(acf_pacf_ids, function(sid) {
  j <- which(sensor_id_from_name == sid)
  ggPacf(X[, j], lag.max = lag_max) + labs(title = paste("Sensor", sid))
})

save_grid_png(
  grid.arrange(grobs = acf_plots_raw, ncol = 3,
               top = textGrob("Autocorrelation functions (ACF) for raw series", gp = gpar(fontsize = 16, fontface = "bold"))),
  "acf_subset_raw",
  width = 14, height = 8
)
save_grid_png(
  grid.arrange(grobs = pacf_plots_raw, ncol = 3,
               top = textGrob("Partial autocorrelation functions (PACF) for raw series", gp = gpar(fontsize = 16, fontface = "bold"))),
  "pacf_subset_raw",
  width = 14, height = 8
)

grid.arrange(
  grobs = acf_plots_raw,
  ncol = 3,
  top = textGrob("Autocorrelation functions (ACF) for raw series", gp = gpar(fontsize=16, fontface = "bold"))
)

grid.arrange(
  grobs = acf_plots_raw,
  ncol = 3,
  top = textGrob("Partial autocorrelation functions (PACF) for raw series", gp = gpar(fontsize=16, fontface = "bold"))
)

acf_mat_raw <- sapply(seq_len(ncol(X)), function(j) stats::acf(X[, j], lag.max = lag_max, plot = FALSE)$acf[,1,1])
pacf_mat_raw <- sapply(seq_len(ncol(X)), function(j) stats::pacf(X[, j], lag.max = lag_max, plot = FALSE)$acf)

df_acf_raw  <- data.frame(lag = 0:lag_max, mean = rowMeans(acf_mat_raw),  sd = apply(acf_mat_raw, 1, sd))
df_pacf_raw <- data.frame(lag = 1:lag_max, mean = rowMeans(pacf_mat_raw), sd = apply(pacf_mat_raw, 1, sd))

p_avg_acf_raw <- ggplot(df_acf_raw, aes(lag, mean)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_vline(xintercept = seq(24, lag_max, by = 24), linetype = "dashed", color = "grey60") +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.15) +
  geom_line(linewidth = 0.9) +
  labs(title = "Average ACF (raw series)", x = "Lag", y = "Mean ACF")

p_avg_pacf_raw <- ggplot(df_pacf_raw, aes(lag, mean)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_vline(xintercept = seq(24, lag_max, by = 24), linetype = "dashed", color = "grey60") +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.15) +
  geom_line(linewidth = 0.9) +
  labs(title = "Average PACF (raw series)", x = "Lag", y = "Mean PACF")

save_grid_png(grid.arrange(p_avg_acf_raw, p_avg_pacf_raw, ncol = 2), "avg_acf_pacf_raw", width = 14, height = 6)
grid.arrange(p_avg_acf_raw, p_avg_pacf_raw, ncol = 2)

# -----------------------------
# 7) STL deseasonalization (frequency = 24), per sensor
# -----------------------------
freq <- 24
Tn <- nrow(X)
K  <- ncol(X)

season_stl <- matrix(NA_real_, Tn, K)
for (j in seq_len(K)) {
  fit <- stl(ts(X[, j], frequency = freq), s.window = "periodic", robust = TRUE)
  season_stl[, j] <- fit$time.series[, "seasonal"]
  if (j %% 50 == 0) cat("STL processed", j, "/", K, "\n")
}

# Deseasonalized VAR input: remove seasonal component only
X_adj <- X - season_stl
colnames(X_adj) <- colnames(X)

# -----------------------------
# 8) ACF/PACF AFTER STL (subset + averages)
# -----------------------------
acf_plots_adj <- lapply(acf_pacf_ids, function(sid) {
  j <- which(sensor_id_from_name == sid)
  ggAcf(X_adj[, j], lag.max = lag_max) + labs(title = paste("Sensor", sid))
})
pacf_plots_adj <- lapply(acf_pacf_ids, function(sid) {
  j <- which(sensor_id_from_name == sid)
  ggPacf(X_adj[, j], lag.max = lag_max) + labs(title = paste("Sensor", sid))
})

save_grid_png(
  grid.arrange(grobs = acf_plots_adj, ncol = 3,
               top = textGrob("Autocorrelation functions (ACF) for deseasonalized series", gp = gpar(fontsize = 16, fontface = "bold"))),
  "acf_subset_stl",
  width = 14, height = 8
)
save_grid_png(
  grid.arrange(grobs = pacf_plots_adj, ncol = 3,
               top = textGrob("Partial autocorrelation functions (PACF) for deseasonalized series", gp = gpar(fontsize = 16, fontface = "bold"))),
  "pacf_subset_stl",
  width = 14, height = 8
)

grid.arrange(grobs = acf_plots_adj, ncol = 3,
             top = textGrob("Autocorrelation functions (ACF) for deseasonalized series", gp = gpar(fontsize = 16, fontface = "bold")))
grid.arrange(grobs = pacf_plots_adj, ncol = 3,
             top = textGrob("Partial autocorrelation functions (PACF) for deseasonalized series", gp = gpar(fontsize = 16, fontface = "bold")))


acf_mat_adj <- sapply(seq_len(ncol(X_adj)), function(j) stats::acf(X_adj[, j], lag.max = lag_max, plot = FALSE)$acf[,1,1])
pacf_mat_adj <- sapply(seq_len(ncol(X_adj)), function(j) stats::pacf(X_adj[, j], lag.max = lag_max, plot = FALSE)$acf)

df_acf_adj  <- data.frame(lag = 0:lag_max, mean = rowMeans(acf_mat_adj),  sd = apply(acf_mat_adj, 1, sd))
df_pacf_adj <- data.frame(lag = 1:lag_max, mean = rowMeans(pacf_mat_adj), sd = apply(pacf_mat_adj, 1, sd))

p_avg_acf_adj <- ggplot(df_acf_adj, aes(lag, mean)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_vline(xintercept = seq(24, lag_max, by = 24), linetype = "dashed", color = "grey60") +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.15) +
  geom_line(linewidth = 0.9) +
  labs(title = "Average ACF (deseasonalized series)", x = "Lag", y = "Mean ACF")

p_avg_pacf_adj <- ggplot(df_pacf_adj, aes(lag, mean)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_vline(xintercept = seq(24, lag_max, by = 24), linetype = "dashed", color = "grey60") +
  geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd), alpha = 0.15) +
  geom_line(linewidth = 0.9) +
  labs(title = "Average PACF (deseasonalized series)", x = "Lag", y = "Mean PACF")

save_grid_png(grid.arrange(p_avg_acf_adj, p_avg_pacf_adj, ncol = 2), "avg_acf_pacf_stl", width = 14, height = 6)
grid.arrange(p_avg_acf_adj, p_avg_pacf_adj, ncol = 2)

# -----------------------------
# 9) 24-hour seasonal shape (hour-of-day) — KEEP ONLY THIS seasonal diagnostic
# -----------------------------
avg_season <- rowMeans(season_stl)

df_season_24 <- data.frame(
  hour = hour(air$datetime),
  seasonal = avg_season
) %>%
  group_by(hour) %>%
  summarise(mean_season = mean(seasonal), sd_season = sd(seasonal), .groups = "drop")

p_season_24 <- ggplot(df_season_24, aes(hour, mean_season)) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey60") +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = 0:23) +
  labs(
    title = "Average daily seasonal pattern (STL)",
    x = "Hour of day", y = "Seasonal effect"
  )

save_plot(p_season_24, "stl_seasonal_pattern_24h", width = 8, height = 7)
print(p_season_24)

# -----------------------------
# 10) ADF p-values (optional quick check) on deseasonalized series
# -----------------------------
adf_lags <- 8
adf_res <- lapply(seq_len(ncol(X_adj)), function(j) {
  out <- suppressWarnings(tryCatch(tseries::adf.test(X_adj[, j], k = adf_lags), error = function(e) NULL))
  data.frame(
    sensor = colnames(X_adj)[j],
    p_value = if (is.null(out)) NA_real_ else unname(out$p.value)
  )
}) %>% bind_rows()

p_adf <- ggplot(adf_res %>% filter(!is.na(p_value)), aes(p_value)) +
  geom_histogram(bins = 30, color = "white") +
  labs(
    title = "ADF p-values (after STL deseasonalization)",
    subtitle = "Uncorrected p-values; Null: unit root",
    x = "p-value", y = "Count"
  )

save_plot(p_adf, "adf_pvalues_after_stl", width = 10, height = 6)
print(p_adf)

# -----------------------------
# 11) AR order selection (LOCAL) on STL-deseasonalized series
# -----------------------------
p_max_ar  <- 10
ar_method <- "yw"   # "yw" or "burg"
criterion <- "BIC"  # "AIC" or "BIC"

p_hat <- rep(NA_integer_, ncol(X_adj))
phi_mat <- matrix(NA_real_, nrow = ncol(X_adj), ncol = p_max_ar)
colnames(phi_mat) <- paste0("lag", 1:p_max_ar)

for (i in seq_len(ncol(X_adj))) {
  x <- as.numeric(X_adj[, i])
  fit <- ar(x, order.max = p_max_ar, aic = TRUE, method = ar_method)
  
  aic_vec <- fit$aic
  k_vec   <- 0:p_max_ar
  n_eff   <- length(x)
  
  score <- if (criterion == "BIC") aic_vec + (log(n_eff) - 2) * k_vec else aic_vec
  p_sel <- which.min(score) - 1
  p_hat[i] <- p_sel
  
  if (p_sel > 0) {
    fit_fix <- ar(x, aic = FALSE, order.max = p_sel, method = ar_method)
    phi_mat[i, 1:p_sel] <- as.numeric(fit_fix$ar)
  }
  
  if (i %% 50 == 0) cat("AR selection", i, "/", ncol(X_adj), "\n")
}

ar_dist <- data.frame(sensor = colnames(X_adj), p_hat = p_hat)

p_hist <- ggplot(ar_dist %>% filter(!is.na(p_hat)), aes(p_hat)) +
  geom_histogram(binwidth = 1, boundary = -0.5, color = "white") +
  scale_x_continuous(breaks = 0:p_max_ar) +
  labs(
    title = paste0("Selected AR order per sensor (", criterion, ")"),
    x = "Selected lag order p", y = "Count"
  )

coef_long <- as.data.frame(phi_mat) %>%
  mutate(sensor = colnames(X_adj), p_hat = p_hat) %>%
  pivot_longer(starts_with("lag"), names_to = "lag", values_to = "phi") %>%
  mutate(lag = as.integer(sub("^lag", "", lag))) %>%
  filter(!is.na(phi), !is.na(p_hat), p_hat > 0, lag <= p_hat)

p_box <- ggplot(coef_long, aes(factor(lag), phi)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_boxplot(outlier.alpha = 0.2) +
  facet_wrap(~ p_hat, scales = "free_x") +
  labs(
    title = "AR coefficient distributions",
    x = "Lag", y = "AR coefficient"
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

save_grid_png(grid.arrange(p_hist, p_box, ncol = 2), "ar_order_and_coeffs_stl", width = 14, height = 6)
grid.arrange(p_hist, p_box, ncol = 2)

# -----------------------------
# 12) Save ONLY cleaned + STL-deseasonalized data (for VAR)
# -----------------------------
saveRDS(
  list(
    X             = X,
    X_adj         = X_adj,
    datetime      = air$datetime,
    sensor_table  = sensor_table,
    season_period = 24,
    note          = "STL deseasonalization (freq=24, s.window='periodic', robust=TRUE); seasonal component removed, then centered per sensor."
  ),
  file = "airbox_stl_deseasonalized.rds"
)

cat("\nSaved clean VAR input to airbox_stl_deseasonalized.rds\n")
