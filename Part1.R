############################################################
# Part 1 — Taiwan AirBox EDA + stationarity/seasonality checks
# Adjusted per your requests:
# - Keep AR order selection LOCAL (per series), no global p
# - ADF p-value plot: standard histogram (no log scale)
# - AR coefficient display: compare coefficients ONLY within the SAME selected p
#   (faceted boxplots by selected order p)
# - Remove outliers V29 and V70 after outlier analysis; keep mapping consistent
# - Save ONLY transformed series + datetime + sensor_table
############################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(lubridate)
  library(forecast)   # ggAcf, ggPacf, na.interp
  library(corrplot)
  library(maps)
  library(gridExtra)
  library(grid)
  library(tseries)    # adf.test
})

# -----------------------------
# 0) Global ggplot theme
# -----------------------------
theme_project <- function(base_size = 13) {
  theme_bw() +
    theme(
      plot.title = element_text(size = base_size + 3, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size, margin = margin(b = 8)),
      axis.title = element_text(size = base_size + 1),
      axis.text = element_text(size = base_size - 1),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 1),
      panel.grid.major = element_line(color = "grey85", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    )
}
theme_set(theme_project())

# -----------------------------
# 1) Load the data
# -----------------------------
air_raw <- fread("TaiwanAirBox032017.csv")  # columns: time, V1..V516
loc <- fread("locations032017.csv")         # columns: V1, latitude, longitude
setnames(loc, old = "V1", new = "sensor_id")

sensor_cols <- setdiff(names(air_raw), "time")
sensor_id_from_name <- as.integer(sub("^V", "", sensor_cols))

# -----------------------------
# 2) Build datetime index
# -----------------------------
start_time <- ymd_hms("2017-03-01 00:00:00", tz = "Asia/Taipei")

air <- air_raw %>%
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
cat("Initial: T =", Tn, "hours; K =", K, "sensors\n")

# -----------------------------
# 4) Locations table
# -----------------------------
sensor_table <- data.frame(
  sensor_id   = sensor_id_from_name,
  sensor_name = paste0("S", sensor_id_from_name)
) %>%
  left_join(loc, by = "sensor_id")

# -----------------------------
# 5) Descriptive summaries + correlations
# -----------------------------
raw_mean <- colMeans(X, na.rm = TRUE)
raw_sd   <- apply(X, 2, sd, na.rm = TRUE)

desc <- sensor_table %>%
  mutate(mean_raw = raw_mean, sd_raw = raw_sd)

print(desc %>% summarise(
  mean_of_means = mean(mean_raw, na.rm = TRUE),
  median_sd     = median(sd_raw, na.rm = TRUE),
  min_sd        = min(sd_raw, na.rm = TRUE),
  max_sd        = max(sd_raw, na.rm = TRUE)
))

X_scaled <- scale(X, center = TRUE, scale = TRUE)
C0_full <- cor(X_scaled, use = "pairwise.complete.obs")

off_diag <- C0_full[upper.tri(C0_full)]
avg_corr <- mean(off_diag, na.rm = TRUE)
cat("Average lag-0 correlation (off-diagonal):", avg_corr, "\n")

avg_corr_per_sensor <- (rowSums(C0_full, na.rm = TRUE) - 1) / (K - 1)
desc <- desc %>% mutate(avg_corr = avg_corr_per_sensor)

# -----------------------------
# 6) Outlier analysis (plot only), then remove V29 and V70
# -----------------------------
outlier_ids <- desc %>%
  arrange(mean_raw) %>%
  slice(1:3) %>%
  pull(sensor_id)

desc_plot <- desc %>% mutate(is_outlier = sensor_id %in% outlier_ids)

ggplot(desc_plot, aes(x = mean_raw, y = sd_raw)) +
  geom_point(aes(shape = is_outlier), size = 3, alpha = 0.6) +
  geom_text(
    data = desc %>% filter(sensor_id %in% outlier_ids),
    aes(label = sensor_id),
    hjust = -0.4, vjust = 0, size = 5
  ) +
  scale_shape_manual(
    values = c(`FALSE` = 19, `TRUE` = 24),
    labels = c("Normal series", "Low-mean series")
  ) +
  labs(
    title = "Mean vs standard deviation across series",
    x = "Mean", y = "Standard deviation", shape = NULL
  ) +
  theme(legend.position = "inside", legend.position.inside = c(0.87, 0.10))

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

ggplot(outlier_df, aes(x = datetime, y = value)) +
  geom_line() +
  facet_wrap(~ sensor_label, scales = "free_y", ncol = 1) +
  labs(title = "Outlier time series", x = NULL, y = "Value")

# --- Remove V29 and V70 after outlier inspection ---
remove_ids  <- c(29, 70)
remove_cols <- paste0("V", remove_ids)

sensor_cols <- setdiff(sensor_cols, remove_cols)
air <- air %>% select(datetime, all_of(sensor_cols))

sensor_id_from_name <- as.integer(sub("^V", "", sensor_cols))

X <- as.matrix(air[, ..sensor_cols])
colnames(X) <- paste0("S", sensor_id_from_name)

Tn <- nrow(X)
K  <- ncol(X)
cat("After removing V29 & V70: T =", Tn, "hours; K =", K, "sensors\n")

sensor_table <- data.frame(
  sensor_id   = sensor_id_from_name,
  sensor_name = paste0("S", sensor_id_from_name)
) %>%
  left_join(loc, by = "sensor_id")

raw_mean <- colMeans(X, na.rm = TRUE)
raw_sd   <- apply(X, 2, sd, na.rm = TRUE)

desc <- sensor_table %>% mutate(mean_raw = raw_mean, sd_raw = raw_sd)

X_scaled <- scale(X, center = TRUE, scale = TRUE)
C0_full <- cor(X_scaled, use = "pairwise.complete.obs")

# -----------------------------
# 7) Time series overlay plots (full + zoom)
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

grid.arrange(p_full, p_zoom, ncol = 2)

# -----------------------------
# 8) Correlation heatmap
# -----------------------------
corrplot(
  C0_full,
  method = "color",
  tl.pos = "n",
  order = "hclust",
  mar = c(0, 0, 2, 0),
  title = "Lag-0 correlation matrix (all sensors, clustered)"
)

# -----------------------------
# 9) Sensor maps (world + Taiwan)
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
    fill = "grey95", color = "grey60", linewidth = 0.2
  ) +
  geom_point(
    data = desc,
    aes(x = longitude, y = latitude, color = mean_raw),
    size = 1.5, alpha = 0.6
  ) +
  coord_quickmap(xlim = lon_lim, ylim = lat_lim) +
  scale_color_viridis_c(option = "C", name = "Mean value") +
  labs(title = "Series coordinates (global)", x = "Longitude", y = "Latitude") +
  theme_minimal()

p_taiwan <- ggplot() +
  geom_polygon(
    data = taiwan_map,
    aes(x = long, y = lat, group = group),
    fill = "grey95", color = "grey40", linewidth = 0.3
  ) +
  geom_point(
    data = desc,
    aes(x = longitude, y = latitude, color = mean_raw),
    size = 1.5, alpha = 0.6
  ) +
  coord_quickmap(xlim = c(119.5, 122), ylim = c(21.5, 25.5)) +
  scale_color_viridis_c(option = "C", name = "Mean value") +
  labs(title = "Series coordinates (zoomed)", x = "Longitude", y = "Latitude") +
  theme_minimal()

grid.arrange(p_global, p_taiwan, ncol = 2)

# -----------------------------
# 10) ACF/PACF subset + average ACF/PACF
# -----------------------------
acf_pacf_ids <- 1:6
lag_max <- 80

acf_plots <- lapply(acf_pacf_ids, function(sid) {
  j <- which(sensor_id_from_name == sid)
  ggAcf(X[, j], lag.max = lag_max) +
    labs(title = paste("Sensor", sid)) +
    theme_project(8)
})

pacf_plots <- lapply(acf_pacf_ids, function(sid) {
  j <- which(sensor_id_from_name == sid)
  ggPacf(X[, j], lag.max = lag_max) +
    labs(title = paste("Sensor", sid)) +
    theme_project(8)
})

grid.arrange(
  grobs = acf_plots,
  ncol = 3,
  top = textGrob("Autocorrelation functions (ACF)", gp = gpar(fontsize = 16, fontface = "bold"))
)

grid.arrange(
  grobs = pacf_plots,
  ncol = 3,
  top = textGrob("Partial autocorrelation functions (PACF)", gp = gpar(fontsize = 16, fontface = "bold"))
)

acf_mat <- sapply(seq_len(ncol(X)), function(j) {
  stats::acf(X[, j], lag.max = lag_max, plot = FALSE)$acf[, 1, 1]
})
pacf_mat <- sapply(seq_len(ncol(X)), function(j) {
  stats::pacf(X[, j], lag.max = lag_max, plot = FALSE)$acf
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

ggplot(df_acf, aes(x = lag, y = mean)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_vline(xintercept = seq(24, max(df_acf$lag), by = 24),
             linetype = "dashed", color = "grey60", linewidth = 1) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15) +
  geom_line(linewidth = 0.9) +
  labs(title = "Average ACF across all sensors", x = "Lag", y = "Mean ACF")

ggplot(df_pacf, aes(x = lag, y = mean)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_vline(xintercept = seq(24, max(df_pacf$lag), by = 24),
             linetype = "dashed", color = "grey60", linewidth = 1) +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15) +
  geom_line(linewidth = 0.9) +
  labs(title = "Average PACF across all sensors", x = "Lag", y = "Mean PACF")

# -----------------------------
# 11) ADF after seasonal differencing (ALL sensors), uncorrected
# -----------------------------
season_period <- 24
adf_lags <- 30

seasonal_adjust_24 <- function(x, m = 24) {
  x <- as.numeric(x)
  if (any(!is.finite(x))) x <- forecast::na.interp(x)
  diff(x, lag = m)
}

adf_pvalue <- function(y, k = 30) {
  out <- tryCatch(tseries::adf.test(y, k = k), error = function(e) e)
  if (inherits(out, "error")) return(list(stat = NA_real_, p = NA_real_, note = out$message))
  list(stat = unname(out$statistic), p = unname(out$p.value), note = NA_character_)
}

sensor_ids <- colnames(X)

unitroot_all <- lapply(seq_len(ncol(X)), function(j) {
  sid <- sensor_ids[j]
  x <- X[, j]
  y <- seasonal_adjust_24(x, m = season_period)
  
  if (length(y) <= (adf_lags + 10)) {
    return(data.frame(
      sensor_id = sid,
      col_index = j,
      n_raw = length(x),
      n_used = length(y),
      adf_stat = NA_real_,
      p_value = NA_real_,
      note = "too_short_after_seasonal_diff"
    ))
  }
  
  tst <- adf_pvalue(y, k = adf_lags)
  
  data.frame(
    sensor_id = sid,
    col_index = j,
    n_raw = length(x),
    n_used = length(y),
    adf_stat = tst$stat,
    p_value = tst$p,
    note = tst$note
  )
}) %>% bind_rows()

unitroot_all <- unitroot_all %>%
  mutate(
    reject_05 = !is.na(p_value) & p_value < 0.05,
    reject_01 = !is.na(p_value) & p_value < 0.01
  )

cat("\n=== ADF after 24h seasonal differencing (all sensors; uncorrected) ===\n")
print(unitroot_all %>%
        summarise(
          K = n(),
          tested = sum(!is.na(p_value)),
          rejected_05 = sum(reject_05, na.rm = TRUE),
          rejected_01 = sum(reject_01, na.rm = TRUE),
          n_too_short = sum(note == "too_short_after_seasonal_diff", na.rm = TRUE),
          n_errors = sum(!is.na(note) & note != "too_short_after_seasonal_diff", na.rm = TRUE)
        )
)

# Standard histogram (no log scaling)
ggplot(unitroot_all %>% filter(!is.na(p_value)), aes(x = p_value)) +
  geom_histogram(bins = 30, color = "white") +
  labs(
    title = "ADF p-values after seasonal differencing (m=24)",
    subtitle = "Uncorrected p-values; Null: unit root.",
    x = "p-value",
    y = "Count"
  )

# -----------------------------
# Build transformed matrix once (used for AR + saved for modeling)
# -----------------------------
X_seas24 <- apply(X, 2, function(x) {
  x <- as.numeric(x)
  if (any(!is.finite(x))) x <- forecast::na.interp(x)
  diff(x, lag = season_period)
})
colnames(X_seas24) <- colnames(X)
datetime_seas24 <- air$datetime[-seq_len(season_period)]

# -----------------------------
# 12) AR order selection (LOCAL) + better coefficient presentation
#     - Select p per sensor by criterion
#     - Show histogram of selected p
#     - Show coefficient distributions CONDITIONAL on selected p (facets)
# -----------------------------
p_max_ar  <- 10
ar_method <- "yw"        # "yw" or "burg"
criterion <- "BIC"       # "AIC" or "BIC"

K2 <- ncol(X_seas24)

p_hat <- rep(NA_integer_, K2)
phi_mat <- matrix(NA_real_, nrow = K2, ncol = p_max_ar)
colnames(phi_mat) <- paste0("lag", 1:p_max_ar)

for (i in seq_len(K2)) {
  x <- as.numeric(X_seas24[, i])
  if (any(!is.finite(x))) x <- forecast::na.interp(x)
  
  if (length(x) <= (p_max_ar + 10)) next
  
  fit <- ar(x, order.max = p_max_ar, aic = TRUE, method = ar_method)
  
  aic_vec <- fit$aic
  k_vec <- 0:p_max_ar
  n_eff <- length(x)
  
  score <- if (criterion == "BIC") {
    aic_vec + (log(n_eff) - 2) * k_vec
  } else {
    aic_vec
  }
  
  p_sel <- which.min(score) - 1
  p_hat[i] <- p_sel
  
  if (!is.na(p_sel) && p_sel > 0) {
    fit_fix <- ar(x, aic = FALSE, order.max = p_sel, method = ar_method)
    phi_mat[i, 1:p_sel] <- as.numeric(fit_fix$ar)
  }
  
  if (i %% 50 == 0) cat("AR selection:", i, "/", K2, "\n")
}

ar_dist <- data.frame(sensor = colnames(X_seas24), p_hat = p_hat)

# Long format coefficients + keep selected p
coef_long <- as.data.frame(phi_mat) %>%
  mutate(
    sensor = colnames(X_seas24),
    p_hat  = p_hat
  ) %>%
  pivot_longer(cols = starts_with("lag"), names_to = "lag", values_to = "phi") %>%
  mutate(lag = as.integer(sub("^lag", "", lag))) %>%
  filter(!is.na(phi), !is.na(p_hat), p_hat > 0, lag <= p_hat)

# Plot 1: histogram of selected orders
p_hist <- ggplot(ar_dist %>% filter(!is.na(p_hat)), aes(x = p_hat)) +
  geom_histogram(binwidth = 1, boundary = -0.5, color = "white") +
  scale_x_continuous(breaks = 0:p_max_ar) +
  labs(
    title = paste0("Selected AR order per sensor (", criterion, ")"),
    x = "Selected lag order p",
    y = "Count"
  )

# Plot 2: coefficients conditional on selected p (comparability within each facet)
# If too many facets: reduce by setting p_max_ar smaller or facet_wrap(ncol = ...)
p_box <- ggplot(coef_long, aes(x = factor(lag), y = phi)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_boxplot(outlier.alpha = 0.2) +
  facet_wrap(~ p_hat, scales = "free_x") +
  labs(
    title = "AR coefficient distributions",
    subtitle = "Within each facet p, coefficients are comparable.",
    x = "Lag",
    y = "AR coefficient"
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

grid.arrange(p_hist, p_box, ncol = 2)

# -----------------------------
# Save ONLY transformed data + metadata for modeling
# -----------------------------
saveRDS(
  list(
    X_seas24      = X_seas24,
    datetime      = datetime_seas24,
    sensor_table  = sensor_table,
    season_period = season_period
  ),
  file = "airbox_seasonally_adjusted.rds"
)

cat("\nSaved clean modeling input to airbox_seasonally_adjusted.rds\n")
