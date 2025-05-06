simulate_dtr_data <- function(N, T = 30, seed = 123) {
  set.seed(seed)
  df <- data.frame(id = 1:N, M = rnorm(N, 0, 1))
  df$S_H1 <- rnorm(N, 6, 1)
  df$S_B1 <- rnorm(N, 30, 2)
  X <- Z <- matrix(0, N, T)
  S_H <- S_B <- matrix(0, N, T + 1)
  S_H[,1] <- df$S_H1; S_B[,1] <- df$S_B1
  alpha0 <- 0; alpha1 <- -0.4; alpha2 <- 0.2
  beta_H0 <- 0.05; beta_H1 <- -0.20; beta_Hm <- -0.10
  beta_B0 <- 0.02; beta_B1 <- 0.10; beta_Bm <- -0.05
  sigma_H <- 0.1; sigma_B <- 0.2; sigma_Z <- 0.5
  for (t in 1:T) {
    Z[,t]   <- df$M * 0.8 + rnorm(N, sd = sigma_Z)
    logit_p <- alpha0 + alpha1 * S_H[,t] + alpha2 * S_B[,t]
    p_t     <- 1 / (1 + exp(-logit_p))
    X[,t]   <- rbinom(N, 1, p_t)
    S_H[,t+1] <- S_H[,t] + beta_H0 + beta_H1 * X[,t] + beta_Hm * df$M + rnorm(N, sigma_H)
    S_B[,t+1] <- S_B[,t] + beta_B0 + beta_B1 * X[,t] + beta_Bm * df$M + rnorm(N, sigma_B)
  }
  sim_df <- df
  for (t in 1:(T+1)) sim_df[[paste0("S_H",t)]] <- S_H[,t]
  for (t in 1:(T+1)) sim_df[[paste0("S_B",t)]] <- S_B[,t]
  for (t in 1:T)    sim_df[[paste0("Z",t)]]  <- Z[,t]
  for (t in 1:T)    sim_df[[paste0("X",t)]]  <- X[,t]
  sim_df$Y <- -sim_df[[paste0("S_H", T+1)]]
  sim_df
}

weighted_qlearning <- function(df, T = 30, gamma = 1.0) {
  n <- nrow(df)
  mu_H <- sd_H <- mu_B <- sd_B <- mu_Z <- sd_Z <- vector("list", T)
  for (t in 1:T) {
    mu_H[[t]] <- mean(df[[paste0("S_H",t)]])
    sd_H[[t]] <- sd(df[[paste0("S_H",t)]])
    mu_B[[t]] <- mean(df[[paste0("S_B",t)]])
    sd_B[[t]] <- sd(df[[paste0("S_B",t)]])
    mu_Z[[t]] <- mean(df[[paste0("Z",t)]])
    sd_Z[[t]] <- sd(df[[paste0("Z",t)]])
    df[[paste0("sH",t)]] <- (df[[paste0("S_H",t)]] - mu_H[[t]]) / sd_H[[t]]
    df[[paste0("sB",t)]] <- (df[[paste0("S_B",t)]] - mu_B[[t]]) / sd_B[[t]]
    df[[paste0("sZ",t)]] <- (df[[paste0("Z",t)]]   - mu_Z[[t]]) / sd_Z[[t]]
  }
  w <- matrix(NA, n, T)
  betas <- vector("list", T); V <- matrix(0, n, T+1)
  for (t in 1:T) {
    ps_mod <- glm(df[[paste0("X",t)]] ~ df[[paste0("sH",t)]] + df[[paste0("sB",t)]], family = binomial)
    p1     <- predict(ps_mod, type = "response")
    p_i    <- ifelse(df[[paste0("X",t)]] == 1, p1, 1 - p1)
    w[,t]  <- 1 / p_i
  }
  for (t in T:1) {
    X_t <- df[[paste0("X",t)]]
    sH  <- df[[paste0("sH",t)]]; sB <- df[[paste0("sB",t)]]; sZ <- df[[paste0("sZ",t)]]
    Phi <- cbind(1, sH, sB, sZ, X_t, sH * X_t, sB * X_t, sZ * X_t)
    y_t <- if (t == T) df$Y else gamma * V[,t+1]
    W   <- diag(w[,t]); beta <- solve(t(Phi) %*% W %*% Phi, t(Phi) %*% W %*% y_t)
    betas[[t]] <- beta; V[,t] <- pmax(Phi[,1:4] %*% beta[1:4], Phi %*% beta)
  }
  list(betas = betas, V = V)
}

library(dplyr)
H_grid <- 4:12
B_grid <- seq(15, 45, 2)
Z_grid <- seq(0, 10, 2)
state_map_df <- expand.grid(H = H_grid, B = B_grid, Z = Z_grid) %>%
  mutate(s_idx = row_number())

make_pi_matrix <- function(result, state_map_df, T = 30) {
  mu_H <- sapply(1:T, function(t) mean(sim_df[[paste0("S_H", t)]]))
  sd_H <- sapply(1:T, function(t) sd(sim_df[[paste0("S_H", t)]]))
  mu_B <- sapply(1:T, function(t) mean(sim_df[[paste0("S_B", t)]]))
  sd_B <- sapply(1:T, function(t) sd(sim_df[[paste0("S_B", t)]]))
  mu_Z <- sapply(1:T, function(t) mean(sim_df[[paste0("Z",   t)]]))
  sd_Z <- sapply(1:T, function(t) sd(sim_df[[paste0("Z",   t)]]))
  pi_matrix <- matrix(NA, nrow(state_map_df), T)
  for (t in 1:T) {
    beta <- result$betas[[t]]
    df_grid <- state_map_df %>%
      mutate(
        sH = (H - mu_H[t]) / sd_H[t],
        sB = (B - mu_B[t]) / sd_B[t],
        sZ = (Z - mu_Z[t]) / sd_Z[t]
      )
    Q0 <- as.numeric(cbind(1, df_grid$sH, df_grid$sB, df_grid$sZ, 0, 0, 0, 0) %*% beta)
    Q1 <- as.numeric(cbind(1, df_grid$sH, df_grid$sB, df_grid$sZ, 1, df_grid$sH, df_grid$sB, df_grid$sZ) %*% beta)
    pi_matrix[,t] <- ifelse(Q1 > Q0, 1, 0)
  }
  pi_matrix
}

library(ggplot2); library(dplyr); library(patchwork)

plot_policy_slice <- function(stage_t, z_idx, pi_matrix, state_map_df,
                              h_bins, b_bins, z_bins) {
  if (stage_t < 1 || stage_t > ncol(pi_matrix)) stop("Invalid stage_t")
  if (z_idx   < 1 || z_idx   > length(z_bins)) stop("Invalid z_idx")
  policy_for_stage <- data.frame(s_idx = seq_len(nrow(pi_matrix)), Action = pi_matrix[, stage_t])
  policy_data <- left_join(policy_for_stage, state_map_df, by = "s_idx")
  policy_slice <- policy_data %>%
    filter(Z == z_bins[z_idx]) %>%
    mutate(
      H_center      = H,
      B_center      = B,
      Action_Factor = factor(Action, levels = c(0,1), labels = c("Lifestyle","Medication"))
    )
  req_levs <- c("Lifestyle","Medication")
  pres_levs <- unique(policy_slice$Action_Factor)
  dummy_rows <- list()
  dh <- min(h_bins) - (max(h_bins)-min(h_bins))
  db <- min(b_bins) - (max(b_bins)-min(b_bins))
  if (!"Lifestyle" %in% pres_levs)  dummy_rows[["Lifestyle"]]  <- data.frame(H_center=dh, B_center=db, Action_Factor=factor("Lifestyle", req_levs))
  if (!"Medication" %in% pres_levs) dummy_rows[["Medication"]] <- data.frame(H_center=dh, B_center=db, Action_Factor=factor("Medication", req_levs))
  plot_data_final <- if (length(dummy_rows)) bind_rows(policy_slice, bind_rows(dummy_rows)) else policy_slice
  
  ggplot(plot_data_final, aes(x = B_center, y = H_center, fill = Action_Factor)) +
    geom_tile(color = "grey90", linewidth = 0.1) +
    scale_fill_manual(values = c("Lifestyle" = "lightblue", "Medication" = "salmon"),
                      name = "Optimal Action", drop = FALSE) +
    scale_x_continuous(expand = c(0,0), breaks = scales::pretty_breaks()) +
    scale_y_continuous(expand = c(0,0), breaks = scales::pretty_breaks()) +
    coord_cartesian(xlim = range(b_bins), ylim = range(h_bins), expand = FALSE) +
    labs(title = sprintf("Stage %d, Z=%.1f", stage_t, z_bins[z_idx]), x = "BMI", y = "HbA1c") +
    theme_minimal(base_size = 10) +
    theme(
      plot.title      = element_text(hjust = 0.5, size = 12),
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y     = element_text(size = 8)
    )
}

stages_to_plot <- c(1, 15, 30)
z_values       <- c(0, 4, 10)
z_indices      <- match(z_values, Z_grid)
sim_df       <- simulate_dtr_data(2000)
result       <- weighted_qlearning(sim_df)
pi_matrix    <- make_pi_matrix(result, state_map_df)

each_plot <- list()
for (t in stages_to_plot) {
  for (zi in z_indices) {
    if (!is.na(zi)) {
      each_plot[[paste0("S",t,"_Z",z_values[which(z_indices==zi)])]] <-
        plot_policy_slice(t, zi, pi_matrix, state_map_df, H_grid, B_grid, Z_grid)
    }
  }
}

combined_plot <- wrap_plots(each_plot,
                            ncol = length(z_values),
                            widths = rep(1.5, length(z_values)),
                            heights = rep(1, length(stages_to_plot)),
                            guides = "collect"
) & theme(legend.position = "bottom")

final_plot <- combined_plot +
  plot_annotation(title = "Optimal Policy at Stages 1,15,30 × Z=0,4,10")

print(final_plot)

beta30 <- betas[[30]]

phi0_x30 <- c(1, 8.087, 0.627, 1.364, 0, 0, 0, 0)
phi1_x30 <- c(1, 8.087, 0.627, 1.364, 1, 8.087, 0.627, 1.364)

Q0_x30 <- sum(phi0_x30 * beta30)
Q1_x30 <- sum(phi1_x30 * beta30)
V30_x  <- max(Q0_x30, Q1_x30)

cat(sprintf("Q30(sx,0) = %.4f\n", Q0_x30))
cat(sprintf("Q30(sx,1) = %.4f\n", Q1_x30))
cat(sprintf("V30(sx)   = %.4f\n\n", V30_x))

phi0_y30 <- c(1, 5.528, 0.494, 3.000, 0, 0, 0, 0)
phi1_y30 <- c(1, 5.528, 0.494, 3.000, 1, 5.528, 0.494, 3.000)

Q0_y30 <- sum(phi0_y30 * beta30)
Q1_y30 <- sum(phi1_y30 * beta30)
V30_y  <- max(Q0_y30, Q1_y30)

cat(sprintf("Q30(sy,0) = %.4f\n", Q0_y30))
cat(sprintf("Q30(sy,1) = %.4f\n", Q1_y30))
cat(sprintf("V30(sy)   = %.4f\n\n", V30_y))

Ystar_29 <- V30_y
cat(sprintf("Y*_{29}   = %.4f\n\n", Ystar_29))

beta29 <- betas[[29]]

phi0_x29 <- c(1, 8.246, 0.596, 1.366, 0, 0, 0, 0)
phi1_x29 <- c(1, 8.246, 0.596, 1.366, 1, 8.246, 0.596, 1.366)

Q0_x29 <- sum(phi0_x29 * beta29)
Q1_x29 <- sum(phi1_x29 * beta29)
V29_x  <- max(Q0_x29, Q1_x29)

cat(sprintf("Q29(sx,0) = %.4f\n", Q0_x29))
cat(sprintf("Q29(sx,1) = %.4f\n", Q1_x29))
cat(sprintf("V29(sx)   = %.4f\n", V29_x))
