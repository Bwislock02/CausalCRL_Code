library(dplyr)
library(tidyr)
library(pbapply)
library(ggplot2)
library(patchwork)

set.seed(123)

n_stages <- 30
gamma_discount <- 0.99
medication_cost <- 1

initial_mean_hba1c <- 8.0
initial_sd_hba1c <- 1.0
initial_mean_bmi <- 30.0
initial_sd_bmi <- 5.0

COEFFS <- list(
  A = list(prob_M1 = 0.8, prob_M0_X1 = 0.4, prob_M0_X0 = 0.6),
  Z = list(beta_M = 3, beta_H = 0.2, beta_B = -0.05, beta_0 = 0, noise_sd = 1.0),
  S_hba1c = list(beta_S = 1.0, delta_XA_11 = -1.5, delta_XA_10 = -0.2,
                 delta_XA_01 = -0.5, delta_XA_00 = +0.1, beta_BMI = 0.12,
                 intercept_offset = -0.05, noise_sd = 0.5, min_val = 4.0),
  S_bmi = list(beta_S = 1.0, delta_X_1 = 0.1, delta_X_0_A1 = -0.3,
               delta_X_0_A0 = 0.0, beta_HbA1c = 0.02,
               intercept_offset = -0.06, noise_sd = 0.5, min_val = 15.0)
)

hba1c_min <- 4.0; hba1c_max <- 12.0; hba1c_step <- 1.0
bmi_min <- 15.0; bmi_max <- 45.0; bmi_step <- 2.0
z_bins   <- c(0, 4, 10)
z_min    <- min(z_bins)
z_max    <- max(z_bins)
n_z_bins <- length(z_bins)

hba1c_bins <- seq(hba1c_min, hba1c_max, by = hba1c_step)
bmi_bins   <- seq(bmi_min, bmi_max, by = bmi_step)

n_h_bins <- length(hba1c_bins)
n_b_bins <- length(bmi_bins)
n_states_total <- n_h_bins * n_b_bins * n_z_bins
n_actions <- 2

cat(sprintf("MDP Setup:\n N Stages: %d\n Gamma: %.2f\n Medication Cost: %.2f\n N States: %d (%d H x %d B x %d Z)\n N Actions: %d\n",
            n_stages, gamma_discount, medication_cost, n_states_total, n_h_bins, n_b_bins, n_z_bins, n_actions))

find_bin_index <- function(value, bins) {
  value_clamped <- pmax(min(bins), pmin(max(bins), value))
  which.min(abs(bins - value_clamped))
}

get_discrete_state_indices <- function(s_hba1c, s_bmi, z, h_bins, b_bins, z_bins) {
  h_idx <- find_bin_index(s_hba1c, h_bins)
  b_idx <- find_bin_index(s_bmi, b_bins)
  z_idx <- find_bin_index(z, z_bins)
  return(list(h = h_idx, b = b_idx, z = z_idx))
}

get_bin_centers <- function(h_idx, b_idx, z_idx, h_bins, b_bins, z_bins) {
  h_idx <- max(1, min(n_h_bins, h_idx))
  b_idx <- max(1, min(n_b_bins, b_idx))
  z_idx <- max(1, min(n_z_bins, z_idx))
  return(list(h = h_bins[h_idx], b = b_bins[b_idx], z = z_bins[z_idx]))
}

get_1d_index <- function(h_idx, b_idx, z_idx, n_h_bins, n_b_bins) {
  h_idx <- max(1, min(n_h_bins, h_idx))
  b_idx <- max(1, min(n_b_bins, b_idx))
  z_idx <- max(1, min(n_z_bins, z_idx))
  return(as.integer((z_idx - 1) * (n_h_bins * n_b_bins) + (b_idx - 1) * n_h_bins + h_idx))
}

state_map <- expand.grid(h_idx = 1:n_h_bins, b_idx = 1:n_b_bins, z_idx = 1:n_z_bins)
state_map$s_idx <- 1:n_states_total

get_3d_indices <- function(s_idx, state_map) {
  s_idx <- max(1, min(nrow(state_map), s_idx))
  return(state_map[s_idx, c("h_idx", "b_idx", "z_idx")])
}

simulate_single_step <- function(current_hba1c, current_bmi, current_z,
                                 action_X, motivation_M, coeffs) {
  if (any(is.na(c(current_hba1c, current_bmi, current_z, action_X, motivation_M)))) {
    stop(sprintf("NA input: H=%.1f, B=%.1f, Z=%.1f, X=%d, M=%d",
                 current_hba1c, current_bmi, current_z, action_X, motivation_M))
  }
  prob_adherence <- if (motivation_M == 1) { coeffs$A$prob_M1 }
  else if (action_X == 1) { coeffs$A$prob_M0_X1 }
  else { coeffs$A$prob_M0_X0 }
  prob_adherence <- max(0, min(1, prob_adherence))
  a_t <- rbinom(1, 1, prob_adherence)
  
  z_normalized <- pmin(1, pmax(0, current_z / z_max))
  z_effect_multiplier <- 1 + 0.4 * (z_normalized - 0.5)
  
  delta_hba1c <- if (action_X == 1 && a_t == 1) {
    coeffs$S_hba1c$delta_XA_11 * z_effect_multiplier
  } else if (action_X == 1 && a_t == 0) {
    coeffs$S_hba1c$delta_XA_10
  } else if (action_X == 0 && a_t == 1) {
    coeffs$S_hba1c$delta_XA_01
  } else {
    coeffs$S_hba1c$delta_XA_00
  }
  
  bmi_effect_on_hba1c <- coeffs$S_hba1c$beta_BMI * (current_bmi - initial_mean_bmi)
  noise_hba1c <- rnorm(1, 0, coeffs$S_hba1c$noise_sd)
  raw_hba1c <- coeffs$S_hba1c$beta_S * current_hba1c + delta_hba1c +
    bmi_effect_on_hba1c + coeffs$S_hba1c$intercept_offset + noise_hba1c
  next_hba1c <- pmax(raw_hba1c, coeffs$S_hba1c$min_val, na.rm = TRUE)
  
  delta_bmi <- if (action_X == 1) { coeffs$S_bmi$delta_X_1 }
  else if (action_X == 0 && a_t == 1) { coeffs$S_bmi$delta_X_0_A1 }
  else { coeffs$S_bmi$delta_X_0_A0 }
  hba1c_effect_on_bmi <- coeffs$S_bmi$beta_HbA1c * (current_hba1c - initial_mean_hba1c)
  noise_bmi <- rnorm(1, 0, coeffs$S_bmi$noise_sd)
  raw_bmi <- coeffs$S_bmi$beta_S * current_bmi + delta_bmi +
    hba1c_effect_on_bmi + coeffs$S_bmi$intercept_offset + noise_bmi
  next_bmi <- pmax(raw_bmi, coeffs$S_bmi$min_val, na.rm = TRUE)
  
  lin_comb_Z <- coeffs$Z$beta_M * motivation_M + coeffs$Z$beta_H * current_hba1c +
    coeffs$Z$beta_B * current_bmi + coeffs$Z$beta_0
  noise_Z <- rnorm(1, 0, coeffs$Z$noise_sd)
  next_Z <- round(lin_comb_Z + noise_Z)
  next_Z <- pmax(z_min, pmin(z_max, next_Z))
  
  next_hba1c <- ifelse(!is.finite(next_hba1c), coeffs$S_hba1c$min_val, next_hba1c)
  next_bmi <- ifelse(!is.finite(next_bmi), coeffs$S_bmi$min_val, next_bmi)
  next_Z <- ifelse(!is.finite(next_Z), (z_min + z_max) / 2.0, next_Z)
  
  return(list(next_S_hba1c=next_hba1c, next_S_bmi=next_bmi, next_Z=next_Z))
}

n_sim_per_state <- 10000
P_estimate_filename <- sprintf("transition_prob_sims%d.rds", n_sim_per_state)

P_file <- P_estimate_filename

if (file.exists(P_file)) {
  P <- readRDS(P_file)
  cat("Loaded P from", P_file, "\n")
  if(!all(dim(P) == c(n_states_total, n_actions, n_states_total))) {
    stop("Dimensions of loaded P matrix do not match current configuration.")
  }
} else {
  P <- array(0, dim = c(n_states_total, n_actions, n_states_total))
  cat(sprintf("Estimating Transition Probabilities (N_sim/state = %d). This will take time...\n", n_sim_per_state))
  pb <- txtProgressBar(min = 0, max = n_states_total * n_actions, style = 3)
  progress_counter <- 0
  
  for (s_idx in 1:n_states_total) {
    current_indices <- get_3d_indices(s_idx, state_map)
    h_idx <- current_indices$h_idx
    b_idx <- current_indices$b_idx
    z_idx <- current_indices$z_idx
    
    s_centers <- get_bin_centers(h_idx, b_idx, z_idx, hba1c_bins, bmi_bins, z_bins)
    s_hba1c_cont <- s_centers$h
    s_bmi_cont    <- s_centers$b
    s_z_cont      <- s_centers$z
    
    for (action in 0:(n_actions - 1)) {
      action_idx <- action + 1
      next_state_counts <- numeric(n_states_total)
      
      for (sim in 1:n_sim_per_state) {
        motivation_M <- rbinom(1, 1, 0.5)
        next_state_cont <- simulate_single_step(s_hba1c_cont, s_bmi_cont, s_z_cont,
                                                action, motivation_M, COEFFS)
        next_state_disc_indices <- get_discrete_state_indices(
          next_state_cont$next_S_hba1c, next_state_cont$next_S_bmi, next_state_cont$next_Z,
          hba1c_bins, bmi_bins, z_bins
        )
        s_prime_idx <- get_1d_index(
          next_state_disc_indices$h, next_state_disc_indices$b,
          next_state_disc_indices$z, n_h_bins, n_b_bins
        )
        
        if (s_prime_idx >= 1 && s_prime_idx <= n_states_total) {
          next_state_counts[s_prime_idx] <- next_state_counts[s_prime_idx] + 1
        } else {
          warning(sprintf("Invalid s_prime_idx %d generated for s_idx %d, action %d",
                          s_prime_idx, s_idx, action))
        }
      }
      
      probabilities <- next_state_counts / n_sim_per_state
      if (sum(probabilities) > 0) {
        probabilities <- probabilities / sum(probabilities)
      } else {
        probabilities[s_idx] <- 1
        warning(sprintf("Zero counts for s_idx %d, action %d. Setting self-loop.",
                        s_idx, action))
      }
      P[s_idx, action_idx, ] <- probabilities
      
      progress_counter <- progress_counter + 1
      setTxtProgressBar(pb, progress_counter)
    }
  }
  close(pb)
  cat("Transition Probability Estimation Complete.\n")
  saveRDS(P, P_estimate_filename)
}

state_map <- expand.grid(h_idx = 1:n_h_bins, b_idx = 1:n_b_bins, z_idx = 1:n_z_bins)
state_map$s_idx <- 1:n_states_total
get_3d_indices <- function(s_idx, state_map) {
  s_idx <- max(1, min(nrow(state_map), s_idx))
  return(state_map[s_idx, c("h_idx", "b_idx", "z_idx")])
}

R_base <- numeric(n_states_total)
for (s_prime_idx in 1:n_states_total) {
  indices <- get_3d_indices(s_prime_idx, state_map)
  centers <- get_bin_centers(indices$h_idx, indices$b_idx, indices$z_idx,
                             hba1c_bins, bmi_bins, z_bins)
  R_base[s_prime_idx] <- -centers$h
}
cat("Base Reward function R(s') defined.\n")

V <- matrix(0, nrow = n_states_total, ncol = n_stages + 1)
Pi <- matrix(NA, nrow = n_states_total, ncol = n_stages)

cat("Starting Value Iteration with Medication Cost...\n")
pb <- txtProgressBar(min = 0, max = n_stages, style = 3)

for (t in n_stages:1) {
  V_next <- V[, t + 1]
  term_in_brackets <- R_base + gamma_discount * V_next
  
  if(any(!is.finite(term_in_brackets))) {
    warning(paste("Non-finite values found in term_in_brackets at stage", t,
                  ". These will propagate. Ensure Q-value handling is robust."))
  }
  
  P_a0 <- P[, 1, ]
  P_a1 <- P[, 2, ]
  
  Q_t_a0 <- P_a0 %*% term_in_brackets
  Q_t_a1 <- (P_a1 %*% term_in_brackets) - medication_cost
  
  if(any(!is.finite(Q_t_a0))) {
    warning(paste("Non-finite values found in Q_t_a0 at stage", t))
    Q_t_a0[!is.finite(Q_t_a0)] <- -Inf
  }
  if(any(!is.finite(Q_t_a1))) {
    warning(paste("Non-finite values found in Q_t_a1 at stage", t))
    Q_t_a1[!is.finite(Q_t_a1)] <- -Inf
  }
  
  V[, t] <- pmax(Q_t_a0, Q_t_a1, na.rm = TRUE)
  valid_q <- is.finite(Q_t_a0) | is.finite(Q_t_a1)
  Pi[valid_q, t] <- ifelse(Q_t_a1[valid_q] > Q_t_a0[valid_q], 1, 0)
  Pi[!valid_q, t] <- 0 
  
  setTxtProgressBar(pb, n_stages - t + 1)
}
close(pb)
cat("Value Iteration Complete.\n")

optimal_policy_stage1 <- Pi[, 1]
cat("\n--- Example: Optimal Policy for Stage 1 ---\n")
example_s_indices <- c(1, round(n_states_total / 2), n_states_total)
example_s_indices <- example_s_indices[example_s_indices >= 1 & example_s_indices <= n_states_total]

for (s_idx in example_s_indices) {
  indices_3d <- get_3d_indices(s_idx, state_map)
  centers_3d <- get_bin_centers(indices_3d$h_idx, indices_3d$b_idx, indices_3d$z_idx,
                                hba1c_bins, bmi_bins, z_bins)
  action <- Pi[s_idx, 1]
  if (is.na(action)) {
    action_desc <- "NA"
    action <- -1
  } else {
    action_desc <- ifelse(action == 1, "Medication", "Lifestyle")
  }
  cat(sprintf(" State s=%d (H~%.1f, B~%.1f, Z~%.1f): Optimal Action at t=1 -> %s (%d)\n",
              s_idx, centers_3d$h, centers_3d$b, centers_3d$z, action_desc, action))
}

optimal_policy_stageN <- Pi[, n_stages]
cat("\n--- Example: Optimal Policy for Stage N ---\n")
for (s_idx in example_s_indices) {
  indices_3d <- get_3d_indices(s_idx, state_map)
  centers_3d <- get_bin_centers(indices_3d$h_idx, indices_3d$b_idx, indices_3d$z_idx,
                                hba1c_bins, bmi_bins, z_bins)
  action <- Pi[s_idx, n_stages]
  if (is.na(action)) {
    action_desc <- "NA"
    action <- -1
  } else {
    action_desc <- ifelse(action == 1, "Medication", "Lifestyle")
  }
  cat(sprintf(" State s=%d (H~%.1f, B~%.1f, Z~%.1f): Optimal Action at t=%d -> %s (%d)\n",
              s_idx, centers_3d$h, centers_3d$b, centers_3d$z, n_stages, action_desc, action))
}

value_function_stage1 <- V[, 1]
cat("\n--- Example: Optimal Value Function for Stage 1 ---\n")
for (s_idx in example_s_indices) {
  indices_3d <- get_3d_indices(s_idx, state_map)
  centers_3d <- get_bin_centers(indices_3d$h_idx, indices_3d$b_idx, indices_3d$z_idx,
                                hba1c_bins, bmi_bins, z_bins)
  value <- V[s_idx, 1]
  cat(sprintf(" State s=%d (H~%.1f, B~%.1f, Z~%.1f): Optimal Value V*(s) at t=1 = %.3f\n",
              s_idx, centers_3d$h, centers_3d$b, centers_3d$z, value))
}

cat("\nDone.\n")

state_map_df <- state_map %>%
  dplyr::select(s_idx, h_idx, b_idx, z_idx)

plot_policy_slice <- function(stage_t, z_idx, pi_matrix, state_map_df,
                              hba1c_bins, bmi_bins, z_bins) {
  if (stage_t < 1 || stage_t > ncol(pi_matrix)) {
    stop("Invalid stage_t provided.")
  }
  if (z_idx < 1 || z_idx > length(z_bins)) {
    stop("Invalid z_idx provided.")
  }
  
  policy_data <- data.frame(
    s_idx = 1:nrow(pi_matrix),
    Action = pi_matrix[, stage_t]
  ) %>%
    left_join(state_map_df, by = "s_idx") %>%
    filter(z_idx == !!z_idx) %>%
    mutate(
      H_center     = hba1c_bins[h_idx],
      B_center     = bmi_bins[b_idx],
      Action_Factor = factor(
        Action,
        levels = c(0, 1),
        labels = c("Lifestyle", "Medication")
      )
    )
  
  required_levels <- c("Lifestyle", "Medication")
  present_levels  <- unique(policy_data$Action_Factor)
  dummy_rows <- list()
  dummy_h <- min(hba1c_bins) - diff(range(hba1c_bins))
  dummy_b <- min(bmi_bins)   - diff(range(bmi_bins))
  
  if (!("Lifestyle" %in% present_levels)) {
    dummy_rows[[1]] <- data.frame(
      H_center = dummy_h,
      B_center = dummy_b,
      Action_Factor = factor("Lifestyle", levels = required_levels)
    )
  }
  if (!("Medication" %in% present_levels)) {
    dummy_rows[[2]] <- data.frame(
      H_center = dummy_h,
      B_center = dummy_b,
      Action_Factor = factor("Medication", levels = required_levels)
    )
  }
  if (length(dummy_rows)) {
    policy_data <- bind_rows(policy_data, bind_rows(dummy_rows))
  }
  
  z_val <- z_bins[z_idx]
  ggplot(policy_data, aes(x = B_center, y = H_center, fill = Action_Factor)) +
    geom_tile(color = "grey90", linewidth = 0.1) +
    scale_fill_manual(
      values = c("Lifestyle" = "lightblue", "Medication" = "salmon"),
      name   = "Action",
      drop   = FALSE
    ) +
    coord_cartesian(
      xlim   = range(bmi_bins),
      ylim   = range(hba1c_bins),
      expand = FALSE
    ) +
    labs(
      title = sprintf("Policy Heatmap — Stage %d, Z = %.1f", stage_t, z_val),
      x     = "BMI",
      y     = "HbA1c"
    ) +
    theme_minimal() +
    theme(
      plot.title       = element_text(hjust = 0.5),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      legend.position  = "none"
    )
}

stages_to_plot    <- c(1, 15, ncol(Pi))
z_indices_to_plot <- seq_len(min(3, length(z_bins)))

plot_list <- list()
for (t in stages_to_plot) {
  for (z in z_indices_to_plot) {
    plot_list[[paste0("t", t, "_z", z)]] <-
      plot_policy_slice(
        stage_t      = t,
        z_idx        = z,
        pi_matrix    = Pi,
        state_map_df = state_map_df,
        hba1c_bins   = hba1c_bins,
        bmi_bins     = bmi_bins,
        z_bins       = z_bins
      )
  }
}

combined_plot <- wrap_plots(plot_list, 
                            ncol   = length(z_indices_to_plot),
                            guides = "collect") +
  plot_annotation(title = "Optimal Policy Across Stages & Z") &
  theme(legend.position = "bottom")

print(combined_plot)

fixed_h_idx <- 7
fixed_b_idx <- round(n_b_bins / 2)
fixed_z_idx <- round(n_z_bins / 2)

cat(sprintf("Extracting probabilities starting from state (h=%d, b=%d, z=%d)\n",
            fixed_h_idx, fixed_b_idx, fixed_z_idx))
cat(sprintf(" -> Corresponding Centers: (H~%.1f, B~%.1f, Z~%.1f)\n",
            hba1c_bins[fixed_h_idx], bmi_bins[fixed_b_idx], z_bins[fixed_z_idx]))

results_list <- list()
start_s_idx <- get_1d_index(fixed_h_idx, fixed_b_idx, fixed_z_idx, n_h_bins, n_b_bins)

if(start_s_idx < 1 || start_s_idx > n_states_total) {
  stop(paste("Invalid start_s_idx calculated:", start_s_idx))
}

cat("Processing...\n")
for (action_a in 0:1) {
  action_idx <- action_a + 1
  
  prob_next_h <- numeric(n_h_bins)
  names(prob_next_h) <- 1:n_h_bins
  
  for (s_prime_idx in 1:n_states_total) {
    s_prime_indices <- get_3d_indices(s_prime_idx, state_map)
    h_prime <- s_prime_indices$h_idx
    
    if (h_prime >= 1 && h_prime <= n_h_bins) {
      prob_s_prime <- P[start_s_idx, action_idx, s_prime_idx]
      prob_next_h[h_prime] <- prob_next_h[h_prime] + prob_s_prime
    } else {
      warning(paste("Invalid h_prime index", h_prime, "encountered for s_prime_idx", s_prime_idx))
    }
  }
  
  total_prob_h = sum(prob_next_h)
  if (total_prob_h > 1e-9) {
    prob_next_h <- prob_next_h / total_prob_h
  } else {
    warning(paste("Total probability is near zero for start state", start_s_idx, "action", action_a))
    prob_next_h[fixed_h_idx] <- 1.0
  }
  
  for (next_h_idx in 1:n_h_bins) {
    results_list[[length(results_list) + 1]] <- data.frame(
      Start_H_Idx = fixed_h_idx,
      Start_H_Center = hba1c_bins[fixed_h_idx],
      Fixed_B_Idx = fixed_b_idx,
      Fixed_Z_Idx = fixed_z_idx,
      Action = action_a,
      Action_Desc = ifelse(action_a == 0, "Lifestyle", "Medication"),
      Next_H_Idx = next_h_idx,
      Next_H_Center = hba1c_bins[next_h_idx],
      Probability = round(prob_next_h[next_h_idx], 5)
    )
  }
}

prob_single_start_df <- bind_rows(results_list)
cat("Probability extraction complete.\n")
print(prob_single_start_df)

table_a0 <- prob_single_start_df %>%
  filter(Action == 0) %>%
  select(Next_H_Center, Probability) %>%
  pivot_wider(names_from = Next_H_Center, values_from = Probability, names_prefix = "H'=")
cat("\nProbabilities for Action 0 (Lifestyle):\n")
print(table_a0)

table_a1 <- prob_single_start_df %>%
  filter(Action == 1) %>%
  select(Next_H_Center, Probability) %>%
  pivot_wider(names_from = Next_H_Center, values_from = Probability, names_prefix = "H'=")
cat("\nProbabilities for Action 1 (Medication):\n")
print(table_a1)