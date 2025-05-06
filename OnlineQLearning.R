library(ggplot2)
library(dplyr)
library(patchwork)

H_bins <- seq(4, 12, 1)
B_bins <- seq(15, 45, 2)
Z_bins <- seq(0, 10, 2)
state_map_df <- expand.grid(H = H_bins, B = B_bins, Z = Z_bins) %>%
  mutate(s_idx = row_number())
n_states <- nrow(state_map_df)
actions <- c(0,1)
tab_key <- with(state_map_df, paste(H, B, Z, sep="_"))
idx_lookup <- setNames(state_map_df$s_idx, tab_key)

T_horizon     <- 30
alpha         <- 0.1
gamma         <- 0.99
epsilon0      <- 1.0
epsilon_min   <- 0.05
N_episodes    <- 3000000
epsilon_decay <- (epsilon_min/epsilon0)^(1/N_episodes)
C_med         <- 1.0

discretise <- function(x,bins) bins[ which.min(abs(bins - x)) ]

Q_list <- lapply(1:(T_horizon+1), function(i) matrix(0, n_states, length(actions)))

set.seed(42)
epsilon <- epsilon0
pb <- txtProgressBar(min=0, max=N_episodes, style=3)
for(ep in seq_len(N_episodes)){
  M <- rnorm(1,0,1)
  H_d <- sample(H_bins,1)
  B_d <- sample(B_bins,1)
  Z_d <- sample(Z_bins,1)
  H_cont <- H_d; B_cont <- B_d
  for(t in 1:T_horizon){
    key <- paste(H_d, B_d, Z_d, sep="_")
    s_idx <- idx_lookup[key]
    if(runif(1)<epsilon) a_t <- sample(actions,1) else {
      qv <- Q_list[[t]][s_idx, ]; a_t <- actions[which.max(qv)]
    }
    Z_next      <- sample(Z_bins,1)
    H_next_cont <- H_cont + 0.05 + (-0.20)*a_t + (-0.10)*M + rnorm(1,0,0.1)
    B_next_cont <- B_cont + 0.02 + (0.10)*a_t + (-0.05)*M + rnorm(1,0,0.2)
    H_next_d    <- discretise(H_next_cont, H_bins)
    B_next_d    <- discretise(B_next_cont, B_bins)
    r_t <- -H_next_cont - ifelse(a_t==1, C_med, 0)
    next_key   <- paste(H_next_d,B_next_d,Z_next,sep="_")
    s_next_idx <- idx_lookup[next_key]
    y_t <- r_t + gamma * max(Q_list[[t+1]][s_next_idx, ])
    old <- Q_list[[t]][s_idx, a_t+1]
    Q_list[[t]][s_idx, a_t+1] <- old + alpha*(y_t-old)
    H_d <- H_next_d; B_d <- B_next_d; Z_d <- Z_next
    H_cont <- H_next_cont; B_cont <- B_next_cont
  }
  epsilon <- max(epsilon*epsilon_decay, epsilon_min)
  if(ep %% 1000==0) setTxtProgressBar(pb,ep)
}
close(pb)

pi_matrix_tab <- sapply(Q_list[1:T_horizon], function(Q) ifelse(Q[,2]>Q[,1],1,0))

library(ggplot2)
library(dplyr)
library(patchwork)

H_bins <- 4:12
B_bins <- seq(15, 45, 2)
Z_bins <- seq(0, 10, 2)
state_map_df <- expand.grid(H = H_bins, B = B_bins, Z = Z_bins) %>%
  mutate(s_idx = row_number())

plot_policy_slice_online <- function(stage_t, z_idx, pi_matrix, state_map_df,
                                     h_bins, b_bins, z_bins) {
  if (stage_t < 1 || stage_t > ncol(pi_matrix)) stop("Invalid stage_t")
  if (z_idx   < 1 || z_idx   > length(z_bins)) stop("Invalid z_idx")
  
  policy_df <- data.frame(s_idx = seq_len(nrow(pi_matrix)),
                          Action = pi_matrix[, stage_t]) %>%
    left_join(state_map_df, by = "s_idx") %>%
    filter(Z == z_bins[z_idx]) %>%
    mutate(
      H_center      = H,
      B_center      = B,
      Action_Factor = factor(Action,
                             levels = c(0,1),
                             labels = c("Lifestyle","Medication"))
    )
  
  req_levs <- c("Lifestyle","Medication")
  pres_levs <- unique(policy_df$Action_Factor)
  dummy_rows <- list()
  dh <- min(h_bins) - (max(h_bins)-min(h_bins))
  db <- min(b_bins) - (max(b_bins)-min(b_bins))
  if (!"Lifestyle"  %in% pres_levs)
    dummy_rows[["Lifestyle" ]] <-
    data.frame(H_center=dh, B_center=db,
               Action_Factor=factor("Lifestyle", req_levs))
  if (!"Medication" %in% pres_levs)
    dummy_rows[["Medication"]] <-
    data.frame(H_center=dh, B_center=db,
               Action_Factor=factor("Medication", req_levs))
  
  plot_data <- if (length(dummy_rows))
    bind_rows(policy_df, bind_rows(dummy_rows))
  else policy_df
  
  ggplot(plot_data,
         aes(x = B_center, y = H_center, fill = Action_Factor)) +
    geom_tile(color = "grey90", linewidth = 0.1) +
    scale_fill_manual(
      name = "Optimal Action",
      values = c("Lifestyle"  = "lightblue",
                 "Medication" = "salmon"),
      drop   = FALSE
    ) +
    coord_cartesian(xlim = range(b_bins), ylim = range(h_bins), expand = FALSE) +
    labs(
      title = sprintf("Stage %d, Z = %.1f", stage_t, z_bins[z_idx]),
      x     = "BMI",
      y     = "HbA1c"
    ) +
    theme_minimal(base_size = 10) +
    theme(
      plot.title      = element_text(hjust = 0.5, size = 12),
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y     = element_text(size = 8)
    )
}

stages_to_plot <- c(1, 15, 30)
z_values       <- c(0, 4, 10)
z_indices      <- match(z_values, Z_bins)

each_plot <- list()
for (t in stages_to_plot) {
  for (zi in z_indices) {
    if (!is.na(zi)) {
      each_plot[[paste0("S",t,"_Z",z_values[which(z_indices==zi)])]] <-
        plot_policy_slice_online(
          stage_t      = t,
          z_idx        = zi,
          pi_matrix    = pi_matrix_tab,
          state_map_df = state_map_df,
          h_bins       = H_bins,
          b_bins       = B_bins,
          z_bins       = Z_bins
        )
    }
  }
}

combined_plot <- wrap_plots(
  each_plot,
  ncol    = length(z_values),
  widths  = rep(1.5, length(z_values)),
  heights = rep(1, length(stages_to_plot)),
  guides  = "collect"
) & theme(legend.position = "bottom")

final_plot <- combined_plot +
  plot_annotation(title = "Online Tabular Q-Learning: Optimal Policy at Stages 1,15,30 × Z=0,4,10")

print(final_plot)

t0      <- 29
alpha   <- 0.1
gamma   <- 0.99

s29_cont <- c(h = 4, b = 8, z = 3)
s30_cont <- c(h = 3, b = 7, z = 2)

H29_d <- discretise(s29_cont["h"], H_bins)
B29_d <- discretise(s29_cont["b"], B_bins)
Z29_d <- discretise(s29_cont["z"], Z_bins)

H30_d <- discretise(s30_cont["h"], H_bins)
B30_d <- discretise(s30_cont["b"], B_bins)
Z30_d <- discretise(s30_cont["z"], Z_bins)

idx29 <- idx_lookup[paste(H29_d, B29_d, Z29_d, sep = "_")]
idx30 <- idx_lookup[paste(H30_d, B30_d, Z30_d, sep = "_")]

a29    <- 1
Q29_old<- Q_list[[t0]][idx29, a29+1]
Q30_0  <- Q_list[[t0+1]][idx30, 1]
Q30_1  <- Q_list[[t0+1]][idx30, 2]
maxQ30 <- max(Q30_0, Q30_1)

r30    <- 0.0
Y29    <- r30 + gamma * maxQ30
delta29<- Y29 - Q29_old
Q29_new<- Q29_old + alpha * delta29

cat("Discrete s29 = (", H29_d, ",", B29_d, ",", Z29_d, ")\n")
cat("Q29_old =", sprintf("%.4f\n", Q29_old))
cat("Q30(s',0) =", sprintf("%.4f\n", Q30_0))
cat("Q30(s',1) =", sprintf("%.4f\n", Q30_1), "\n")
cat("max Q30   =", sprintf("%.4f\n", maxQ30))
cat("Y29       =", sprintf("%.4f\n", Y29))
cat("delta29   =", sprintf("%.4f\n", delta29))
cat("Q29_new   =", sprintf("%.4f\n", Q29_new), "\n")

Q_list[[t0]][idx29, a29+1] <- Q29_new
