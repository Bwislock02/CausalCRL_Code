# Causal Reinforcement Learning in Sequential Decision-Making for Healthcare

**Author:** Benjamin Wislocki
**Project:** Causal Reinforcement Learning in Sequential Decision-Making (May 2025, Durham University)

## Project Overview

This repository contains R scripts developed as part of a dissertation focusing on the application of Reinforcement Learning (RL) and Causal Inference to optimise sequential decision-making, particularly for Dynamic Treatment Regimes (DTRs) in healthcare. The primary case study revolves around managing Type 2 Diabetes, where the goal is to learn optimal strategies that adapt to patient states (e.g., HbA1c, BMI, motivation) over time.

The project explores different RL methodologies:
1.  **Model-Based RL:** Using a known or estimated model of the environment (Markov Decision Process).
2.  **Offline (Batch) RL:** Learning from a fixed dataset of previously collected experiences, addressing challenges like confounding using techniques such as inverse probability weighting.
3.  **Online RL:** Learning through direct interaction with an environment (or a simulation of it) relying on the experimental 'NUC' condition.

The scripts herein provide implementations for these approaches within a simulated healthcare context.

## Core Problem Addressed

Standard RL methods can falter when learning from observational data due to confounding variables that affect both the treatment (action) received and the outcome. This project integrates causal reasoning to mitigate such biases, aiming for more reliable policy evaluation and optimisation in real-world healthcare settings.

## Scripts in this Repository

This repository includes three main R scripts, each demonstrating a different RL algorithm:

### 1. `MDPValueIteration.R` - Model-Based Reinforcement Learning

* **Purpose:** Implements a model-based RL approach where the environment dynamics are explicitly defined and used to find an optimal policy.
* **Method:**
    * Defines a Markov Decision Process (MDP) with states (discretised HbA1c, BMI, and a proxy for patient motivation Z), actions (lifestyle intervention vs. medication), and a reward function (penalising high HbA1c and accounting for medication costs).
    * Simulates patient state transitions based on defined dynamics and action effects.
    * Estimates the transition probability matrix $P(s' | s, a)$ via Monte Carlo simulation.
    * Applies Value Iteration to compute the optimal state-value function ($V^*$) and derive the optimal policy ($\pi^*$) over a finite horizon.
* **Key Features:**
    * Detailed simulation of patient physiological changes and motivation.
    * Discretisation of continuous state variables.
    * Explicit estimation and use of the transition model.
    * Visualisation of the learned policy as heatmaps for different patient states and stages.
    * Analysis of transition probabilities for specific states.
* **Relevance:** Serves as a foundational example of solving an MDP when the model is known or can be accurately simulated. This is a benchmark against which model-free methods can be compared.

### 2. `OfflineWeightedQLearning.R` - Offline Reinforcement Learning

* **Purpose:** Implements an offline (batch) RL algorithm, specifically Weighted Q-Learning, to learn an optimal policy from a fixed dataset of simulated patient trajectories.
* **Method:**
    * Simulates a dataset of patient trajectories under a (potentially suboptimal or confounded) behavioural policy. This data mimics observational data.
    * Applies Weighted Q-Learning, which uses Inverse Probability Weighting (IPW) to account for confounding bias. Propensity scores (probability of taking an action given the state) are estimated and used to re-weight the contributions of observed transitions.
    * Learns Q-functions for each stage by fitting weighted linear models.
* **Key Features:**
    * Data simulation function to generate longitudinal data with a known underlying structure.
    * Implementation of Weighted Q-Learning for DTR estimation.
    * Propensity score estimation for weighting.
    * Standardisation of state variables.
    * Visualisation of the learned policy.
* **Relevance:** Demonstrates how to learn effective policies from observational data, a crucial capability in healthcare where experimental data collection can be costly or unethical. Addresses confounding, a key theme of the dissertation.

### 3. `OnlineQLearning.R` - Online Reinforcement Learning

* **Purpose:** Implements an online RL algorithm, specifically tabular Q-Learning, to learn an optimal policy through simulated direct interaction with the environment.
* **Method:**
    * Defines a discretised state space (HbA1c, BMI, Z) and actions.
    * Uses Q-Learning with an epsilon-greedy exploration strategy to learn the optimal action-value function ($Q^*(s,a)$) for each state and action at each stage.
    * Updates Q-values iteratively based on experienced transitions and rewards during simulated episodes.
    * The simulation model for patient transitions is similar to that in `MDPValueIteration.R` but is interacted with step-by-step.
* **Key Features:**
    * Tabular Q-Learning implementation for a finite horizon.
    * Epsilon-greedy exploration strategy with decay.
    * Online learning through simulated patient episodes.
    * Discretisation of state variables for tabular representation.
    * Visualisation of the learned policy.
* **Relevance:** Illustrates a model-free approach where the agent learns by trial and error, suitable when an explicit model of the environment is not available but interaction is possible (e.g., through a high-fidelity simulator or, cautiously, in real-world adaptive interventions).

## Technical Details

* **Language:** R
* **Key R Libraries Used:**
    * `dplyr` (data manipulation)
    * `tidyr` (data tidying)
    * `ggplot2` (visualisation)
    * `patchwork` (combining plots)
    * `pbapply` (progress bars for apply functions - in `MDPValueIteration.R`)

## How to Use

1.  **Prerequisites:** Ensure you have R installed, along with the R libraries listed above. You can install them using `install.packages(c("dplyr", "tidyr", "ggplot2", "patchwork", "pbapply"))`.
2.  **Running the Scripts:**
    * Open R or an R IDE (like RStudio).
    * Source the desired script, e.g., `source("MDPValueIteration.R")`.
    * The scripts are generally self-contained and will run the simulations, learn the policies, and produce output (console messages, plots).
3.  **Outputs:**
    * **Console Output:** Progress updates, summaries of learned policies/values for example states.
    * **Plots:** Heatmaps of optimal policies, typically showing recommended actions (e.g., "Lifestyle" vs. "Medication") based on patient HbA1c, BMI, and motivation (Z) at different decision stages.
    * **Saved Data (for `MDPValueIteration.R`):** The estimated transition probability matrix (`transition_prob_sims<N>.rds`) is saved to avoid re-computation on subsequent runs.

## Simulation Environment (Conceptual)

The scripts simulate a simplified model of Type 2 Diabetes progression:
* **States:** Patient's HbA1c, BMI, and a latent "motivation" factor (Z).
* **Actions:** Typically binary – e.g., (0) lifestyle advice or (1) medication.
* **Transitions:** How HbA1c and BMI change depends on the current state, the action taken, adherence (influenced by motivation and action), and some stochasticity. Motivation (Z) itself can evolve.
* **Rewards:** Primarily based on achieving lower HbA1c levels, with costs associated with certain actions (e.g., medication).

## Connection to Dissertation

These scripts provide the practical implementation and empirical results for the methodologies discussed in the dissertation "Causal Reinforcement Learning in Sequential Decision-Making." They illustrate:
* The foundational principles of MDPs and value iteration (`MDPValueIteration.R`).
* The challenges of learning from observational data and how causal inference techniques (like IPW in `OfflineWeightedQLearning.R`) can be integrated into RL.
* The mechanics of online, model-free learning (`OnlineQLearning.R`).

The results from these simulations (e.g., learned policies, estimated value functions) are used to support the theoretical discussions and methodological contributions of the dissertation.
