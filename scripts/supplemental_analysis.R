rm(list = ls())

# setup
library(tidyverse)
library(ggbeeswarm)
library(here)
library(effsize)
library(lme4)
source(here("scripts", "utils.R"))

# set plot theme
plot_theme <- theme_light() + 
  theme(panel.grid = element_blank(),
        legend.position = 'none')
theme_set(plot_theme)
cond_cols <- c('salmon2','dodgerblue2')

data_dir <- "clean_data"
raw_data_dir <- "raw_data"

# ----------------
# AISLE REPEATS (RECON TEST)
# ----------------

# EXPERIMENT 1

recon <- read.csv(here(data_dir, "exp1_recon_clean.csv"), stringsAsFactors = F)
enc <- read.csv(here(raw_data_dir, "exp1_enc.csv"), stringsAsFactors = F)

# get duplicated positions in each event
pos_duplicate <- enc %>%
  filter(startsWith(display, "motor_trial")) %>%
  group_by(id, list, event) %>%
  summarise(aisle_dup1 = get_mode(spatial_pos)[1],
            aisle_dup2 = get_mode(spatial_pos)[2])

recon_dup <- inner_join(recon, pos_duplicate, by = c('id','list','event')) %>%
  mutate(aisle_dup = ifelse(spatial_pos == aisle_dup1 | spatial_pos == aisle_dup2, 1, 0))

# acc by condition and aisle duplicate
recon_acc <- recon_dup %>%
  group_by(id, condition, aisle_dup) %>%
  summarise(mean_acc = mean(accuracy)) %>%
  ungroup()
recon_acc_group <- recon_acc %>%
  group_by(condition, aisle_dup) %>%
  summarise(group_mean_acc = mean(mean_acc), group_sem_acc = sem(mean_acc)) %>%
  ungroup()

# stats - ttest
t.test(mean_acc ~ condition, data = filter(recon_acc, aisle_dup == 1), paired = T)
effsize::cohen.d(mean_acc ~ condition | Subject(id), data = filter(recon_acc, aisle_dup == 1), paired = T)
t.test(mean_acc ~ condition, data = filter(recon_acc, aisle_dup == 0), paired = T)
effsize::cohen.d(mean_acc ~ condition | Subject(id), data = filter(recon_acc, aisle_dup == 0), paired = T)

# stats - regression model
recon_model_data <- mutate(recon_dup, condition_e = ifelse(condition == 'pred', 0.5, -0.5), # effect coded
                           dup_e = ifelse(aisle_dup == 0, -0.5, 0.5), # effect coded
                           seq_pos_z = scale(seq_pos, center = T, scale = F)) # center pos = 0
recon_model <- glmer(accuracy ~ condition_e + seq_pos_z + dup_e + (condition_e + seq_pos_z + dup_e || id),
                     data = recon_model_data, family = 'binomial')
summary(recon_model)
confint(recon_model, method = "Wald")

# EXPERIMENT 2

recon <- read.csv(here(data_dir, "exp2_recon_clean.csv"), stringsAsFactors = F)
enc <- read.csv(here(raw_data_dir, "exp2_enc.csv"), stringsAsFactors = F)

# get duplicated positions in each event
pos_duplicate <- enc %>%
  filter(startsWith(display, "motor_trial")) %>%
  group_by(id, list, event) %>%
  summarise(aisle_dup1 = get_mode(spatial_pos)[1],
            aisle_dup2 = get_mode(spatial_pos)[2])

recon_dup <- inner_join(recon, pos_duplicate, by = c('id','list','event')) %>%
  mutate(aisle_dup = ifelse(spatial_pos == aisle_dup1 | spatial_pos == aisle_dup2, 1, 0))

# acc by condition and aisle duplicate
recon_acc <- recon_dup %>%
  group_by(id, condition, aisle_dup) %>%
  summarise(mean_acc = mean(accuracy)) %>%
  ungroup()
recon_acc_group <- recon_acc %>%
  group_by(condition, aisle_dup) %>%
  summarise(group_mean_acc = mean(mean_acc), group_sem_acc = sem(mean_acc)) %>%
  ungroup()

# stats - ttest
t.test(mean_acc ~ condition, data = filter(recon_acc, aisle_dup == 1), paired = T)
effsize::cohen.d(mean_acc ~ condition | Subject(id), data = filter(recon_acc, aisle_dup == 1), paired = T)
t.test(mean_acc ~ condition, data = filter(recon_acc, aisle_dup == 0), paired = T)
effsize::cohen.d(mean_acc ~ condition | Subject(id), data = filter(recon_acc, aisle_dup == 0), paired = T)

# stats - regression model
recon_model_data <- mutate(recon_dup, condition_e = ifelse(condition == 'pred', 0.5, -0.5), # effect coded
                           dup_e = ifelse(aisle_dup == 0, -0.5, 0.5), # effect coded
                           seq_pos_z = scale(seq_pos, center = T, scale = F)) # center pos = 0
recon_model <- glmer(accuracy ~ condition_e + seq_pos_z + dup_e + (condition_e + seq_pos_z + dup_e || id),
                     data = recon_model_data, family = 'binomial')
summary(recon_model)
confint(recon_model, method = "Wald")

# EXPERIMENT 3
recon <- read.csv(here(data_dir, "exp3_recon_clean.csv"), stringsAsFactors = F)
enc <- read.csv(here(raw_data_dir, "exp3_enc.csv"), stringsAsFactors = F)

# get duplicated positions in each event
pos_duplicate <- enc %>%
  filter(startsWith(display, "motor_trial")) %>%
  group_by(id, list, event) %>%
  summarise(aisle_dup1 = get_mode(spatial_pos)[1],
            aisle_dup2 = get_mode(spatial_pos)[2])

recon_dup <- inner_join(recon, pos_duplicate, by = c('id','list','event')) %>%
  mutate(aisle_dup = ifelse(spatial_pos == aisle_dup1 | spatial_pos == aisle_dup2, 1, 0))

# acc by condition and aisle duplicate
recon_acc <- recon_dup %>%
  group_by(id, condition, aisle_dup) %>%
  summarise(mean_acc = mean(accuracy)) %>%
  ungroup()
recon_acc_group <- recon_acc %>%
  group_by(condition, aisle_dup) %>%
  summarise(group_mean_acc = mean(mean_acc), group_sem_acc = sem(mean_acc)) %>%
  ungroup()

# stats - ttest
t.test(mean_acc ~ condition, data = filter(recon_acc, aisle_dup == 1), paired = T)
effsize::cohen.d(mean_acc ~ condition | Subject(id), data = filter(recon_acc, aisle_dup == 1), paired = T)
t.test(mean_acc ~ condition, data = filter(recon_acc, aisle_dup == 0), paired = T)
effsize::cohen.d(mean_acc ~ condition | Subject(id), data = filter(recon_acc, aisle_dup == 0), paired = T)

# stats - regression model
recon_model_data <- mutate(recon_dup, condition_e = ifelse(condition == 'pred', 0.5, -0.5), # effect coded
                           dup_e = ifelse(aisle_dup == 0, -0.5, 0.5), # effect coded
                           seq_pos_z = scale(seq_pos, center = T, scale = F)) # center pos = 0
recon_model <- glmer(accuracy ~ condition_e + seq_pos_z + dup_e + (condition_e + seq_pos_z + dup_e || id),
                     data = recon_model_data, family = 'binomial')
summary(recon_model)
confint(recon_model, method = "Wald") 

# ----------------
# NEXT THING
# ----------------
