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

# ----------------
# PRETRAINING
# ----------------

## 1: pretraining study RTs for pred v. rand sequences

# load data
pretrain_study <- read.csv(here(data_dir, "exp1_pretrain_study_clean.csv"), stringsAsFactors = F)

# remove outlier RTs
all_rts <- filter(pretrain_study, !is.na(response))$rt
sd_cutoff <- mean(all_rts) + (3 * sd(all_rts))
pretrain_study <- filter(pretrain_study, rt <= sd_cutoff)

# mean RTs by condition
pretrain_study_rt <- pretrain_study %>%
  filter(block == 3) %>%
  group_by(id, condition) %>%
  summarise(mean_rt = mean(rt)) %>%
  ungroup()
pretrain_study_rt_group <- pretrain_study_rt %>%
  group_by(condition) %>%
  summarise(group_mean_rt = mean(mean_rt), group_sem_rt = sem(mean_rt)) %>%
  ungroup()

# plot
ggplot(pretrain_study_rt, aes(x = condition, y = mean_rt, fill = condition)) +
  geom_bar(data = pretrain_study_rt_group, aes(x = condition, y = group_mean_rt),
           stat = 'identity', width = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1.2, color = 'white') +
  geom_errorbar(data = pretrain_study_rt_group, width = 0.1,
                aes(x = condition, y = group_mean_rt,
                    ymin = group_mean_rt - group_sem_rt, ymax = group_mean_rt + group_sem_rt)) +
  scale_fill_manual(values = cond_cols) +
  labs(x = 'condition', y = 'RT (ms)')

# stats
wilcox.test(mean_rt ~ condition, data = pretrain_study_rt, paired = T)
effsize::cliff.delta(mean_rt ~ condition, data = pretrain_study_rt, paired = T)

## 2: pretraining test accuracy for pred sequences

pretrain_test <- read.csv(here(data_dir, "exp1_pretrain_test_clean.csv"), stringsAsFactors = F)

# mean acc by condition and test repetition
pretrain_test_acc <- pretrain_test %>%
  group_by(id, test_rep) %>%
  summarise(mean_acc = mean(accuracy)) %>%
  ungroup()
pretrain_test_acc_group <- pretrain_test_acc %>%
  group_by(test_rep) %>%
  summarise(group_mean_acc = mean(mean_acc),
            group_sem_acc = sem(mean_acc)) %>%
  ungroup()

# plot
ggplot(pretrain_test_acc, aes(x = test_rep, y = mean_acc)) +
  geom_jitter(width = 0.2, height = 0.02, fill = 'dodgerblue1', color = 'white',
              alpha = 0.4, shape = 21, size = 2.3) +
  geom_line(data = pretrain_test_acc_group, aes(x = test_rep, y = group_mean_acc),
            size = 1.5, color = 'dodgerblue4') +
  geom_errorbar(data = pretrain_test_acc_group,
                aes(x = test_rep, y = group_mean_acc,
                    ymin = group_mean_acc - group_sem_acc, ymax = group_mean_acc + group_sem_acc),
                size = 1.5, width = 0, color = 'dodgerblue4') +
  scale_x_continuous(breaks = 1:6) +
  labs(x = 'sequence test repetition', y = 'accuracy')

# show get descriptive stats (how many/what proportion of participants reach each accuracy level)
pretrain_test_acc %>%
  group_by(test_rep) %>%
  summarise(n_80 = sum(mean_acc > 0.8),
            n_100 = sum(mean_acc == 1),
            prop_100 = sum(mean_acc == 1) / length(mean_acc),
            prop_80 = sum(mean_acc > 0.8) / length(mean_acc))


# ----------------
# ENCODING
# ----------------

## 1: encoding accuracy for pred v. rand events

enc <- read.csv(here(data_dir, "exp1_enc_clean.csv"), stringsAsFactors = F)

# mean acc by condition
enc_resp_acc <- enc %>%
  filter(startsWith(display, "motor_trial")) %>%
  group_by(id, condition) %>%
  summarise(mean_acc = mean(accuracy)) %>%
  ungroup()
enc_resp_acc_group <- enc_resp_acc %>%
  group_by(condition) %>%
  summarise(group_mean_acc = mean(mean_acc), group_sem_acc = sem(mean_acc)) %>%
  ungroup()

# plot
ggplot(enc_resp_acc, aes(x = condition, y = mean_acc, fill = condition)) +
  geom_bar(data = enc_resp_acc_group, aes(x = condition, y = group_mean_acc),
           stat = 'identity', width = 0.8) +
  geom_beeswarm(aes(x = condition, y = mean_acc),
                dodge.width = 1, groupOnX = T, size = 3, pch = 21, color = 'white') +
  geom_errorbar(data = enc_resp_acc_group, width = 0.1,
                aes(x = condition, y = group_mean_acc,
                    ymin = group_mean_acc - group_sem_acc, ymax = group_mean_acc + group_sem_acc)) +
  scale_fill_manual(values = cond_cols) +
  labs(x = 'condition', y = 'aisle response accuracy')
  
# stats
t.test(mean_acc ~ condition, data = enc_resp_acc, paired = T)
effsize::cohen.d(mean_acc ~ condition | Subject(id), data = enc_resp_acc, paired = T)

# mean acc by condition and block/list number
enc_resp_acc <- enc %>%
  filter(startsWith(display, "motor_trial")) %>%
  group_by(id, condition, list) %>%
  summarise(mean_acc = mean(accuracy)) %>%
  ungroup()
enc_resp_acc_group <- enc_resp_acc %>%
  group_by(condition, list) %>%
  summarise(group_mean_acc = mean(mean_acc), group_sem_acc = sem(mean_acc)) %>%
  ungroup()

# plot
ggplot(enc_resp_acc_group, aes(x = list, y = group_mean_acc, color = condition)) +
  geom_jitter(data = enc_resp_acc, aes(x = list, y = mean_acc, fill = condition),
              width = 0.1, height = 0.02, color = 'white', alpha = 0.4, shape = 21, size = 2.3) +
  geom_line(size = 1.5) +
  geom_errorbar(width = 0, size = 1,
                aes(x = list,  ymin = group_mean_acc - group_sem_acc,
                    ymax = group_mean_acc + group_sem_acc)) +
  scale_color_manual(values = cond_cols) + scale_fill_manual(values = cond_cols) +
  labs(x = 'block number', y = 'aisle response accuracy')

# stats
ttest_multiple_fdr(enc_resp_acc, equation = 'mean_acc ~ condition', grouping = 'list')

# 2. encoding RTs for pred v. rand events

# remove outlier RTs
all_rts <- enc$rt
sd_cutoff <- mean(all_rts) + (3 * sd(all_rts))
enc <- filter(enc, rt <= sd_cutoff)

# mean rt by condition
enc_resp_rt <- enc %>%
  filter(startsWith(display, "motor_trial")) %>%
  group_by(id, condition) %>%
  summarise(mean_rt = mean(rt)) %>%
  ungroup()

# stats
wilcox.test(mean_rt ~ condition, data = enc_resp_rt, paired = T)

# ----------------
# ORDER RECONSTRUCTION
# ----------------

recon <- read.csv(here(data_dir, "exp1_recon_clean.csv"), stringsAsFactors = F)

# 1. order recon accuracy for pred v. rand events — ordinal accuracy

# order memory by condition
recon_acc <- recon %>%
  group_by(id, condition) %>%
  summarise(mean_acc = mean(accuracy)) %>%
  ungroup()
recon_acc_group <- recon_acc %>%
  group_by(condition) %>%
  summarise(group_mean_acc = mean(mean_acc), group_sem_acc = sem(mean_acc)) %>%
  ungroup()

# plot
ggplot(recon_acc, aes(x = condition, y = mean_acc, fill = condition)) +
  geom_bar(data = recon_acc_group, aes(x = condition, y = group_mean_acc),
           stat = 'identity', width = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1, color = 'white') +
  geom_errorbar(data = recon_acc_group, width = 0.1,
                aes(x = condition, y = group_mean_acc,
                    ymin = group_mean_acc - group_sem_acc, ymax = group_mean_acc + group_sem_acc)) +
  geom_hline(yintercept = 1/6, color = 'grey20', linetype = 'dashed') +
  scale_fill_manual(values = cond_cols) +
  labs(x = 'condition', y = 'ordinal accuracy')

t.test(mean_acc ~ condition, data = recon_acc, paired = T)
effsize::cohen.d(mean_acc ~ condition | Subject(id), data = recon_acc, paired = T)

# 2. order recon accuracy for pred v. rand events — levenshtein distance

# get dataframe with this measure (implemented in utils.R script for readability)
recon_leven_data <- suppressMessages(get_leven_dist(recon))
# without message suppression, this will print out a bunch of warnings, because participants
# don't always select all 6 of the items during each recon trial

# levenshtein distance by condition
recon_leven <- recon_leven_data %>%
  group_by(id, condition) %>%
  summarise(mean_leven_dist = mean(leven_dist)) %>%
  ungroup()
recon_leven_group <- recon_leven %>%
  group_by(condition) %>%
  summarise(group_mean_leven_dist = mean(mean_leven_dist), group_sem_leven_dist = sem(mean_leven_dist)) %>%
  ungroup()

# stats
t.test(mean_leven_dist ~ condition, data = recon_leven, paired = T)
effsize::cohen.d(mean_leven_dist ~ condition | Subject(id), data = recon_leven, paired = T)

# 2. proportion of items remembered in correct order within pred v. rand events

recon_event <- recon %>%
  group_by(id, list, event) %>%
  summarise(n_correct = sum(accuracy), condition = get_mode(condition)) %>%
  group_by(id, condition) %>%
  summarise(prop_6 = sum(n_correct == 6) / length(n_correct),
            prop_5 = sum(n_correct == 5) / length(n_correct),
            prop_4 = sum(n_correct == 4) / length(n_correct),
            prop_3 = sum(n_correct == 3) / length(n_correct),
            prop_2 = sum(n_correct == 2) / length(n_correct),
            prop_1 = sum(n_correct == 1) / length(n_correct),
            prop_0 = sum(n_correct == 0) / length(n_correct)) %>%
  pivot_longer(cols = starts_with('prop_'), names_to = 'num_correct', values_to = 'proportion') %>%
  mutate(num_correct = as.numeric(gsub('prop_', '', num_correct)),
         condition = as.factor(condition)) %>%
  ungroup()
recon_event_group <- recon_event %>%
  group_by(condition, num_correct) %>%
  summarise(group_mean_prop = mean(proportion), group_sem_prop = sem(proportion)) %>%
  ungroup()

# plot
ggplot(recon_event_group, aes(x = num_correct, y = group_mean_prop, fill = condition)) +
  geom_bar(stat = 'identity', position = position_dodge(0.9)) +
  geom_point(data = recon_event, aes(x = num_correct, y = proportion, color = condition),
             position = position_jitterdodge(jitter.width = 0.3, jitter.height = 0.01, dodge.width = 0.9),
             pch = 21, color = 'white', size = 2) +
  geom_errorbar(aes(ymin = group_mean_prop - group_sem_prop,
                    ymax = group_mean_prop + group_sem_prop),
                position = position_dodge(0.9), width = 0) +
  labs(x = 'number of items selected in correct order (per event)', y = 'proportion of events') +
  scale_fill_manual(values = cond_cols) + scale_color_manual(values = cond_cols)

# stats
ttest_multiple_fdr(recon_event, equation = "proportion ~ condition", grouping = "num_correct")

# 3. order recon for pred v. rand events as a function of sequence position

# order memory by condition and seq pos
recon_acc <- recon %>%
  group_by(id, condition, seq_pos) %>%
  summarise(mean_acc = mean(accuracy)) %>%
  ungroup()
recon_acc_group <- recon_acc %>%
  group_by(condition, seq_pos) %>%
  summarise(group_mean_acc = mean(mean_acc), group_sem_acc = sem(mean_acc)) %>%
  ungroup()

# plot
ggplot(recon_acc_group, aes(x = seq_pos, y = group_mean_acc, color = condition)) +
  geom_hline(yintercept = 1/6, color = 'grey20', linetype = 'dashed') +
  geom_jitter(data = recon_acc, aes(x = seq_pos, y = mean_acc, fill = condition),
              width = 0.1, height = 0.02, color = 'white', alpha = 0.4, shape = 21, size = 3) +
  geom_line(size = 1.5) +
  geom_errorbar(aes(x = seq_pos,
                    ymin = group_mean_acc - group_sem_acc, ymax = group_mean_acc + group_sem_acc),
                width = 0, size = 1) +
  scale_color_manual(values = cond_cols) + scale_fill_manual(values = cond_cols) +
  labs(x = 'sequence position', y = 'ordinal accuracy')

# stats - ttest
ttest_multiple_fdr(recon_acc, equation = "mean_acc ~ condition", grouping = "seq_pos")

# stats - logistic regression model
recon_model_data <- mutate(recon, condition_e = ifelse(condition == 'pred', 0.5, -0.5), # effect coded
                           seq_pos_z = scale(seq_pos, center = T, scale = F)) # center pos = 0
recon_model <- glmer(accuracy ~ condition_e * seq_pos_z + (condition_d + seq_pos_z || id),
                     data = recon_model_data, family = 'binomial')
summary(recon_model)
confint(recon_model) # takes awhile to run


# ----------------
# SPATIAL MEMORY
# ----------------

spatial <- read.csv(here(data_dir, "exp1_spatial_clean.csv"), stringsAsFactors = F)

# 1. spatial memory for pred v. rand events

# mean acc by condition
spatial_acc <- spatial %>%
  group_by(id, condition) %>%
  summarise(mean_acc = mean(accuracy)) %>%
  ungroup()
spatial_acc_group <- spatial_acc %>%
  group_by(condition) %>%
  summarise(group_mean_acc = mean(mean_acc), group_sem_acc = sem(mean_acc)) %>%
  ungroup()

# plot
ggplot(spatial_acc, aes(x = condition, y = mean_acc, fill = condition)) +
  geom_bar(data = spatial_acc_group, aes(x = condition, y = group_mean_acc),
           stat = 'identity', width = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1, color = 'white') +
  geom_errorbar(width = 0.1, data = spatial_acc_group,
                aes(x = condition, y = group_mean_acc,
                    ymin = group_mean_acc - group_sem_acc, ymax = group_mean_acc + group_sem_acc)) +
  geom_hline(yintercept = 1/4, color = 'grey20', linetype = 'dashed') +
  scale_fill_manual(values = cond_cols) +
  labs(x = 'condition', y = 'spatial memory accuracy')

# stats
t.test(mean_acc ~ condition, data = spatial_acc, paired = T)
effsize::cohen.d(mean_acc ~ condition | Subject(id), data = spatial_acc, paired = T)

# 2. spatial memory for pred v. rand events as a function of sequence position

# order memory by condition and seq pos
spatial_acc <- spatial %>%
  group_by(id, condition, seq_pos) %>%
  summarise(mean_acc = mean(accuracy)) %>%
  ungroup()
spatial_acc_group <- spatial_acc %>%
  group_by(condition, seq_pos) %>%
  summarise(group_mean_acc = mean(mean_acc), group_sem_acc = sem(mean_acc)) %>%
  ungroup()

# plot
ggplot(spatial_acc_group, aes(x = seq_pos, y = group_mean_acc, color = condition)) +
  geom_hline(yintercept = 1/4, color = 'grey20', linetype = 'dashed') +
  geom_jitter(data = spatial_acc, aes(x = seq_pos, y = mean_acc, fill = condition),
              width = 0.1, height = 0.02, color = 'white', alpha = 0.4, shape = 21, size = 3) +
  geom_line(size = 1.5) +
  geom_errorbar(aes(x = seq_pos,
                    ymin = group_mean_acc - group_sem_acc, ymax = group_mean_acc + group_sem_acc),
                width = 0, size = 1) +
  scale_color_manual(values = cond_cols) + scale_fill_manual(values = cond_cols) +
  labs(x = 'sequence position', y = 'spatial memory accuracy')

# stats - ttest
ttest_multiple_fdr(spatial_acc, equation = "mean_acc ~ condition", grouping = "seq_pos")
