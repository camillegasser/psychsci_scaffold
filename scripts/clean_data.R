rm(list = ls())

# change to project directory
setwd("/Volumes/GoogleDrive/My Drive/grad-school/scaffold_behav/manuscript/psychsci_submission/psychsci_scaffold")

# script setup
library(tidyverse)
library(here)

raw_data_dir <- "raw_data"
clean_data_dir <- "clean_data"

## ==============================
## EXPERIMENT 1
## ==============================

# PRETRAINING
# ----------------

pretrain_data <- read.csv(here(raw_data_dir, "/exp1_pretrain.csv"), stringsAsFactors = F)

# split into study + test trials
pretrain_study <- filter(pretrain_data, startsWith(display, "seq_trial"))
pretrain_test <- filter(pretrain_data, display == "test_trial")

# exclude participants who don't make at least 75% of study responses
pretrain_study_check <- pretrain_study %>%
  group_by(id) %>%
  summarise(elim_trials = sum(is.na(response)),
            elim_prop = elim_trials / length(response)) %>%
  ungroup()
pretrain_exclusions <- pretrain_study_check$id[pretrain_study_check$elim_prop > 0.25]

# exclude study trials with outlier RTs (> 3 SDs above the mean)
study_rts <- filter(pretrain_study, !is.na(response))$rt
sd_cutoff <- mean(study_rts) + (3 * sd(study_rts))

# create clean study data file
pretrain_study_clean <- filter(pretrain_study, !(id %in% pretrain_exclusions), !is.na(response))
write.csv(pretrain_study_clean, file = here(clean_data_dir, "exp1_pretrain_study_clean.csv"), row.names = F)

# create clean test data file
pretrain_test_clean <- filter(pretrain_test, !(id %in% pretrain_exclusions))
write.csv(pretrain_test_clean, file = here(clean_data_dir, "exp1_pretrain_test_clean.csv"), row.names = F)

# ENCODING
# ----------------

enc_data <- read.csv(here(raw_data_dir, "/exp1_enc.csv"), stringsAsFactors = F)

# continue excluding participants that were previously excluded
enc_exclusions <- pretrain_exclusions

# exclude participants who don't make at least 75% of encoding responses
enc_check <- enc_data %>%
  filter(startsWith(display, "motor_trial")) %>%
  group_by(id) %>%
  summarise(elim_trials = sum(is.na(response)),
            elim_prop = elim_trials / length(response)) %>%
  ungroup()
enc_exclusions <- c(enc_exclusions, enc_check$id[enc_check$elim_prop > 0.25])

# exclude participants with poor pretraining performance (< 80% correct on final test)
pretrain_final_perf <- pretrain_test_clean %>%
  filter(test_rep == 6) %>%
  group_by(id) %>%
  summarise(final_acc = mean(accuracy)) %>%
  ungroup()
enc_exclusions <- c(enc_exclusions, pretrain_final_perf$id[pretrain_final_perf$final_acc < 0.8])
enc_exclusions <- unique(enc_exclusions)

# create clean encoding data file
enc_clean <- enc_data %>%
  filter(!(id %in% enc_exclusions), !(is.na(response)),
         !(display %in% c("instructions", "fixation", "event_start")))
write.csv(enc_clean, file = here(clean_data_dir, "exp1_enc_clean.csv"), row.names = F)

# ORDER RECONSTRUCTION
# ----------------

recon_data <- read.csv(here(raw_data_dir, "/exp1_recon.csv"), stringsAsFactors = F)

# continue excluding participants that were previously excluded
recon_exclusions <- enc_exclusions

# exclude participants with poor encoding performance (for predictable/from memory trials)
enc_pred <- enc_clean %>%
  filter(startsWith(display, 'motor_trial'), condition == "pred") %>%
  group_by(id) %>%
  summarise(pred_acc = mean(accuracy)) %>%
  ungroup()
recon_exclusions <- c(recon_exclusions, enc_pred$id[enc_pred$pred_acc < 0.8])

# exclude participants not above chance performance on this memory test
recon_binom <- recon_data %>%
  group_by(id) %>%
  summarise(n_correct = sum(accuracy == 1, na.rm = T),
            n_total = sum(!is.na(response)))

recon_binom[, c("binom_stat", "binom_pval")] <- NA
for (sub in recon_binom$id) {
  bt <- binom.test(filter(recon_binom, id == sub)$n_correct, filter(recon_binom, id == sub)$n_total,
                   p = 1/6, alternative = 'g')
  recon_binom$binom_stat[recon_binom$id == sub] <- bt$statistic
  recon_binom$binom_pval[recon_binom$id == sub] <- bt$p.value
}
recon_exclusions <- c(recon_exclusions, recon_binom$id[recon_binom$binom_pval > 0.05])
recon_exclusions <- unique(recon_exclusions)

# create clean recon file
recon_clean <- recon_data %>%
  filter(!(id %in% recon_exclusions), rt >= 100, !is.na(response))
write.csv(recon_clean, file = here(clean_data_dir, "exp1_recon_clean.csv"), row.names = F)

# SPATIAL MEMORY
# ----------------

spatial_data <- read.csv(here(raw_data_dir, "/exp1_spatial.csv"), stringsAsFactors = F)

# continue excluding participants that were previously excluded
spatial_exclusions <- enc_exclusions

# exclude participants with poor encoding performance (for predictable/from memory trials)
spatial_exclusions <- c(spatial_exclusions, enc_pred$id[enc_pred$pred_acc < 0.8])

# exclude participants not above chance performance on this memory test
spatial_binom <- spatial_data %>%
  group_by(id) %>%
  summarise(n_correct = sum(accuracy == 1, na.rm = T),
            n_total = sum(!is.na(response)))

spatial_binom[, c("binom_stat", "binom_pval")] <- NA
for (sub in spatial_binom$id) {
  bt <- binom.test(filter(spatial_binom, id == sub)$n_correct, filter(spatial_binom, id == sub)$n_total,
                   p = 1/4, alternative = 'g')
  spatial_binom$binom_stat[spatial_binom$id == sub] <- bt$statistic
  spatial_binom$binom_pval[spatial_binom$id == sub] <- bt$p.value
}
spatial_exclusions <- c(spatial_exclusions, spatial_binom$id[spatial_binom$binom_pval > 0.05])
spatial_exclusions <- unique(spatial_exclusions)

# create clean spatial file
spatial_clean <- spatial_data %>%
  filter(!(id %in% spatial_exclusions), rt >= 100, !is.na(response))
write.csv(spatial_clean, file = here(clean_data_dir, "exp1_spatial_clean.csv"), row.names = F)

# REMINDER
# ----------------
remind_data <- read.csv(here(raw_data_dir, "/exp1_reminder.csv"), stringsAsFactors = F)

# continue excluding participants that were previously excluded (prior to memory test analysis)
remind_exclusions <- c(enc_exclusions, enc_pred$id[enc_pred$pred_acc < 0.8])

# create clean reminder file (and add repetition info)
remind_clean <- remind_data %>%
  filter(!(id %in% remind_exclusions)) %>%
  group_by(id) %>%
  mutate(test_rep = rep(c(1,2), each = 6)) %>%
  ungroup()

write.csv(remind_clean, file = here(clean_data_dir, "exp1_remind_clean.csv"), row.names = F)

## ==============================
## EXPERIMENT 2
## ==============================

# PRETRAINING
# ----------------

pretrain_data <- read.csv(here(raw_data_dir, "/exp2_pretrain.csv"), stringsAsFactors = F)

# split into study + test trials
pretrain_study <- filter(pretrain_data, startsWith(display, "seq_trial"))
pretrain_test <- filter(pretrain_data, display == "test_trial")

# exclude participants who don't make at least 75% of study responses
pretrain_study_check <- pretrain_study %>%
  group_by(id) %>%
  summarise(elim_trials = sum(is.na(response)),
            elim_prop = elim_trials / length(response)) %>%
  ungroup()
pretrain_exclusions <- pretrain_study_check$id[pretrain_study_check$elim_prop > 0.25]

# exclude study trials with outlier RTs (> 3 SDs above the mean)
study_rts <- filter(pretrain_study, !is.na(response))$rt
sd_cutoff <- mean(study_rts) + (3 * sd(study_rts))

# create clean study data file
pretrain_study_clean <- filter(pretrain_study, !(id %in% pretrain_exclusions), !is.na(response))
write.csv(pretrain_study_clean, file = here(clean_data_dir, "exp2_pretrain_study_clean.csv"), row.names = F)

# create clean test data file
pretrain_test_clean <- filter(pretrain_test, !(id %in% pretrain_exclusions))
write.csv(pretrain_test_clean, file = here(clean_data_dir, "exp2_pretrain_test_clean.csv"), row.names = F)

# ENCODING
# ----------------

enc_data <- read.csv(here(raw_data_dir, "/exp2_enc.csv"), stringsAsFactors = F)

# continue excluding participants that were previously excluded
enc_exclusions <- pretrain_exclusions

# exclude participants who don't make at least 75% of encoding responses
enc_check <- enc_data %>%
  filter(startsWith(display, "motor_trial")) %>%
  group_by(id) %>%
  summarise(elim_trials = sum(is.na(response)),
            elim_prop = elim_trials / length(response)) %>%
  ungroup()
enc_exclusions <- c(enc_exclusions, enc_check$id[enc_check$elim_prop > 0.25])

# exclude participants with poor pretraining performance (< 80% correct on final test)
pretrain_final_perf <- pretrain_test_clean %>%
  filter(test_rep == 6) %>%
  group_by(id) %>%
  summarise(final_acc = mean(accuracy)) %>%
  ungroup()
enc_exclusions <- c(enc_exclusions, pretrain_final_perf$id[pretrain_final_perf$final_acc < 0.8])
enc_exclusions <- unique(enc_exclusions)

# create clean encoding data file
enc_clean <- enc_data %>%
  filter(!(id %in% enc_exclusions), !(is.na(response)),
         !(display %in% c("instructions", "fixation", "event_start")))
write.csv(enc_clean, file = here(clean_data_dir, "exp2_enc_clean.csv"), row.names = F)

# ORDER RECONSTRUCTION
# ----------------

recon_data <- read.csv(here(raw_data_dir, "/exp2_recon.csv"), stringsAsFactors = F)

# continue excluding participants that were previously excluded
recon_exclusions <- enc_exclusions

# exclude participants not above chance performance on this memory test
recon_binom <- recon_data %>%
  group_by(id) %>%
  summarise(n_correct = sum(accuracy == 1, na.rm = T),
            n_total = sum(!is.na(response)))

recon_binom[, c("binom_stat", "binom_pval")] <- NA
for (sub in recon_binom$id) {
  bt <- binom.test(filter(recon_binom, id == sub)$n_correct, filter(recon_binom, id == sub)$n_total,
                   p = 1/6, alternative = 'g')
  recon_binom$binom_stat[recon_binom$id == sub] <- bt$statistic
  recon_binom$binom_pval[recon_binom$id == sub] <- bt$p.value
}
recon_exclusions <- c(recon_exclusions, recon_binom$id[recon_binom$binom_pval > 0.05])
recon_exclusions <- unique(recon_exclusions)

# create clean recon file
recon_clean <- recon_data %>%
  filter(!(id %in% recon_exclusions), rt >= 100, !is.na(response))
write.csv(recon_clean, file = here(clean_data_dir, "exp2_recon_clean.csv"), row.names = F)

sprintf('%d/80 subjects remaining for reconstruction analyses', length(unique(recon_clean$id)))

# SPATIAL MEMORY
# ----------------

spatial_data <- read.csv(here(raw_data_dir, "/exp2_spatial.csv"), stringsAsFactors = F)

# continue excluding participants that were previously excluded
spatial_exclusions <- enc_exclusions

# exclude participants not above chance performance on this memory test
spatial_binom <- spatial_data %>%
  group_by(id) %>%
  summarise(n_correct = sum(accuracy == 1, na.rm = T),
            n_total = sum(!is.na(response)))

spatial_binom[, c("binom_stat", "binom_pval")] <- NA
for (sub in spatial_binom$id) {
  bt <- binom.test(filter(spatial_binom, id == sub)$n_correct, filter(spatial_binom, id == sub)$n_total,
                   p = 1/4, alternative = 'g')
  spatial_binom$binom_stat[spatial_binom$id == sub] <- bt$statistic
  spatial_binom$binom_pval[spatial_binom$id == sub] <- bt$p.value
}
spatial_exclusions <- c(spatial_exclusions, spatial_binom$id[spatial_binom$binom_pval > 0.05])
spatial_exclusions <- unique(spatial_exclusions)

# create clean spatial file
spatial_clean <- spatial_data %>%
  filter(!(id %in% spatial_exclusions), rt >= 100, !is.na(response))
write.csv(spatial_clean, file = here(clean_data_dir, "exp2_spatial_clean.csv"), row.names = F)

# FINAL RECON
# ----------------

## ==============================
## EXPERIMENT 3
## ==============================

# PRETRAINING
# ----------------

pretrain_data <- read.csv(here(raw_data_dir, "/exp3_pretrain.csv"), stringsAsFactors = F)

# split into study + test trials
pretrain_study <- filter(pretrain_data, startsWith(display, "seq_trial"))
pretrain_test <- filter(pretrain_data, display == "test_trial")

# exclude participants who don't make at least 75% of study responses
pretrain_study_check <- pretrain_study %>%
  group_by(id) %>%
  summarise(elim_trials = sum(is.na(response)),
            elim_prop = elim_trials / length(response)) %>%
  ungroup()
pretrain_exclusions <- pretrain_study_check$id[pretrain_study_check$elim_prop > 0.25]

# exclude two RTs that were impossibly long
# trial response limit was 6s, so these must be due to a glitch in the exp. website
pretrain_study <- filter(pretrain_study, rt < 8000)

# exclude study trials with outlier RTs (> 3 SDs above the mean)
study_rts <- filter(pretrain_study, !is.na(response))$rt
sd_cutoff <- mean(study_rts) + (3 * sd(study_rts))

# create clean study data file
pretrain_study_clean <- filter(pretrain_study, !(id %in% pretrain_exclusions), !is.na(response))
write.csv(pretrain_study_clean, file = here(clean_data_dir, "exp3_pretrain_study_clean.csv"), row.names = F)

# create clean test data file
pretrain_test_clean <- filter(pretrain_test, !(id %in% pretrain_exclusions))
write.csv(pretrain_test_clean, file = here(clean_data_dir, "exp3_pretrain_test_clean.csv"), row.names = F)

# ENCODING
# ----------------

enc_data <- read.csv(here(raw_data_dir, "/exp3_enc.csv"), stringsAsFactors = F)

# continue excluding participants that were previously excluded
enc_exclusions <- pretrain_exclusions

# exclude participants who don't make at least 75% of encoding responses
enc_check <- enc_data %>%
  filter(startsWith(display, "motor_trial")) %>%
  group_by(id) %>%
  summarise(elim_trials = sum(is.na(response)),
            elim_prop = elim_trials / length(response)) %>%
  ungroup()
enc_exclusions <- c(enc_exclusions, enc_check$id[enc_check$elim_prop > 0.25])

# exclude participants with poor pretraining performance (< 80% correct on final test)
pretrain_final_perf <- pretrain_test_clean %>%
  filter(test_rep == 6) %>%
  group_by(id) %>%
  summarise(final_acc = mean(accuracy)) %>%
  ungroup()
enc_exclusions <- c(enc_exclusions, pretrain_final_perf$id[pretrain_final_perf$final_acc < 0.8])
enc_exclusions <- unique(enc_exclusions)

# create clean encoding data file
enc_clean <- enc_data %>%
  filter(!(id %in% enc_exclusions), !(is.na(response)),
         !(display %in% c("instructions", "fixation", "event_start")))
write.csv(enc_clean, file = here(clean_data_dir, "exp3_enc_clean.csv"), row.names = F)

# ORDER RECONSTRUCTION
# ----------------

recon_data <- read.csv(here(raw_data_dir, "/exp3_recon.csv"), stringsAsFactors = F)

# continue excluding participants that were previously excluded
recon_exclusions <- enc_exclusions

# exclude participants with poor encoding performance (for predictable/from memory trials)
enc_pred <- enc_clean %>%
  filter(startsWith(display, 'motor_trial'), condition == "pred") %>%
  group_by(id) %>%
  summarise(pred_acc = mean(accuracy)) %>%
  ungroup()
recon_exclusions <- c(recon_exclusions, enc_pred$id[enc_pred$pred_acc < 0.8])

# exclude participants not above chance performance on this memory test
recon_binom <- recon_data %>%
  group_by(id) %>%
  summarise(n_correct = sum(accuracy == 1, na.rm = T),
            n_total = sum(!is.na(response)))

recon_binom[, c("binom_stat", "binom_pval")] <- NA
for (sub in recon_binom$id) {
  bt <- binom.test(filter(recon_binom, id == sub)$n_correct, filter(recon_binom, id == sub)$n_total,
                   p = 1/6, alternative = 'g')
  recon_binom$binom_stat[recon_binom$id == sub] <- bt$statistic
  recon_binom$binom_pval[recon_binom$id == sub] <- bt$p.value
}
recon_exclusions <- c(recon_exclusions, recon_binom$id[recon_binom$binom_pval > 0.05])
recon_exclusions <- unique(recon_exclusions)

# create clean recon file
recon_clean <- recon_data %>%
  filter(!(id %in% recon_exclusions), rt >= 100, !is.na(response))
write.csv(recon_clean, file = here(clean_data_dir, "exp3_recon_clean.csv"), row.names = F)

# ITEM RECOGNITION
# ----------------

recog_data <- read.csv(here(raw_data_dir, "/exp3_recog.csv"), stringsAsFactors = F)

# continue excluding participants that were previously excluded
recog_exclusions <- enc_exclusions

# exclude participants with poor encoding performance (for predictable/from memory trials)
recog_exclusions <- c(recog_exclusions, enc_pred$id[enc_pred$pred_acc < 0.8])

# exclude participants not above chance performance on this memory test
recog_binom <- recog_data %>%
  group_by(id) %>%
  summarise(n_correct = sum(accuracy == 1, na.rm = T),
            n_total = sum(!is.na(response)))

recog_binom[, c("binom_stat", "binom_pval")] <- NA
for (sub in recog_binom$id) {
  bt <- binom.test(filter(recog_binom, id == sub)$n_correct, filter(recog_binom, id == sub)$n_total,
                   p = 1/2, alternative = 'g')
  recog_binom$binom_stat[recog_binom$id == sub] <- bt$statistic
  recog_binom$binom_pval[recog_binom$id == sub] <- bt$p.value
}
recog_exclusions <- c(recog_exclusions, recog_binom$id[recog_binom$binom_pval > 0.05])
recog_exclusions <- unique(recog_exclusions)

# create clean recog file
recog_clean <- recog_data %>%
  filter(!(id %in% recog_exclusions), rt >= 100, !is.na(response))
write.csv(recog_clean, file = here(clean_data_dir, "exp3_recog_clean.csv"), row.names = F)

# REMINDER
# ----------------
remind_data <- read.csv(here(raw_data_dir, "/exp3_reminder.csv"), stringsAsFactors = F)

# continue excluding participants that were previously excluded (prior to memory test analysis)
remind_exclusions <- c(enc_exclusions, enc_pred$id[enc_pred$pred_acc < 0.8])

# create clean reminder file (and add repetition info)
remind_clean <- remind_data %>%
  filter(!(id %in% remind_exclusions)) %>%
  group_by(id) %>%
  mutate(test_rep = rep(c(1,2), each = 6)) %>%
  ungroup()

write.csv(remind_clean, file = here(clean_data_dir, "exp3_remind_clean.csv"), row.names = F)
