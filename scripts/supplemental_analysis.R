rm(list = ls())

# change to project directory
setwd("/Volumes/GoogleDrive/My Drive/grad-school/scaffold_behav/manuscript/psychsci_submission/psychsci_scaffold")

# setup
library(Rmisc)
library(tidyverse)
library(ggbeeswarm)
library(effsize)
library(lme4)
library(psycho)
library(rstatix)
source("scripts/utils.R")

# set plot theme
plot_theme <- theme_light() + 
  theme(panel.grid = element_blank(),
        legend.position = 'none')
theme_set(plot_theme)
cond_cols <- c('salmon2','dodgerblue2')

data_dir <- "clean_data"
raw_data_dir <- "raw_data"
plot_dir <- "plots"

# ----------------
# PRETRAINING RTs BY SEQUENCE POSITION
# ----------------
# supplemental fig. 1

# EXPERIMENT 1

# load data and remove outlier RTs
pretrain_study <- read.csv(paste0(data_dir, "/exp1_pretrain_study_clean.csv"), stringsAsFactors = F)
all_rts <- filter(pretrain_study, !is.na(response))$rt
sd_cutoff <- mean(all_rts) + (3 * sd(all_rts))
pretrain_study <- filter(pretrain_study, rt <= sd_cutoff)

# mean RTs by condition
pretrain_study_rt <- pretrain_study %>%
  filter(block == 3) %>%
  group_by(id, condition, seq_pos) %>%
  summarise(mean_rt = mean(rt)) %>%
  ungroup()

# summarize across participants
pretrain_study_rt_group <- pretrain_study_rt %>%
  group_by(condition, seq_pos) %>%
  summarise(group_mean_rt = mean(mean_rt)) %>%
  ungroup()
pretrain_study_rt_group$group_sem_rt <- summarySEwithin(data = pretrain_study, measurevar = "rt",
                                                        withinvars = c("condition", "seq_pos"),
                                                        idvar = "id")$se

# plot
ggplot(pretrain_study_rt_group, aes(x = seq_pos, y = group_mean_rt, color = condition)) +
  geom_jitter(data = pretrain_study_rt, aes(x = seq_pos, y = mean_rt, fill = condition),
              width = 0.1, height = 0.05, color = 'white',
              alpha = 0.4, shape = 21, size = 3) +
  geom_line(size = 1.2) +
  geom_errorbar(width = 0, size = 1, aes(x = seq_pos,
                                         ymin = group_mean_rt - group_sem_rt,
                                         ymax = group_mean_rt + group_sem_rt)) +
  scale_color_manual(values = cond_cols) +
  scale_fill_manual(values = cond_cols) +
  scale_x_continuous(breaks = 1:6) +
  ylim(0, 2200) +
  labs(x = 'sequence position', y = 'RT (ms)')
ggsave(paste0(plot_dir, "/supp_exp1_pretrainRT_seqpos.pdf"), width = 4.5, height = 3.5)

# stats
# some subs excluded bc after exclusions, they don't have enough trials in each condition/seq pos bin
usable_subs <- pretrain_study_rt %>% group_by(id) %>% summarise(count = length(id))
usable_subs <- usable_subs$id[usable_subs$count == max(usable_subs$count)]
wilcox_multiple_fdr(subset(pretrain_study_rt, id %in% usable_subs),
                    equation = 'mean_rt ~ condition', grouping = 'seq_pos')

# EXPERIMENT 2

# load data and remove outlier RTs
pretrain_study <- read.csv(paste0(data_dir, "/exp2_pretrain_study_clean.csv"), stringsAsFactors = F)
all_rts <- filter(pretrain_study, !is.na(response))$rt
sd_cutoff <- mean(all_rts) + (3 * sd(all_rts))
pretrain_study <- filter(pretrain_study, rt <= sd_cutoff)

# mean RTs by condition
pretrain_study_rt <- pretrain_study %>%
  filter(block == 3) %>%
  group_by(id, condition, seq_pos) %>%
  summarise(mean_rt = mean(rt)) %>%
  ungroup()

# summarize across participants
pretrain_study_rt_group <- pretrain_study_rt %>%
  group_by(condition, seq_pos) %>%
  summarise(group_mean_rt = mean(mean_rt)) %>%
  ungroup()
pretrain_study_rt_group$group_sem_rt <- summarySEwithin(data = pretrain_study, measurevar = "rt",
                                                        withinvars = c("condition", "seq_pos"),
                                                        idvar = "id")$se

# plot
ggplot(pretrain_study_rt_group, aes(x = seq_pos, y = group_mean_rt, color = condition)) +
  geom_jitter(data = pretrain_study_rt, aes(x = seq_pos, y = mean_rt, fill = condition),
              width = 0.1, height = 0.05, color = 'white',
              alpha = 0.4, shape = 21, size = 3) +
  geom_line(size = 1.2) +
  geom_errorbar(width = 0, size = 1, aes(x = seq_pos,
                                         ymin = group_mean_rt - group_sem_rt,
                                         ymax = group_mean_rt + group_sem_rt)) +
  scale_color_manual(values = cond_cols) +
  scale_fill_manual(values = cond_cols) +
  scale_x_continuous(breaks = 1:6) +
  ylim(0,2200) +
  labs(x = 'sequence position', y = 'RT (ms)')
ggsave(paste0(plot_dir, "/supp_exp2_pretrainRT_seqpos.pdf"), width = 4.5, height = 3.5)

# stats
# some subs excluded bc after exclusions, they don't have enough trials in each condition/seq pos bin
usable_subs <- pretrain_study_rt %>% group_by(id) %>% summarise(count = length(id))
usable_subs <- usable_subs$id[usable_subs$count == max(usable_subs$count)]
wilcox_multiple_fdr(subset(pretrain_study_rt, id %in% usable_subs),
                    equation = 'mean_rt ~ condition', grouping = 'seq_pos')

# EXPERIMENT 3

# load data and remove outlier RTs
pretrain_study <- read.csv(paste0(data_dir, "/exp3_pretrain_study_clean.csv"), stringsAsFactors = F)
all_rts <- filter(pretrain_study, !is.na(response))$rt
sd_cutoff <- mean(all_rts) + (3 * sd(all_rts))
pretrain_study <- filter(pretrain_study, rt <= sd_cutoff)

# mean RTs by condition
pretrain_study_rt <- pretrain_study %>%
  filter(block == 3) %>%
  group_by(id, condition, seq_pos) %>%
  summarise(mean_rt = mean(rt)) %>%
  ungroup()

# summarize across participants
pretrain_study_rt_group <- pretrain_study_rt %>%
  group_by(condition, seq_pos) %>%
  summarise(group_mean_rt = mean(mean_rt)) %>%
  ungroup()
pretrain_study_rt_group$group_sem_rt <- summarySEwithin(data = pretrain_study, measurevar = "rt",
                                                        withinvars = c("condition", "seq_pos"),
                                                        idvar = "id")$se

# plot
ggplot(pretrain_study_rt_group, aes(x = seq_pos, y = group_mean_rt, color = condition)) +
  geom_jitter(data = pretrain_study_rt, aes(x = seq_pos, y = mean_rt, fill = condition),
              width = 0.1, height = 0.05, color = 'white',
              alpha = 0.4, shape = 21, size = 3) +
  geom_line(size = 1.2) +
  geom_errorbar(width = 0, size = 1, aes(x = seq_pos,
                                         ymin = group_mean_rt - group_sem_rt,
                                         ymax = group_mean_rt + group_sem_rt)) +
  scale_color_manual(values = cond_cols) +
  scale_fill_manual(values = cond_cols) +
  scale_x_continuous(breaks = 1:6) +
  ylim(0,2200) +
  labs(x = 'sequence position', y = 'RT (ms)')
ggsave(paste0(plot_dir, "/supp_exp3_pretrainRT_seqpos.pdf"), width = 4.5, height = 3.5)

# stats
# some subs excluded bc after exclusions, they don't have enough trials in each condition/seq pos bin
usable_subs <- pretrain_study_rt %>% group_by(id) %>% summarise(count = length(id))
usable_subs <- usable_subs$id[usable_subs$count == max(usable_subs$count)]
wilcox_multiple_fdr(subset(pretrain_study_rt, id %in% usable_subs),
                    equation = 'mean_rt ~ condition', grouping = 'seq_pos')


# ----------------
# ENCODING ACCURACY
# ----------------
# supplemental fig. 2

# EXPERIMENT 1

enc <- read.csv(paste0(data_dir, "/exp1_enc_clean.csv"), stringsAsFactors = F)

# mean accuracy by condition
enc_acc <- enc %>%
  group_by(id, condition) %>%
  summarise(mean_acc = mean(accuracy)) %>%
  ungroup()

# summarize across participants
enc_acc_group <- enc_acc %>%
  group_by(condition) %>%
  summarise(group_mean_acc = mean(mean_acc)) %>%
  ungroup()
enc_acc_group$group_sem_acc <- summarySEwithin(data = enc_acc, measurevar = "mean_acc",
                                               withinvars = "condition", idvar = "id")$se

# plot
ggplot(enc_acc, aes(x = condition, y = mean_acc, fill = condition)) +
  geom_bar(data = enc_acc_group, aes(x = condition, y = group_mean_acc),
           stat = 'identity', width = 0.8) +
  geom_beeswarm(data = enc_acc, aes(x = condition, y = mean_acc),
                dodge.width = 1, groupOnX = T, size = 4,
                pch = 21, color = 'white') +
  geom_errorbar(data = enc_acc_group,
                aes(x = condition, y = group_mean_acc,
                    ymin = group_mean_acc - group_sem_acc,
                    ymax = group_mean_acc + group_sem_acc),
                width = 0.1) +
  scale_fill_manual(values = cond_cols) +
  labs(x = 'condition', y = 'aisle response accuracy') +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0,1.1))
ggsave(paste0(plot_dir, "/supp_exp1_encAcc.pdf"), width = 4.5, height = 3.5)

# stats
t.test(mean_acc ~ condition, data = enc_acc, paired = T)
effsize::cohen.d(mean_acc ~ condition  | Subject(id), data = enc_acc, paired = T)

# mean accuracy by condition and block number ("list")
enc_acc <- enc %>%
  group_by(id, condition, list) %>%
  summarise(mean_acc = mean(accuracy)) %>%
  ungroup()

# summarize across participants
enc_acc_group <- enc_acc %>%
  group_by(condition, list) %>%
  summarise(group_mean_acc = mean(mean_acc)) %>%
  ungroup()
enc_acc_group$group_sem_acc <- summarySEwithin(data = enc_acc, measurevar = "mean_acc",
                                               withinvars = c("condition", "list"), idvar = "id")$se

# plot
ggplot(enc_acc_group, aes(x = list, y = group_mean_acc, color = condition)) +
  geom_jitter(data = enc_acc, aes(x = list, y = mean_acc, fill = condition),
              width = 0.1, height = 0.02, color = 'white',
              alpha = 0.4, shape = 21, size = 2.3) +
  geom_line(size = 1.5) +
  geom_errorbar(aes(x = list,
                    ymin = group_mean_acc - group_sem_acc,
                    ymax = group_mean_acc + group_sem_acc),
                width = 0, size = 1) +
  scale_color_manual(values = cond_cols) +
  scale_fill_manual(values = cond_cols) +
  scale_x_continuous(breaks = 1:8) +
  labs(x = 'block number', y = 'aisle response accuracy')
ggsave(paste0(plot_dir, "/supp_exp1_encAcc_block.pdf"), width = 5, height = 3.5)

# stats
ttest_multiple_fdr(enc_acc, equation = 'mean_acc ~ condition', grouping = 'list')

# EXPERIMENT 3

enc <- read.csv(paste0(data_dir, "/exp3_enc_clean.csv"), stringsAsFactors = F)

# mean accuracy by condition
enc_acc <- enc %>%
  group_by(id, condition) %>%
  summarise(mean_acc = mean(accuracy)) %>%
  ungroup()

# summarize across participants
enc_acc_group <- enc_acc %>%
  group_by(condition) %>%
  summarise(group_mean_acc = mean(mean_acc)) %>%
  ungroup()
enc_acc_group$group_sem_acc <- summarySEwithin(data = enc_acc, measurevar = "mean_acc",
                                               withinvars = "condition", idvar = "id")$se

# plot
ggplot(enc_acc, aes(x = condition, y = mean_acc, fill = condition)) +
  geom_bar(data = enc_acc_group, aes(x = condition, y = group_mean_acc),
           stat = 'identity', width = 0.8) +
  geom_beeswarm(data = enc_acc, aes(x = condition, y = mean_acc),
                dodge.width = 1, groupOnX = T, size = 4,
                pch = 21, color = 'white') +
  geom_errorbar(data = enc_acc_group,
                aes(x = condition, y = group_mean_acc,
                    ymin = group_mean_acc - group_sem_acc,
                    ymax = group_mean_acc + group_sem_acc),
                width = 0.1) +
  scale_fill_manual(values = cond_cols) +
  labs(x = 'condition', y = 'aisle response accuracy') +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0,1.1))
ggsave(paste0(plot_dir, "/supp_exp3_encAcc.pdf"), width = 4.5, height = 3.5)

# stats
t.test(mean_acc ~ condition, data = enc_acc, paired = T)
effsize::cohen.d(mean_acc ~ condition  | Subject(id), data = enc_acc, paired = T)

# mean accuracy by condition and block number ("list")
enc_acc <- enc %>%
  group_by(id, condition, list) %>%
  summarise(mean_acc = mean(accuracy)) %>%
  ungroup()

# summarize across participants
enc_acc_group <- enc_acc %>%
  group_by(condition, list) %>%
  summarise(group_mean_acc = mean(mean_acc)) %>%
  ungroup()
enc_acc_group$group_sem_acc <- summarySEwithin(data = enc_acc, measurevar = "mean_acc",
                                               withinvars = c("condition", "list"), idvar = "id")$se

# plot
ggplot(enc_acc_group, aes(x = list, y = group_mean_acc, color = condition)) +
  geom_jitter(data = enc_acc, aes(x = list, y = mean_acc, fill = condition),
              width = 0.1, height = 0.02, color = 'white',
              alpha = 0.4, shape = 21, size = 2.3) +
  geom_line(size = 1.5) +
  geom_errorbar(aes(x = list,
                    ymin = group_mean_acc - group_sem_acc,
                    ymax = group_mean_acc + group_sem_acc),
                width = 0, size = 1) +
  scale_color_manual(values = cond_cols) +
  scale_fill_manual(values = cond_cols) +
  scale_x_continuous(breaks = 1:8) +
  labs(x = 'block number', y = 'aisle response accuracy')
ggsave(paste0(plot_dir, "/supp_exp3_encAcc_block.pdf"), width = 5, height = 3.5)

# stats
ttest_multiple_fdr(enc_acc, equation = 'mean_acc ~ condition', grouping = 'list')

# ----------------
# ORDER MEMORY: MULTILEVEL MODEL CONTROLLING FOR CONFOUNDS
# ----------------
# supplemental tables 1-3

n_iter <- 8000
bayes_output_dir <- 'bmm_fits'

# EXPERIMENT 1

recon <- read.csv(paste0(data_dir, "/exp1_recon_clean.csv"), stringsAsFactors = F)

bayes_fn <- paste0(bayes_output_dir, '/exp1_recon_model_bayes.rda')
if (file.exists(bayes_fn)) {
  load(file = paste0(bayes_fn))
  
} else {
  recon_model <- recon %>%
    select(id, list, event, event_name, condition, seq_pos, accuracy) %>%
    mutate(condition_d = ifelse(condition == 'pred', 0.5, -0.5), # pred = 0.5
           cat_e = ifelse(event_name == 'grocery store.png', 0.5, -0.5), # food = 0.5
           seq_pos_c = scale(seq_pos, center = T, scale = F), # centered
           list_c = scale(list, center = T, scale = F), # centered
           event_c = scale(event, center = T, scale = F)) # centered
  
  recon_bm <- brm(accuracy ~ condition_d + list_c + cat_e + seq_pos_c + event_c +
                    (condition_d + list_c + cat_e + seq_pos_c + event_c || id),
                  data = recon_model, family = bernoulli(link = 'logit'),
                  iter = n_iter, seed = 123)
  save(recon_bm, file = paste0(bayes_output_dir, '/exp1_recon_model_bayes.rda'))
}
print(recon_bm, digits = 3)

# EXPERIMENT 2

recon <- read.csv(paste0(data_dir, "/exp2_recon_clean.csv"), stringsAsFactors = F)

bayes_fn <- paste0(bayes_output_dir, '/exp2_recon_model_bayes.rda')
if (file.exists(bayes_fn)) {
  load(file = paste0(bayes_fn))
  
} else {
  recon_model <- recon %>%
    select(id, list, event, event_name, condition, seq_pos, accuracy) %>%
    mutate(condition_d = ifelse(condition == 'pred', 0.5, -0.5), # pred = 0.5
           cat_e = ifelse(event_name == 'grocery store.png', 0.5, -0.5), # food = 0.5
           seq_pos_c = scale(seq_pos, center = T, scale = F), # centered
           list_c = scale(list, center = T, scale = F), # centered
           event_c = scale(event, center = T, scale = F)) # centered
  recon_bm <- brm(accuracy ~ condition_d + list_c + cat_e + seq_pos_c + event_c +
                    (condition_d + list_c + cat_e + seq_pos_c + event_c || id),
                  data = recon_model, family = bernoulli(link = 'logit'),
                  iter = n_iter, seed = 123)
  save(recon_bm, file = paste0(bayes_output_dir, '/exp2_recon_model_bayes.rda'))
}

print(recon_bm, digits = 3)

# EXPERIMENT 3

recon <- read.csv(paste0(data_dir, "/exp3_recon_clean.csv"), stringsAsFactors = F)

bayes_fn <- paste0(bayes_output_dir, '/exp3_recon_model_bayes.rda')
if (file.exists(bayes_fn)) {
  load(file = paste0(bayes_fn))
  
} else {
  recon_model <- recon %>%
    select(id, list, event, event_name, condition, seq_pos, accuracy) %>%
    mutate(condition_d = ifelse(condition == 'pred', 0.5, -0.5), # pred = 0.5
           cat_e = ifelse(event_name == 'grocery store.png', 0.5, -0.5), # food = 0.5
           seq_pos_c = scale(seq_pos, center = T, scale = F), # centered
           list_c = scale(list, center = T, scale = F), # centered
           event_c = scale(event, center = T, scale = F)) # centered
  recon_bm <- brm(accuracy ~ condition_d + list_c + cat_e + seq_pos_c + event_c +
                    (condition_d + list_c + cat_e + seq_pos_c + event_c || id),
                  data = recon_model, family = bernoulli(link = 'logit'),
                  iter = n_iter, seed = 123)
  save(recon_bm, file = paste0(bayes_output_dir, '/exp3_recon_model_bayes.rda'))
}

print(recon_bm, digits = 3)

# ----------------
# ENCODING ACCURACY & ORDER RECON MEMORY
# ----------------
# supplemental analysis: "Order memory and aisle response accuracy during encoding"

# EXPERIMENT 1

recon <- read.csv(paste0(data_dir, "/exp1_recon_clean.csv"), stringsAsFactors = F)
enc <- read.csv(paste0(data_dir, "/exp1_enc_clean.csv"), stringsAsFactors = F)

# add encoding accuracy to data frame
enc_resp <- enc %>%
  select(id, list, event, condition, stim, accuracy) %>%
  rename(enc_acc = accuracy, recon_correct = stim)

recon_enc <- left_join(recon, enc_resp, by = c('id', 'list', 'recon_correct','event', 'condition'))

# mean accuracy by condition, correct encoding responses only
recon_enc <- recon_enc %>%
  filter(enc_acc == 1) %>%
  group_by(id, condition) %>%
  summarise(mean_acc = mean(accuracy)) %>%
  ungroup()
recon_enc_acc_group <- recon_enc %>%
  group_by(condition) %>%
  summarise(group_mean_acc = mean(mean_acc)) %>%
  ungroup()
recon_enc_acc_group$group_sem_acc <- summarySEwithin(data = recon_enc, measurevar = "mean_acc",
                                                     withinvars = "condition", idvar = "id")$se

# recompute condition difference
t.test(mean_acc ~ condition, data = recon_enc, paired = T)
effsize::cohen.d(mean_acc ~ condition | Subject(id), data = recon_enc, paired = T)

# EXPERIMENT 3

recon <- read.csv(paste0(data_dir, "/exp3_recon_clean.csv"), stringsAsFactors = F)
enc <- read.csv(paste0(data_dir, "/exp3_enc_clean.csv"), stringsAsFactors = F)

# add encoding accuracy to data frame
enc_resp <- enc %>%
  select(id, list, event, condition, stim, accuracy) %>%
  rename(enc_acc = accuracy, recon_correct = stim)

recon_enc <- left_join(recon, enc_resp, by = c('id', 'list', 'recon_correct','event', 'condition'))

# mean accuracy by condition, correct encoding responses only
recon_enc <- recon_enc %>%
  filter(enc_acc == 1) %>%
  group_by(id, condition) %>%
  summarise(mean_acc = mean(accuracy)) %>%
  ungroup()
recon_enc_acc_group <- recon_enc %>%
  group_by(condition) %>%
  summarise(group_mean_acc = mean(mean_acc)) %>%
  ungroup()
recon_enc_acc_group$group_sem_acc <- summarySEwithin(data = recon_enc, measurevar = "mean_acc",
                                                     withinvars = "condition", idvar = "id")$se

# recompute condition difference
t.test(mean_acc ~ condition, data = recon_enc, paired = T)
effsize::cohen.d(mean_acc ~ condition | Subject(id), data = recon_enc, paired = T)

# ----------------
# AISLE REPETITION & ORDER RECON MEMORY
# ----------------
# supplemental analysis: "Order memory and aisle repetition during encoding"

# EXPERIMENT 1

recon <- read.csv(paste0(data_dir, "/exp1_recon_clean.csv"), stringsAsFactors = F)
enc <- read.csv(paste0(raw_data_dir, "/exp1_enc.csv"), stringsAsFactors = F)
# loading non-cleaned encoding file because we need all trials to get aisle repetition info

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

recon <- read.csv(paste0(data_dir, "/exp2_recon_clean.csv"), stringsAsFactors = F)
enc <- read.csv(paste0(raw_data_dir, "/exp2_enc.csv"), stringsAsFactors = F)
# loading non-cleaned encoding file because we need all trials to get aisle repetition info

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
recon <- read.csv(paste0(data_dir, "/exp3_recon_clean.csv"), stringsAsFactors = F)
enc <- read.csv(paste0(raw_data_dir, "/exp3_enc.csv"), stringsAsFactors = F)
# loading non-cleaned encoding file because we need all trials to get aisle repetition info

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
  summarise(group_mean_acc = mean(mean_acc)) %>%
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
# SPATIAL MEMORY & AISLE LOCATION
# ----------------
# supplemental analysis: "Spatial memory performance and aisle location"

# EXPERIMENT 1

spatial <- read.csv(paste0(data_dir, "/exp1_spatial_clean.csv"), stringsAsFactors = F)

# accuracy by aisle location
spatial_acc <- spatial %>%
  group_by(id, spatial_correct) %>%
  summarise(mean_acc = mean(accuracy)) %>%
  ungroup() %>%
  convert_as_factor(id, spatial_correct) %>%
  select(id, spatial_correct, mean_acc) %>%
  as.data.frame()

# run anova + post-hoc ttests
av <- anova_test(data = spatial_acc, dv = mean_acc, wid = id,
                 within = spatial_correct, effect.size = 'pes')
av_tab <- get_anova_table(av)
av_tab

# pes (partial eta sq) confidence interval
get.ci.partial.eta.squared(av_tab$F, av_tab$DFn, av_tab$DFd, conf.level = 0.95)

# post hoc t tests  
ps <- pairwise_t_test(formula = mean_acc ~ spatial_correct, paired = T,
                      p.adjust.method = 'fdr', data = spatial_acc)
print(ps)

# get adjusted CI + effect size for significant comparisons
# compute corrected CIs (following Benjamini et al., 2005)
ps_adj <- ps$p.adj
q <- 0.05
m <- length(ps_adj)
R <- sum(ps_adj < q)
ci_level <- 1 - R * q/m

t.test(mean_acc ~ spatial_correct, paired = T, conf.level = ci_level,
       data = subset(spatial_acc, spatial_correct %in% c(1, 2)))
effsize::cohen.d(subset(spatial_acc, spatial_correct == 1)$mean_acc,
                 subset(spatial_acc, spatial_correct == 2)$mean_acc, paired = T)

t.test(mean_acc ~ spatial_correct, paired = T, conf.level = ci_level,
       data = subset(spatial_acc, spatial_correct %in% c(1, 3)))
effsize::cohen.d(subset(spatial_acc, spatial_correct == 1)$mean_acc,
                 subset(spatial_acc, spatial_correct == 3)$mean_acc, paired = T)

t.test(mean_acc ~ spatial_correct, paired = T, conf.level = ci_level,
       data = subset(spatial_acc, spatial_correct %in% c(1, 4)))
effsize::cohen.d(subset(spatial_acc, spatial_correct == 1)$mean_acc,
                 subset(spatial_acc, spatial_correct == 4)$mean_acc, paired = T)

# EXPERIMENT 2

spatial <- read.csv(paste0(data_dir, "/exp2_spatial_clean.csv"), stringsAsFactors = F)

# accuracy by aisle location
spatial_acc <- spatial %>%
  group_by(id, spatial_correct) %>%
  summarise(mean_acc = mean(accuracy)) %>%
  ungroup() %>%
  convert_as_factor(id, spatial_correct) %>%
  select(id, spatial_correct, mean_acc) %>%
  as.data.frame()

# run anova + post-hoc ttests
av <- anova_test(data = spatial_acc, dv = mean_acc, wid = id,
                 within = spatial_correct, effect.size = 'pes')
av_tab <- get_anova_table(av)
av_tab

# pes (partial eta sq) confidence interval
get.ci.partial.eta.squared(av_tab$F, av_tab$DFn, av_tab$DFd, conf.level = 0.95)

# post hoc t tests  
ps <- pairwise_t_test(formula = mean_acc ~ spatial_correct, paired = T,
                      p.adjust.method = 'fdr', data = spatial_acc)
print(ps)

# get adjusted CI (+ effect size) for significant comparisons
# compute corrected CIs (following Benjamini et al., 2005)
ps_adj <- ps$p.adj
q <- 0.05
m <- length(ps_adj)
R <- sum(ps_adj < q)
ci_level <- 1 - R * q/m

t.test(mean_acc ~ spatial_correct, paired = T, conf.level = ci_level,
       data = subset(spatial_acc, spatial_correct %in% c(1, 2)))
effsize::cohen.d(subset(spatial_acc, spatial_correct == 1)$mean_acc,
                 subset(spatial_acc, spatial_correct == 2)$mean_acc, paired = T)

t.test(mean_acc ~ spatial_correct, paired = T, conf.level = ci_level,
       data = subset(spatial_acc, spatial_correct %in% c(1, 3)))
effsize::cohen.d(subset(spatial_acc, spatial_correct == 1)$mean_acc,
                 subset(spatial_acc, spatial_correct == 3)$mean_acc, paired = T)

t.test(mean_acc ~ spatial_correct, paired = T, conf.level = ci_level,
       data = subset(spatial_acc, spatial_correct %in% c(1, 4)))
effsize::cohen.d(subset(spatial_acc, spatial_correct == 1)$mean_acc,
                 subset(spatial_acc, spatial_correct == 4)$mean_acc, paired = T)

t.test(mean_acc ~ spatial_correct, paired = T, conf.level = ci_level,
       data = subset(spatial_acc, spatial_correct %in% c(2, 4)))
effsize::cohen.d(subset(spatial_acc, spatial_correct == 2)$mean_acc,
                 subset(spatial_acc, spatial_correct == 4)$mean_acc, paired = T)

t.test(mean_acc ~ spatial_correct, paired = T, conf.level = ci_level,
       data = subset(spatial_acc, spatial_correct %in% c(3, 4)))
effsize::cohen.d(subset(spatial_acc, spatial_correct == 3)$mean_acc,
                 subset(spatial_acc, spatial_correct == 4)$mean_acc, paired = T)
# ----------------
# ITEM RECOGNITION BY CONFIDENCE
# ----------------
# supplemental figure 3a

recog <- read.csv(paste0(data_dir, "/exp3_recog_clean.csv"), stringsAsFactors = F)

# dprime by condition and confidence
recog_acc <- recog %>%
  group_by(id, condition, confidence) %>%
  summarise(n_hits = sum(recog_correct == 'old' & response_bin == 'old'),
            n_fa = sum(recog_correct == 'new' & response_bin == 'old'),
            n_cr = sum(recog_correct == 'new' & response_bin == 'new'),
            n_miss = sum(recog_correct == 'old' & response_bin == 'new'),
            n_old = sum(recog_correct == 'old'),
            n_new = sum(recog_correct == 'new')) %>%
  mutate(hit_rate = n_hits / n_old,
         fa_rate = n_fa / n_new,
         cr_rate = n_cr / n_new,
         miss_rate = n_miss / n_old) %>%
  filter(!any(is.na(hit_rate), is.na(fa_rate), is.na(cr_rate), is.na(miss_rate))) %>%
  mutate(dprime = dprime(n_hit = n_hits, n_fa = n_fa, n_cr = n_cr, n_miss = n_miss,
                         adjusted = T)$dprime) %>%
  ungroup()
recog_acc_group <- recog_acc %>%
  group_by(condition, confidence) %>%
  summarise(group_mean_dprime = mean(dprime))
recog_acc_group$group_sem_dprime <- summarySEwithin(data = recog_acc, measurevar = "dprime",
                                                    withinvars = c("condition", "confidence"),
                                                    idvar = "id")$se
# plot
ggplot(recog_acc_group, aes(x = condition, y = group_mean_dprime, fill = condition)) +
  geom_bar(stat = 'identity', position = position_dodge(), width = 0.8) +
  geom_hline(yintercept = 0, color = 'grey20') +
  geom_dotplot(data = recog_acc, aes(x = condition, y = dprime), dotsize = 1.3,
               binaxis = 'y', stackdir = 'center', alpha = 1, color = 'white') + 
  geom_errorbar(aes(x = condition, ymin = group_mean_dprime - group_sem_dprime,
                    ymax = group_mean_dprime + group_sem_dprime),
                width = 0.1) +
  labs(x = 'condition', y = 'dprime') +
  scale_fill_manual(values = cond_cols) +
  facet_wrap(~ confidence)
ggsave(paste0(plot_dir, "/supp_exp3_recogAcc_conf.pdf"), width = 5.5, height = 3.5)

# stats
# some subs excluded bc they don't have enough trials in each condition/conf bin
usable_subs <- recog_acc %>% group_by(id) %>% summarise(count = length(id))
usable_subs <- usable_subs$id[usable_subs$count == max(usable_subs$count)]
t.test(dprime ~ condition, paired = T,
       data = subset(recog_acc, id %in% usable_subs & confidence == "high"))
t.test(dprime ~ condition, paired = T,
       data = subset(recog_acc, id %in% usable_subs & confidence == "low"))

# ----------------
# ITEM RECOGNITION BY SEQUENCE POSITION
# ----------------
# supplemental figure 3b

recog <- read.csv(paste0(data_dir, "/exp3_recog_clean.csv"), stringsAsFactors = F)

# dprime by condition and seq_pos
recog_acc <- recog %>%
  group_by(id, condition, seq_pos) %>%
  summarise(n_hits = sum(recog_correct == 'old' & response_bin == 'old'),
            n_fa = sum(recog_correct == 'new' & response_bin == 'old'),
            n_cr = sum(recog_correct == 'new' & response_bin == 'new'),
            n_miss = sum(recog_correct == 'old' & response_bin == 'new'),
            n_old = sum(recog_correct == 'old'),
            n_new = sum(recog_correct == 'new')) %>%
  mutate(hit_rate = n_hits / n_old,
         fa_rate = n_fa / n_new,
         cr_rate = n_cr / n_new,
         miss_rate = n_miss / n_old) %>%
  filter(!any(is.na(hit_rate), is.na(fa_rate), is.na(cr_rate), is.na(miss_rate))) %>%
  mutate(dprime = dprime(n_hit = n_hits, n_fa = n_fa, n_cr = n_cr, n_miss = n_miss,
                         adjusted = T)$dprime) %>%
  ungroup()
recog_acc_group <- recog_acc %>%
  group_by(condition, seq_pos) %>%
  summarise(group_mean_dprime = mean(dprime))
recog_acc_group$group_sem_dprime <- summarySEwithin(data = recog_acc, measurevar = "dprime",
                                                    withinvars = c("condition", "seq_pos"),
                                                    idvar = "id")$se
# plot
ggplot(recog_acc_group, aes(x = seq_pos, y = group_mean_dprime, color = condition)) +
  geom_hline(yintercept = 0, color = 'grey20') +
  geom_jitter(data = recog_acc, aes(x = seq_pos, y = dprime, fill = condition),
              width = 0.1, height = 0.05, color = 'white', alpha = 0.4, shape = 21, size = 3) +
  geom_line(size = 1.5) +
  geom_errorbar(width = 0, size = 1, aes(x = seq_pos,
                                         ymin = group_mean_dprime - group_sem_dprime,
                                         ymax = group_mean_dprime + group_sem_dprime)) +
  scale_color_manual(values = cond_cols) +
  scale_fill_manual(values = cond_cols) +
  scale_x_continuous(breaks = 1:6) +
  labs(x = 'sequence position', y = 'dprime')
ggsave(paste0(plot_dir, "/supp_exp3_recogAcc_seqpos.pdf"), width = 5, height = 3.5)

# stats
ttest_multiple_fdr(recog_acc, equation = 'dprime ~ condition', grouping = 'seq_pos')

# ----------------
# REMINDER PERFORMANCE
# ----------------
# supplemental analysis: "Reminder task performance"

# EXPERIMENT 1

reminder <- read.csv(paste0(data_dir, "/exp1_remind_clean.csv"), stringsAsFactors = F)

# accuracy during the second/last reminder repetition
remind_acc <- reminder %>%
  filter(test_rep == 2) %>%
  group_by(id) %>%
  summarise(mean_acc = mean(accuracy, na.rm = T)) %>%
  ungroup()

remind_group <- remind_acc %>%
  summarise(group_mean_acc = mean(mean_acc),
            group_sd_acc = sd(mean_acc))
remind_group

# EXPERIMENT 3

reminder <- read.csv(paste0(data_dir, "/exp3_remind_clean.csv"), stringsAsFactors = F)

# accuracy during the second/last reminder repetition
remind_acc <- reminder %>%
  filter(test_rep == 2) %>%
  group_by(id) %>%
  summarise(mean_acc = mean(accuracy, na.rm = T)) %>%
  ungroup()

remind_group <- remind_acc %>%
  summarise(group_mean_acc = mean(mean_acc),
            group_sd_acc = sd(mean_acc))
remind_group

# ----------------
# CROSS-EXPERIMENT COMPARISONS
# ----------------
# supplemental table 4

# order memory test comparisons

e1_recon <- read.csv(paste0(data_dir, "/exp1_recon_clean.csv"), stringsAsFactors = F)
e2_recon <- read.csv(paste0(data_dir, "/exp2_recon_clean.csv"), stringsAsFactors = F)
e3_recon <- read.csv(paste0(data_dir, "/exp3_recon_clean.csv"), stringsAsFactors = F)

# within each experiment, get accuracy diff (pred - rand) for each participant
e1_recon_diff <- e1_recon %>%
  group_by(id, condition) %>%
  summarise(mean_acc = mean(accuracy)) %>%
  pivot_wider(names_from = 'condition', values_from = 'mean_acc',
              names_prefix = 'acc_', id_cols = 'id') %>%
  mutate(acc_diff = acc_pred - acc_rand)

e2_recon_diff <- e2_recon %>%
  group_by(id, condition) %>%
  summarise(mean_acc = mean(accuracy)) %>%
  pivot_wider(names_from = 'condition', values_from = 'mean_acc',
              names_prefix = 'acc_', id_cols = 'id') %>%
  mutate(acc_diff = acc_pred - acc_rand)

e3_recon_diff <- e3_recon %>%
  group_by(id, condition) %>%
  summarise(mean_acc = mean(accuracy)) %>%
  pivot_wider(names_from = 'condition', values_from = 'mean_acc',
              names_prefix = 'acc_', id_cols = 'id') %>%
  mutate(acc_diff = acc_pred - acc_rand)

# EXP 1 v 2
var.test(e1_recon_diff$acc_diff, e2_recon_diff$acc_diff)
t.test(e1_recon_diff$acc_diff, e2_recon_diff$acc_diff, paired = F, var.equal = T)
effsize::cohen.d(e1_recon_diff$acc_diff, e2_recon_diff$acc_diff, paired = F)

# EXP 1 v 3
var.test(e1_recon_diff$acc_diff, e3_recon_diff$acc_diff)
t.test(e1_recon_diff$acc_diff, e3_recon_diff$acc_diff, paired = F, var.equal = T)
effsize::cohen.d(e1_recon_diff$acc_diff, e3_recon_diff$acc_diff, paired = F)

# EXP 2 v 3
var.test(e2_recon_diff$acc_diff, e3_recon_diff$acc_diff)
t.test(e2_recon_diff$acc_diff, e3_recon_diff$acc_diff, paired = F, var.equal = T)
effsize::cohen.d(e2_recon_diff$acc_diff, e3_recon_diff$acc_diff, paired = F)

# spatial memory test comparisons

e1_spatial <- read.csv(paste0(data_dir, "/exp1_spatial_clean.csv"), stringsAsFactors = F)
e2_spatial <- read.csv(paste0(data_dir, "/exp2_spatial_clean.csv"), stringsAsFactors = F)

# within each experiment, get accuracy diff (pred - rand) for each participant
e1_spatial_diff <- e1_spatial %>%
  group_by(id, condition) %>%
  summarise(mean_acc = mean(accuracy)) %>%
  pivot_wider(names_from = 'condition', values_from = 'mean_acc',
              names_prefix = 'acc_', id_cols = 'id') %>%
  mutate(acc_diff = acc_pred - acc_rand)

e2_spatial_diff <- e2_spatial %>%
  group_by(id, condition) %>%
  summarise(mean_acc = mean(accuracy)) %>%
  pivot_wider(names_from = 'condition', values_from = 'mean_acc',
              names_prefix = 'acc_', id_cols = 'id') %>%
  mutate(acc_diff = acc_pred - acc_rand)

# EXP 1 v 2
var.test(e1_spatial_diff$acc_diff, e2_spatial_diff$acc_diff)
t.test(e1_spatial_diff$acc_diff, e2_spatial_diff$acc_diff, paired = F, var.equal = T)
effsize::cohen.d(e1_spatial_diff$acc_diff, e2_spatial_diff$acc_diff, paired = F)

# ----------------
# FINAL RECON TEST
# ----------------
# analyses not detailed in manuscript

final_recon <- read.csv(paste0(raw_data_dir, "/exp2_final_recon.csv"), stringsAsFactors = F)

# get participants included in main experiment 2 recon analyses to remove exclusions
exp2_recon <- read.csv(paste0(data_dir, "/exp2_recon_clean.csv"), stringsAsFactors = F)
final_recon <- final_recon %>%
  filter(id %in% unique(exp2_recon$id),
         rt >= 100, !is.na(response))

# exclude those who got below-chance accuracy on this final test too
final_recon_binom <- final_recon %>%
  group_by(id) %>%
  summarise(n_correct = sum(accuracy == 1, na.rm = T),
            n_total = sum(!is.na(response)))
final_recon_binom[, c("binom_stat", "binom_pval")] <- NA

for (sub in final_recon_binom$id) {
  bt <- binom.test(filter(final_recon_binom, id == sub)$n_correct,
                   filter(final_recon_binom, id == sub)$n_total,
                   p = 1/6, alternative = 'g')
  final_recon_binom$binom_stat[final_recon_binom$id == sub] <- bt$statistic
  final_recon_binom$binom_pval[final_recon_binom$id == sub] <- bt$p.value
}
final_recon <- filter(final_recon, !(id %in% final_recon_binom$id[final_recon_binom$binom_pval > 0.05]))

# order memory by condition
final_recon_acc <- final_recon %>%
  group_by(id, condition) %>%
  summarise(mean_acc = mean(accuracy)) %>%
  ungroup()
final_recon_acc_group <- final_recon_acc %>%
  group_by(condition) %>%
  summarise(group_mean_acc = mean(mean_acc)) %>%
  ungroup()
final_recon_acc_group$group_sem_acc <- summarySEwithin(data = final_recon_acc, measurevar = "mean_acc",
                                                       withinvars = "condition", idvar = "id")$se

# plot
ggplot(final_recon_acc, aes(x = condition, y = mean_acc, fill = condition)) +
  geom_bar(data = final_recon_acc_group, aes(x = condition, y = group_mean_acc),
           stat = 'identity', width = 0.8) +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 1.5, color = 'white') +
  geom_errorbar(data = final_recon_acc_group, width = 0.1,
                aes(x = condition, y = group_mean_acc,
                    ymin = group_mean_acc - group_sem_acc, ymax = group_mean_acc + group_sem_acc)) +
  geom_hline(yintercept = 1/6, color = 'grey20', linetype = 'dashed') +
  scale_fill_manual(values = cond_cols) +
  labs(x = 'condition', y = 'ordinal accuracy (final test)') +
  theme(text= element_text(size = 15))

t.test(mean_acc ~ condition, data = final_recon_acc, paired = T)
effsize::cohen.d(mean_acc ~ condition | Subject(id), data = final_recon_acc, paired = T)

# order memory by event
final_recon_event <- final_recon %>%
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
final_recon_event_group <- final_recon_event %>%
  group_by(condition, num_correct) %>%
  summarise(group_mean_prop = mean(proportion)) %>%
  ungroup()
final_recon_event_group$group_sem_prop <- summarySEwithin(data = final_recon_event,
                                                          measurevar = "proportion",
                                                          withinvars = c("condition", "num_correct"),
                                                          idvar = "id")$se

# stats
ttest_multiple_fdr(final_recon_event, equation = "proportion ~ condition", grouping = "num_correct")

# order memory by condition and seq pos
final_recon_acc <- final_recon %>%
  group_by(id, condition, seq_pos) %>%
  summarise(mean_acc = mean(accuracy)) %>%
  ungroup()
final_recon_acc_group <- final_recon_acc %>%
  group_by(condition, seq_pos) %>%
  summarise(group_mean_acc = mean(mean_acc)) %>%
  ungroup()
final_recon_acc_group$group_sem_acc <- summarySEwithin(data = final_recon_acc, measurevar = "mean_acc",
                                                       withinvars = c("condition", "seq_pos"),
                                                       idvar = "id")$se
# stats - ttest
ttest_multiple_fdr(final_recon_acc, equation = "mean_acc ~ condition", grouping = "seq_pos")

# stats - logistic regression model
final_recon_model_data <- mutate(final_recon,
                                 condition_e = ifelse(condition == 'pred', 0.5, -0.5), # effect coded
                                 seq_pos_z = scale(seq_pos, center = T, scale = F)) # center pos = 0
final_recon_model <- glmer(accuracy ~ condition_e * seq_pos_z + (condition_e + seq_pos_z || id),
                           data = final_recon_model_data, family = 'binomial')
summary(final_recon_model)
confint(final_recon_model) # takes awhile to run
