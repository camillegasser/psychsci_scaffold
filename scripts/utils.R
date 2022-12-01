# useful custom functions

# ============== 1 ==============
#'@title get_mode
#'@description find the mode of an array; handles ties & non-numeric data
#'@param x vector of numbers/strings to get mode of
#'@return mode(s) of the input array

get_mode <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  return(ux[tab == max(tab)])
}

# ============== 2 ==============
#'@title ttest_multiple_fdr
#'@description custom function for computing FDR-corrected t-tests & CIs
#'@param df data in the format of one row per participant and each condition
#'@param equation equation in the format y ~ x
#'@param grouping name of grouping variable used partition into multiple t-tests
#'@param q threshold for FDR correction (default = 0.05)

ttest_multiple_fdr <- function(df, equation, grouping, q = 0.05) {
  all_groups <- sort(unique(df[[grouping]]))
  n_groups <- length(all_groups)
  pvals <- c()
  
  # conduct t-tests
  for (i in all_groups) {
    t <- t.test(formula = eval(parse(text = equation)),
                data = filter(df, .data[[grouping]] == i),
                paired = TRUE)
    pvals <- c(pvals, t$p.value)
  }
  
  # correct pvals for multiple comparisons
  pvals_adj <- p.adjust(pvals, method = "fdr")
  
  # compute FDR-corrected CIs (following Benjamini et al., 2005)
  m <- length(pvals_adj)
  R <- sum(pvals_adj < q)
  ci_level <- 1 - R * q/m
  
  if (is.infinite(ci_level) | ci_level == 1) {
    ci_level = 0.95
    print('Inf CI level; using 95%')
  }
  
  # show t-test, cohen's d, and corrected p-vals + CIs for all tests
  idx = 1
  for (i in all_groups) {
    cat("\n")
    print(sprintf('TEST: %s, %s = %d', equation, grouping, i))
    
    t <- t.test(formula = eval(parse(text = equation)),
                data = filter(df, .data[[grouping]] == i),
                paired = TRUE,
                conf.level = ci_level)
    print(t)
    print(sprintf('adjusted pval: %.4f', pvals_adj[idx]))
    
    d <- effsize::cohen.d(formula = eval(parse(text = paste0(equation, '| Subject(id)'))),
                          data = filter(df, .data[[grouping]] == i),
                          paired = TRUE)
    print(d)
    #print(sprintf('cohens d: %.2f', d$estimate))
    idx = idx + 1
  }
}

# ============== 3 ==============
#'@title wilcox_multiple_fdr
#'@description custom function for computing FDR-corrected Wilcoxon signed rank tests & CIs
#'@param df data in the format of one row per participant and each condition
#'@param equation equation in the format y ~ x
#'@param grouping name of grouping variable used partition into multiple tests
#'@param q threshold for FDR correction (default = 0.05)

wilcox_multiple_fdr <- function(df, equation, grouping, q = 0.05) {
  all_groups <- sort(unique(df[[grouping]]))
  n_groups <- length(all_groups)
  pvals <- c()
  
  # conduct wilcox tests
  for (i in all_groups) {
    w <- wilcox.test(formula = eval(parse(text = equation)),
                     data = filter(df, .data[[grouping]] == i),
                     paired = TRUE)
    pvals <- c(pvals, w$p.value)
  }
  
  # correct pvals for multiple comparisons
  pvals_adj <- p.adjust(pvals, method = "fdr")
  
  # show wilcox test, cliff's d, and corrected p-vals for all tests
  idx = 1
  for (i in all_groups) {
    print(sprintf('TEST: %s, %s = %d', equation, grouping, i))
    
    w <- wilcox.test(formula = eval(parse(text = equation)),
                     data = filter(df, .data[[grouping]] == i),
                     paired = TRUE)
    print(w)
    print(sprintf('adjusted pval: %.4f', pvals_adj[idx]))
    
    d <- effsize::cliff.delta(formula = eval(parse(text = equation)),
                              data = filter(df, .data[[grouping]] == i),
                              paired = TRUE)
    print(d)
    #print(sprintf('cliffs d: %.2f', d$estimate))
    idx = idx + 1
  }
}

# ============== 4 ==============
#'@title levenshtein
#'@description computes minimum number of edits (insertions, deletions, or substitutions) needed to
#              transform participants' order into correct order
#'@param x first array
#'@param y second array
#'@return levenshtein distance between x and y
#'credit: https://stackoverflow.com/questions/46794378/levenshtein-distance-in-r-on-sentence-level

levenshtein <- function(x, y){
  unique_words <- unique(c(x,y))
  letter_x <- plyr::mapvalues(x,
                              from = unique_words,
                              to = letters[1:length(unique_words)])
  letter_y <- plyr::mapvalues(y,
                              from = unique_words,
                              to = letters[1:length(unique_words)])
  return(adist(paste0(letter_x, collapse = ''), paste0(letter_y, collapse = '')))
}

# ============== 5 ==============
#'@title get_leven_dist
#'@description custom function for creating dataframe to compute levenshtein distance
#'@param recon_df trial-level dataframe for order reconstruction test
#'@return event-level dataframe with levenshtein distance included

get_leven_dist <- function(recon_df) {
  recon_event <- recon_df %>%
    select(id, list, event, condition, seq_pos, response, recon_correct) %>%
    pivot_wider(names_from = seq_pos, values_from = c(response, recon_correct), names_prefix = 'pos') %>%
    lapply(function(x) { gsub('.png','',x) }) %>%
    as.data.frame()
  
  # merge list of items into a single dataframe cell
  recon_event <- recon_event %>%
    group_by(id, list, event, condition) %>%
    unite(resp_seq, response_pos1, response_pos2, response_pos3,
          response_pos4, response_pos5, response_pos6, sep = ', ', remove = F) %>%
    unite(corr_seq,  recon_correct_pos1, recon_correct_pos2, recon_correct_pos3,
          recon_correct_pos4, recon_correct_pos5, recon_correct_pos6, sep = ', ', remove = F) %>%
    select(c(id, list, event, condition, resp_seq, corr_seq)) %>%
    ungroup()
  
  # compute levenshtein distance (not very fast)
  for (i in 1:nrow(recon_event)) {
    s1 <- unlist(strsplit(recon_event[i,]$resp_seq, ', '))
    s2 <- unlist(strsplit(recon_event[i,]$corr_seq, ', '))
    # leven_dist <- levenshtein(s1, s2)[1,1]
    leven_dist <- levenshtein(s1[s1 != "NA"], s2)[1,1]
    recon_event[i,'leven_dist'] <- leven_dist
  }
  
  return(recon_event)
}

