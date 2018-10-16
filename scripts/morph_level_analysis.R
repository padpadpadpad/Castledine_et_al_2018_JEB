#--------------------------------------------------#
# Castledine et al. 2018 - morph-level analysis ####
#--------------------------------------------------#

# clear workspace
mise::mise(vars = TRUE, figs = TRUE, pkgs = TRUE, console = TRUE)

# load in packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)
library(stringr)
library(lme4)
library(MicrobioUoE) # to install, run: devtools::install_github('padpadpadpad/MicrobioUoE')
library(emmeans)
library(gridExtra)
library(readr)

# figure path
path_fig <- 'figures'

# read in data ####
d <- read.csv('data/all_trials_morphs_counts.csv', stringsAsFactors = FALSE)

#-------------------------------#
# transform and wrangle data ####
#-------------------------------#

# 1. filter out rows where T0 and T1 = NA
d <- filter(d, !is.na(T0_count_ml)) %>%
  group_by(., comb, morph) %>%
  # 2. count the number of each morph per trial
  mutate(., morph_num = n()) %>%
  # 3. calculate number of other morphs per trial
  ungroup() %>%
  group_by(., comb) %>%
  # 4. add column for total clones pre trial
  mutate(., total_clones = n()) %>%
  ungroup() %>%
  # 5. filter for where total clones = 1 or 3
  filter(., total_clones < 4 & total_clones >= 2) %>%
  # 6. calculate number of other clones of other morph per trial
  group_by(., comb, morph) %>%
  mutate(., T0_other = sum(T0_count_ml) - T0_count_ml,
         T1_other = sum(T1_count_ml) - T1_count_ml) %>%
  ungroup() %>%
  # 7. filter for instances where competing morph is the same as the morph counts
  filter(., T0_other > 0) %>%
  # 8. rename columns
  rename(., T0_LZ = T0_count_ml, 
         T1_LZ = T1_count_ml,
         T0_WT = T0_other,
         T1_WT = T1_other) %>%
  # 9. filter out so just LZ types are kept
  filter(., type == 'LZ') %>%
  # 10. calculate fitness - relative fitness of LZ compared to WT of same morph
  group_by(., comb) %>%
  mutate(., rel_fitness = log(T1_LZ/T0_LZ) / log(T1_WT/T0_WT)) %>%
  select(., comb, treatment, treatment2, morph, WT_div, LZ_div, WT_coev, LZ_coev, total_clones, rel_fitness) %>%
  ungroup() %>%
  # 11. rename columns of coevolved
  rename(., resident = LZ_coev, invader = WT_coev) %>%
  # 12. add columns for the coevolutionary history of the resident (lacz) and the invader (wild-type)
  mutate(., resident = case_when(resident == 'N' & LZ_div == 1 ~ 'alone',
                              resident == 'N' & LZ_div == 2 ~ 'uncoev',
                              TRUE ~ 'coev'),
         invader = case_when(invader == 'N' & WT_div == 1 ~ 'alone',
                             invader == 'N' & WT_div == 2 ~ 'uncoev',
                             TRUE ~ 'coev'),
        rel_fitness = ifelse(is.infinite(rel_fitness), NA, rel_fitness))

# check if all trials are present
d_check <- group_by(d, treatment2) %>%
  summarise(count = n()) %>%
  arrange(count)

#----------------------------------#
# load in information of blocks ####
#----------------------------------#

# read in and transform block data
d_pairs <- read.csv('data/all_combs.csv', stringsAsFactors = FALSE) %>%
  janitor::clean_names() %>%
  # filter for combinations in d
  filter(., combination %in% unique(d$comb)) %>%
  # grab number from WT column and add column for blocking
  mutate(wt_comm = readr::parse_number(wt),
         block = case_when(wt_comm %in% c(1, 7) ~ '1',
                           wt_comm %in% c(9, 10) ~ '2',
                           TRUE ~ '3')) %>%
  rename(., comb = combination) %>%
  separate(., lacz, c('lacz1', 'lacz2'), sep = '_') %>%
  gather(., 'lacz', 'morph', c(lacz1, lacz2)) %>%
  select(., comb, block, morph) %>% 
  mutate(., lacz_comm = readr::parse_number(morph),
         morph = stringr::str_extract(morph, '[:alpha:]{1,}'))

# merge d and d_pairs
d <- merge(d, d_pairs, by = c('comb', 'morph'), all.x = TRUE)

# create column for LZ_div vs WT_div
d <- unite(d, id, c(LZ_div, WT_div), sep = '_', remove = FALSE)

# create means of each treatment combination to create a balanced analysis
d_means <- group_by(d, block, morph, lacz_comm, id, WT_div, LZ_div, invader, resident, total_clones) %>%
  summarise(., rel_fitness = mean(rel_fitness)) %>%
  ungroup()

# add column for facet labels
d_means <- mutate(d_means, facet_labs = ifelse(morph == 'SM', '(a) smooth morph', '(b) wrinkly-spreader morph'))

#----------------------------------------------------------#
# model 1 - does addition of extra morph change outcome ####
#----------------------------------------------------------#

mod1 <- lmer(rel_fitness ~ id * morph + (1|block), d_means, na.action = na.fail, REML = 'FALSE')
mod2 <- lmer(rel_fitness ~ id + morph + (1|block), d_means, na.action = na.fail, REML = 'FALSE')
anova(mod1, mod2)
# keep model 1

# pairwise contrasts per morph
emmeans::emmeans(mod1, pairwise ~ id|morph)

# average over all other effects
emmeans::emmeans(mod1, ~morph)

#----------------------------------------------#
# test for effect of coevolution on outcome ####
#----------------------------------------------#

# Have to test separately for lacz and wild-type otherwise would have impossible factor combinations 

# LZ coevolutionary history
LZ_two <- filter(d_means, total_clones == 3 & LZ_div == 2)
LZ_mod1 <- lmer(rel_fitness ~ resident * morph + (1|block), LZ_two)
LZ_mod2 <- lmer(rel_fitness ~ resident + morph + (1|block), LZ_two)
LZ_mod3 <- lmer(rel_fitness ~ morph + (1|block), LZ_two)
LZ_mod4 <- lmer(rel_fitness ~ 1 + (1|block), LZ_two)
anova(LZ_mod1, LZ_mod2)
anova(LZ_mod2, LZ_mod3)
anova(LZ_mod3, LZ_mod4)
# relative between morphs is significant but nothing else is

# WT coevolutionary history
WT_two <- filter(d_means, total_clones == 3 & WT_div == 2)
WT_mod1 <- lmer(rel_fitness ~ invader*morph + (1|block), WT_two)
WT_mod2 <- lmer(rel_fitness ~ invader+morph + (1|block), WT_two)
WT_mod3 <- lmer(rel_fitness ~ morph + (1|block), WT_two)
WT_mod4 <- lmer(rel_fitness ~ 1 + (1|block), WT_two)
anova(WT_mod1, WT_mod2)
anova(WT_mod2, WT_mod3)
anova(WT_mod3, WT_mod4)
# No morph effect for wild-type fitness

#------------------#
# make Figure 2 ####
#------------------#

# set up new text for x labels
x_lab_text_WS <- c(expression(atop('LZ WS', 'WT WS')),
                expression(atop('LZ SM and WS', 'WT WS')),
                expression(atop('LZ WS', 'WT SM and WS')))
x_lab_text_SM <- c(expression(atop('LZ SM', 'WT SM')),
                   expression(atop('LZ SM and WS', 'WT SM')),
                   expression(atop('LZ SM', 'WT SM and WS')))

# add column for coevolution shape
d_means <- mutate(d_means, shape = case_when(WT_div == 1 & LZ_div == 1 ~ 'alone',
                                             WT_div == 1 & LZ_div == 2 & resident == 'coev' ~ 'coev',
                                             LZ_div == 1 & WT_div == 2 & invader == 'coev' ~ 'coev',
                                             TRUE ~ 'uncoev'))

d_means <- mutate(d_means, letter = case_when(shape == 'alone' ~ '',
                                              shape == 'coev' ~ 'C',
                                              TRUE ~ 'R'),
                  x_axis = as.factor(paste(LZ_div, WT_div, sep = '_')))
d_means$x_axis <- forcats::fct_relevel(d_means$x_axis, '1_2', after = Inf)
d_means$x_axis2 = jitter(as.numeric(d_means$x_axis))


plot_WS <- ggplot(filter(d_means, morph == 'WS'), aes(x_axis, rel_fitness)) +
  scale_x_discrete(labels = x_lab_text_WS) +
  geom_hline(aes(yintercept = 1), linetype = 2) +
  geom_pretty_boxplot(col = 'black', fill = 'black') +
  geom_point(aes(x = x_axis2), shape = 21, fill = 'white', size = 5) +
  geom_text(aes(x = x_axis2, label = letter), size = 3) +
  ylab(expression(atop('Relative fitness of LacZ morph', 'compared to same wild-type morph'))) +
  xlab('Functional group combination') +
  ylim(c(0.4, 2.25)) +
  theme_bw(base_family = 'Times', base_size = 14) +
  scale_shape_manual(values = c(21, 22, 24)) +
  theme(legend.position = 'none') +
  ggtitle('(b) wrinkly-spreader morph')

plot_SM <- ggplot(filter(d_means, morph == 'SM'), aes(x_axis, rel_fitness)) +
  scale_x_discrete(labels = x_lab_text_WS) +
  geom_hline(aes(yintercept = 1), linetype = 2) +
  geom_pretty_boxplot(col = 'black', fill = 'black') +
  geom_point(aes(x = x_axis2), shape = 21, fill = 'white', size = 5) +
  geom_text(aes(x = x_axis2, label = letter), size = 3) +
  ylab(expression(atop('Relative fitness of LacZ morph', 'compared to same wild-type morph'))) +
  xlab('Functional group combination') +
  ylim(c(0.4, 2.25)) +
  scale_shape_manual(values = c(21, 22, 24)) +
  theme_bw(base_family = 'Times', base_size = 14) +
  theme(legend.position = 'none') +
  ggtitle('(a) smooth morph')

plot_both <- grid.arrange(plot_SM + ylab('') + xlab(''), plot_WS + ylab('') + xlab(''), ncol = 2)

# save plot
ggsave(file.path(path_fig, 'Figure_2.pdf'), plot_both, height = 6, width = 11)
ggsave(file.path(path_fig, 'Figure_2.png'), plot_both, height = 6, width = 11)


