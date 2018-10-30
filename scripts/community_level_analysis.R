#------------------------------------------------------#
# Castledine et al. 2018 - community-level analysis ####
#------------------------------------------------------#

# clear workspace
mise::mise(vars = TRUE, figs = TRUE, pkgs = TRUE, console = TRUE)

# load in packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(magrittr)
library(lme4)
library(MicrobioUoE) # to install, run: devtools::install_github('padpadpadpad/MicrobioUoE')
library(emmeans)
library(readr)

# figure path
path_fig <- 'figures'

# read in data ####
d <- read.csv('data/all_trials_morphs_counts.csv', stringsAsFactors = FALSE)

#-------------------------------#
# transform and wrangle data ####
#-------------------------------#

# LacZ is the resident community, WT is the invader community

# add columns to d
# 1. filter out rows where T0 and T1 = NA
d <- filter(d, !is.na(T0_count_ml)) %>%
  group_by(., comb) %>%
  # 2. add column for total clones in each microcosm
  dplyr::mutate(., total_clones = dplyr::n()) %>%
  ungroup() %>%
  # 3. filter out instances where there are four clones
  filter(., total_clones == 4) %>%
  # 4. sum abundances at the type (lacz or WT levels)
  group_by(., comb, type, WT_coev, LZ_coev) %>%
  summarise(., T0 = sum(T0_count_ml, na.rm = TRUE),
         T1 = sum(T1_count_ml, na.rm = TRUE)) %>%
  ungroup() %>%
  # 5. add new columns for T0_resident and T0_invader
  #    and T1_invader and T1_resident
  group_by(comb) %>%
  mutate(., T0_resident = T0,
         T0_invader = sum(T0) - T0_resident,
         T1_resident = T1,
         T1_invader = sum(T1) - T1_resident) %>%
  ungroup() %>%
  # 6. filter for just LacZ so only the "resident" community values are left
  filter(., type == 'LZ') %>%
  # 7. calculate fitness - relative fitness of LZ (resident) compared to WT (invader) of whole community
  group_by(., comb) %>%
  mutate(., rel_fitness = log(T1_resident/T0_resident) / log(T1_invader/T0_invader)) %>%
  ungroup() %>%
  # 8. select only desired columns
  select(., -starts_with('T')) %>%
  # 9. rename columns of coevolved
  rename(., resident = LZ_coev, invader = WT_coev) %>%
  # 10. change values of coevolved columns
  mutate_at(., c('resident', 'invader'), function(x)ifelse(x == 'Y', 'coev', 'uncoev'))

# check if all trials are present
d_check <- dplyr::group_by(d, invader, resident) %>%
  dplyr::summarise(count = dplyr::n()) %>%
  arrange(count)

#----------------------------------#
# load in information of blocks ####
#----------------------------------#

# load in information of the block effect
d_pairs <- read.csv('data/all_combs.csv', stringsAsFactors = FALSE) %>%
  janitor::clean_names() %>%
  # 1. filter for combinations in d_four
  filter(., combination %in% unique(d$comb)) %>%
  # 2. grab number from WT column
  mutate(wt_comm = parse_number(wt),
  # 3. add column for blocking based on WT community
         block = case_when(wt_comm %in% c(1, 7) ~ 1,
                           wt_comm %in% c(9, 10) ~ 2,
                           TRUE ~ 3)) %>%
  select(., combination, block) %>%
  rename(., comb = combination)

# bind data together
d <- merge(d, d_pairs, by = 'comb')

#--------------------------------------------------------------------------#
# model 2 - does coevolution alter the outcome of community coalescence ####
#--------------------------------------------------------------------------#

# run a linear mixed effect model with a random effect of block
mod1 <- lmer(rel_fitness ~ invader*resident + (1|block), d, REML = FALSE)
mod2 <- lmer(rel_fitness ~ invader + resident + (1|block), d, REML = FALSE)
anova(mod1, mod2)
# the interaction is significant

# do pairwise contrasts of this model
emmeans::emmeans(mod1, pairwise~invader*resident)

# there is no significant difference, but the most complex model is the one that is favoured through model selection
# model with no terms (rel_fitness ~ 1) favoured with all model combinations but will not get there using traditional model simplification

# average over all other effects to look at mean communtiy performance
emmeans::emmeans(mod1, ~1)

#------------------#
# make Figure 3 ####
#------------------#

# set up new text for x labels
x_lab_text <- c(expression(atop('LZ'[coev], '\nWT'[coev])),
                expression(atop('LZ'[random], '\nWT'[coev])),
                expression(atop('LZ'[coev], '\nWT'[random])),
                expression(atop('LZ'[random], '\nWT'[random])))

ggplot(d, aes(interaction(resident, invader), rel_fitness)) +
  geom_hline(aes(yintercept = 1), linetype = 2) +
  geom_pretty_boxplot(fill = 'black', col = 'black') +
  geom_point(position = position_jitter(width = 0.1), size = 5, shape = 21, fill = 'white') +
  theme_bw(base_family = 'Times', base_size = 14) +
  ylab('Relative success of LacZ compared to wild-type community') +
  xlab('Coevolutionary history of LacZ and wild-type communities') +
  scale_x_discrete(labels = x_lab_text) +
  ylim(c(0.6, 2))

# save plot
ggsave(file.path(path_fig, 'Figure_3.pdf'), last_plot(), height = 6, width = 7)
ggsave(file.path(path_fig, 'Figure_3.png'), last_plot(), height = 6, width = 7)
