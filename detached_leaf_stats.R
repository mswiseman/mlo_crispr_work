
library(openxlsx)
library(dplyr)
library(stringr)
library(ggplot2)
library(lme4)           
library(lmerTest)     
library(emmeans)        
library(multcomp)      
library(broom)          
library(broom.mixed)    

# ---------------------------
# 1) Load + clean
# ---------------------------
file_path <- "2025_detached_leaf_ratings_Nov_19.xlsx"
dat <- read.xlsx(file_path, sheet = "Sheet1")

dat <- dat %>%
  mutate(across(where(is.character), ~ gsub("[ \t]", "", .)))

names(dat) <- tolower(names(dat))

rename_if_present <- function(df, pattern, new_name) {
  cand <- grep(pattern, names(df), ignore.case = TRUE, value = TRUE)
  if (length(cand) == 1 && !(new_name %in% names(df))) {
    names(df)[names(df) == cand] <- new_name
  }
  df
}

dat <- dat %>%
  rename_if_present("^basal[_ ]?line$", "basal_line") %>%
  rename_if_present("^clonal[_ ]?line$|^clone$|^clonal$", "clonal_line") %>%
  rename_if_present("^run$", "run") %>%
  rename_if_present("^leaf[_ ]?no$|^leafno$|^leaf$", "leaf_no") %>%
  rename_if_present("^disease[_ ]?severity$", "disease_severity")

req <- c("basal_line", "clonal_line", "date_rated", "leaf_no", "disease_severity")
missing_cols <- setdiff(req, names(dat))
stopifnot(length(missing_cols) == 0)

dat <- dat %>%
  mutate(
    basal_line       = factor(basal_line),
    clonal_line      = factor(clonal_line),
    run              = factor(run),
    date_rated       = factor(date_rated),
    leaf_no          = factor(leaf_no),
    disease_severity = as.numeric(disease_severity)
  )

# ---------------------------
# 2) Filter once
# ---------------------------
dat2 <- dat %>%
  filter(!basal_line %in% c("Symphony")) %>%
  group_by(clonal_line) %>%
    filter(n_distinct(date_rated) >= 2) %>%  
  ungroup() %>%
  droplevels()

# Helper for consistent plots
plot_means_with_letters <- function(emm_obj, title, out_svg) {
  cld_tbl <- multcomp::cld(emm_obj, adjust = "tukey", Letters = letters) %>%
    as.data.frame() %>%
    as_tibble() %>%
    transmute(
      basal_line = as.character(basal_line),
      letter = gsub("\\s+", "", .group)
    )

  emm_tbl <- as.data.frame(emm_obj) %>%
    as_tibble() %>%
    mutate(
      basal_line = as.character(basal_line),
      lower = emmean - qt(0.975, df = df) * SE,
      upper = emmean + qt(0.975, df = df) * SE
    ) %>%
    left_join(cld_tbl, by = "basal_line") %>%
    arrange(emmean) %>%
    mutate(basal_line = factor(basal_line, levels = basal_line))

  p <- ggplot(emm_tbl, aes(x = emmean, y = basal_line)) +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
    geom_point(size = 2.5) +
    geom_text(aes(label = letter), nudge_x = 2, size = 3) +
    theme_classic() +
    labs(
      x = "Estimated mean disease severity (± 95% CI)",
      y = "Basal line",
      title = title
    )

  print(p)
  ggsave(out_svg, p, height = 5, width = 7, dpi = 300)
  invisible(list(emm_tbl = emm_tbl, plot = p))
}

plot_diff_vs_wt <- function(ct_obj, title, out_svg) {
  df <- as.data.frame(ct_obj) %>%
    as_tibble() %>%
    mutate(
      basal_line = str_split_fixed(contrast, " - ", 2)[, 1],
      control    = str_split_fixed(contrast, " - ", 2)[, 2],
      lower      = estimate - qt(0.975, df = df) * SE,
      upper      = estimate + qt(0.975, df = df) * SE
    ) %>%
    arrange(estimate) %>%
    mutate(basal_line = factor(basal_line, levels = unique(basal_line)))

  p <- ggplot(df, aes(x = estimate, y = basal_line)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
    geom_point(size = 2.5) +
    theme_classic() +
    labs(
      x = "Difference vs WT (Estimate ± 95% CI)",
      y = "Basal line",
      title = title
    )

  print(p)
  ggsave(out_svg, p, height = 5, width = 7, dpi = 300)
  invisible(list(pairs_df = df, plot = p))
}

# ---------------------------
# 3) Split into TWO datasets (no cross-background comparisons)
# ---------------------------
dat_fuggle <- dat2 %>% filter(str_detect(as.character(basal_line), "Fuggle")) %>% droplevels()
dat_tett   <- dat2 %>% filter(str_detect(as.character(basal_line), "Tettnanger")) %>% droplevels()

# ---------------------------
# 4) Fit TWO separate (L)MMs
# Note: This is NOT a GLMM (because disease_severity is numeric/continuous).
# It's a linear mixed model. If you truly want GLMM, tell me the distribution/link.
# ---------------------------
m_fuggle <- lmer(disease_severity ~ basal_line + (1 | date_rated) + (1 | clonal_line),
                 data = dat_fuggle, REML = TRUE)

m_tett <- lmer(disease_severity ~ basal_line + (1 | date_rated) + (1 | clonal_line),
               data = dat_tett, REML = TRUE)

print(summary(m_fuggle))
print(summary(m_tett))

# ---------------------------
# 5) EMMs + within-background letters
# ---------------------------
emm_fuggle <- emmeans(m_fuggle, ~ basal_line)
emm_tett   <- emmeans(m_tett,   ~ basal_line)

res_f_means <- plot_means_with_letters(
  emm_fuggle,
  title = "Fuggle basal lines: estimated means + Tukey letters (within Fuggle only)",
  out_svg = "basal_fuggle_means_letters.svg"
)

res_t_means <- plot_means_with_letters(
  emm_tett,
  title = "Tettnanger basal lines: estimated means + Tukey letters (within Tettnanger only)",
  out_svg = "basal_tettnanger_means_letters.svg"
)

# ---------------------------
# 6) Within-background WT comparisons (separately)
# ---------------------------
ctrl_F <- "Fuggle_wt"
stopifnot(ctrl_F %in% levels(dat_fuggle$basal_line))
ct_fuggle <- contrast(emm_fuggle, method = "trt.vs.ctrl", ref = ctrl_F, adjust = "holm")

# pick Tettnanger WT robustly
ctrl_T_candidates <- levels(dat_tett$basal_line)[
  str_detect(levels(dat_tett$basal_line), "Tettnanger") & str_detect(levels(dat_tett$basal_line), "wt|WT")
]
stopifnot(length(ctrl_T_candidates) > 0)
ctrl_T <- ctrl_T_candidates[1]
ct_tett <- contrast(emm_tett, method = "trt.vs.ctrl", ref = ctrl_T, adjust = "holm")

print(ct_fuggle)
print(ct_tett)

res_f_diff <- plot_diff_vs_wt(
  ct_fuggle,
  title = "Fuggle basal lines vs Fuggle WT (Holm-adjusted; within Fuggle only)",
  out_svg = "basal_fuggle_diff_vs_wt.svg"
)

res_t_diff <- plot_diff_vs_wt(
  ct_tett,
  title = "Tettnanger basal lines vs Tettnanger WT (Holm-adjusted; within Tettnanger only)",
  out_svg = "basal_tettnanger_diff_vs_wt.svg"
)

#list(
#  m_fuggle = m_fuggle,
#  m_tett   = m_tett,
#  emm_fuggle = emm_fuggle,
#  emm_tett   = emm_tett,
#  fuggle_means_tbl = res_f_means$emm_tbl,
#  tett_means_tbl   = res_t_means$emm_tbl,
#  fuggle_vs_wt_tbl = res_f_diff$pairs_df,
#  tett_vs_wt_tbl   = res_t_diff$pairs_df
#)

######################### CLONAL ANALYSIS ######################### 

dat_clonal <- dat %>%
  filter(!basal_line %in% c("Symphony")) %>%
  group_by(clonal_line) %>%
  filter(n_distinct(date_rated) >= 2) %>%  
  mutate(
    clone_id  = factor(clonal_line),
    ctrl_group = if_else(str_detect(as.character(clonal_line), "^G14_Tettnanger"), "Tettnanger", "Fuggle")
  )

# ---------------------------
# 3) Fit mixed model
# ---------------------------
m_clone <- lmer(disease_severity ~ clone_id + (1 | date_rated), data = dat_clonal, REML = TRUE)
emm_clone <- emmeans(m_clone, ~ clone_id)

# ---------------------------
# 4) Two-control contrasts (robust via emmeans grid)
# ---------------------------
grid <- as.data.frame(emm_clone)
grid_ids <- as.character(grid$clone_id)

ctrl_F <- "Fuggle_wt"
ctrl_T <- "Tettnanger_wt"

# Rule: clones that start with G14-Tettnanger OR G14_Tettnanger use Tettnanger WT
is_tet_clone <- function(x) str_detect(x, "^G14[-_]Tettnanger")

# Build contrast vectors: each clone minus its assigned WT
is_ctrl  <- grid_ids %in% c(ctrl_F, ctrl_T)
test_ids <- grid_ids[!is_ctrl]

contrast_list <- setNames(vector("list", length(test_ids)), test_ids)

for (cl in test_ids) {
  v <- rep(0, length(grid_ids))
  names(v) <- grid_ids
  v[cl] <- 1
  v[ if (is_tet_clone(cl)) ctrl_T else ctrl_F ] <- -1
  contrast_list[[cl]] <- v
}

pairs_vs_wt <- contrast(emm_clone, method = contrast_list, adjust = "holm")
pairs_df <- as.data.frame(pairs_vs_wt) %>%
  as_tibble() %>%
  mutate(
    # parse safely: split on " - " (space-hyphen-space), not internal hyphens
    clone   = str_split_fixed(contrast, " - ", 2)[, 1],
    control = str_split_fixed(contrast, " - ", 2)[, 2],
    lower   = estimate - qt(0.975, df = df) * SE,
    upper   = estimate + qt(0.975, df = df) * SE,
    sig = case_when(
      p.value < 0.001 ~ "***",
      p.value < 0.01  ~ "**",
      p.value < 0.05  ~ "*",
      TRUE            ~ ""
    )
  )

# quick sanity check: any duplicate clone rows?
dup_check <- pairs_df %>% count(clone) %>% filter(n > 1)
if (nrow(dup_check) > 0) {
  warning("Some clone names appear multiple times in contrasts:\n",
          paste(dup_check$clone, collapse = ", "))
}

pairs_df

# ---------------------------
# 5) Plot
# ---------------------------
plot_contr <- pairs_df %>%
  arrange(estimate) %>%
  mutate(clone = factor(clone, levels = unique(clone)))

ggplot(plot_contr, aes(x = estimate, y = clone)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  geom_point(size = 2.5) +
  geom_text(aes(label = sig), nudge_x = 1, size = 3) +
  theme_classic() +
  labs(
    x = "Difference in Disease Severity vs assigned WT (Estimate ± 95% CI)",
    y = "Clonal Replicate ID",
    title = "Clone Replicate Lines Compared to WT (Mixed Model; Holm-adjusted)"
  )

ggsave("grouped_diff_to_wt_clonal.svg", height = 5, width = 7, dpi = 300)

plot_contr