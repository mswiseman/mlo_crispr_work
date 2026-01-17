library(c(
  "dplyr",
  "ggplot2",
  "emmeans",
  "multcomp",
  "patchwork",
  "survival",
  "survminer",
  "broom",
  "tidyverse",
  "janitor",
  "openxlsx",
  "glmmTMB"
))


# ---- 1) Read + prep ----
medium_levels <- paste0("SIM-", 1:16)

dat <- read.xlsx("fuggle_shoot_regen_data_for_dissertation.xlsx") %>%
  clean_names() %>%
  mutate(
    medium = factor(treatment, levels = medium_levels),
    run    = factor(run),
    week   = as.integer(week),
    plate  = if ("plate" %in% names(.)) factor(plate) else factor(interaction(medium, run, block, drop = TRUE))
  ) %>%
  filter(week %in% c(4, 8, 12)) %>%
  mutate(
    shoot_freq = if_else(total_callus_on_plate > 0,
                         calli_producing_shoots / total_callus_on_plate, NA_real_)
  )

dat12 <- dat %>% filter(week == 12)

# structural-zero media at week 12
zero_media <- dat12 %>%
  group_by(medium) %>%
  summarise(successes = sum(calli_producing_shoots, na.rm = TRUE), .groups = "drop") %>%
  filter(successes == 0) %>%
  pull(medium)

# ---- 2) Panel A model (binomial GLMM) ----
# bounded binomial data -> binomial
mA <- glmmTMB(
  cbind(calli_producing_shoots, total_callus_on_plate - calli_producing_shoots) ~
    medium + run + (1 | plate),
  family = binomial(),
  data = dat12
)

summary(mA)

# emmeans for non-zero media only
emmA_nonzero <- emmeans(
  mA, ~ medium,
  type = "response",
  at = list(medium = setdiff(medium_levels, as.character(zero_media)))
)

A_df <- as.data.frame(emmA_nonzero) %>%
  transmute(
    medium = factor(medium, levels = medium_levels),
    prob_pct = prob * 100,
    LCL_pct  = asymp.LCL * 100,
    UCL_pct  = asymp.UCL * 100
  )


A_letters <- cld(emmA_nonzero, Letters = letters, adjust = "sidak") %>%
  as.data.frame() %>%
  transmute(
    medium = factor(medium, levels = medium_levels),
    group  = gsub(" ", "", .group)
  )

# add structural zeros back in for plotting
A_zero <- tibble(
  medium   = factor(zero_media, levels = medium_levels),
  prob_pct = 0, LCL_pct = 0, UCL_pct = 0
)

A_zero_letters <- tibble(
  medium = factor(zero_media, levels = medium_levels),
  group  = "a"
)

A_plot <- bind_rows(A_df, A_zero) %>%
  left_join(bind_rows(A_letters, A_zero_letters), by = "medium") %>%
  arrange(medium)


pA <- ggplot(A_plot, aes(medium, prob_pct, color = medium)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = LCL_pct, ymax = UCL_pct), width = 0.15) +
  geom_text(aes(label = group, y = UCL_pct + 4), color = "black", size = 4) +
  labs(
    title = "A. Estimated Percent of Callus Producing Shoots\n(Per Plate at Week 12; GLMM)",
    x = "Shoot Induction Medium",
    y = "Estimated shoot regeneration (%) ± 95% CI"
  ) +
  coord_cartesian(ylim = c(0, max(A_plot$UCL_pct, na.rm = TRUE) + 12)) +
  scale_color_viridis_d(guide = "none") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ---- 3) Panel B model (NB GLM) ----
mB <- glmmTMB(
  total_num_shoots ~ medium + run,
  family = nbinom2(),
  data = dat12
)

emmB <- emmeans(mB, ~ medium, type = "response")

B_df <- as.data.frame(emmB) %>%
  transmute(
    medium = factor(medium, levels = medium_levels),
    mean_count = response,
    LCL = asymp.LCL,
    UCL = asymp.UCL
  ) %>%
  filter(is.finite(UCL))

B_letters <- cld(emmB, Letters = letters, adjust = "sidak") %>%
  as.data.frame() %>%
  transmute(
    medium = factor(medium, levels = medium_levels),
    group  = gsub(" ", "", .group)
  )



B_plot <- left_join(B_df, B_letters, by = "medium")

pB <- ggplot(dat12, aes(medium, total_num_shoots)) +
  geom_violin(aes(fill = medium), trim = FALSE, alpha = 0.85, show.legend = FALSE) +
  geom_jitter(aes(color = medium), width = 0.2, height = 0,
              alpha = 0.65, size = 1.2, show.legend = FALSE) +
  geom_pointrange(
    data = B_plot,
    aes(y = mean_count, ymin = LCL, ymax = UCL),
    inherit.aes = TRUE,
    size = 0.25
  ) +
  geom_text(
    data = B_plot,
    aes(y = UCL + 0.8, label = group),
    inherit.aes = TRUE,
    color = "black",
    size = 4
  ) +
  scale_fill_viridis_d(guide = "none") +
  scale_color_viridis_d(guide = "none") +
  labs(
    title = "B. Shoots per plate (Week 12; Neg. Bin. GLM)",
    x = "Shoot Induction Medium",
    y = "Shoots per plate (raw) + model mean ± 95% CI"
  ) +
  coord_cartesian(ylim = c(0, max(dat12$total_num_shoots, na.rm = TRUE) + 3)) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ---- 4) Combine + save ----
fig_shoot <- pA | pB
fig_shoot

ggsave("Shoot_multiplot_week12_v2.pdf", fig_shoot, width = 12, height = 6)

# ----------------------------
# Helper: consistent medium labels
# ----------------------------
medium_levels <- sort(unique(dat2$Medium))
dat2$Medium <- factor(dat2$Medium, levels = medium_levels)
dat_surv$Medium <- factor(dat_surv$Medium, levels = medium_levels)

# ----------------------------
# A) Rooting success panel (model-estimated probabilities)
# ----------------------------
emm_success <- emmeans(m_success, ~ Medium, type = "response")
succ_df <- as.data.frame(emm_success) %>%
  mutate(
    Medium = factor(Medium, levels = medium_levels),
    prob_pct = prob * 100,
    LCL_pct  = asymp.LCL * 100,
    UCL_pct  = asymp.UCL * 100
  )

# Letter groups (Tukey-adjusted)
succ_letters <- cld(emm_success, Letters = letters, adjust = "tukey") %>%
  as.data.frame() %>%
  mutate(
    Medium = factor(Medium, levels = medium_levels),
    group = gsub(" ", "", .group)
  ) %>%
  dplyr::select(Medium, group)

succ_df <- left_join(succ_df, succ_letters, by = "Medium")

pA <- ggplot(succ_df, aes(x = Medium, y = prob_pct, color = Medium)) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = LCL_pct, ymax = UCL_pct), width = 0.15) +
  geom_text(aes(label = group, y = UCL_pct + 4), size = 4, color = "black") +
  labs(
    title = "A. Rooting Success (GLMM)",
    x = "Medium",
    y = "Estimated Rooting (%) ± 95% CI"
  ) +
  coord_cartesian(ylim = c(0, max(succ_df$UCL_pct) + 12)) +
  scale_color_viridis_d() +
  theme_classic()

pA <- pA +
  scale_color_viridis_d(name = "Medium")

# ----------------------------
# B) Root count panel
# Raw distribution + model-estimated mean ±95% CI
# ----------------------------
emm_roots <- emmeans(m_roots, ~ Medium, type = "response")
roots_df <- as.data.frame(emm_roots) %>%
  mutate(
    Medium = factor(Medium, levels = medium_levels),
    mean_count = response,
    LCL = asymp.LCL,
    UCL = asymp.UCL
  )

roots_letters <- cld(emm_roots, Letters = letters, adjust = "tukey") %>%
  as.data.frame() %>%
  mutate(
    Medium = factor(Medium, levels = medium_levels),
    group = gsub(" ", "", .group)
  ) %>%
  dplyr::select(Medium, group)

roots_df <- left_join(roots_df, roots_letters, by = "Medium")

pB <- ggplot(dat2, aes(x = Medium, y = Root_Count)) +
  geom_violin(aes(fill = Medium), trim = FALSE, alpha = 0.85, show.legend = FALSE) +
  geom_jitter(aes(color = Medium), width = 0.12, height = 0, alpha = 0.35, size = 1.2, show.legend = FALSE) +
  geom_pointrange(
    data = roots_df,
    aes(x = Medium, y = mean_count, ymin = LCL, ymax = UCL, color = Medium),
    inherit.aes = FALSE,
    size = 0.7
  ) +
  geom_text(
    data = roots_df,
    aes(x = Medium, y = UCL + 1.2, label = group),
    inherit.aes = FALSE,
    color = "black",
    size = 4,
    vjust = -0.5
  ) +
  scale_fill_viridis_d(guide = "none") +
  scale_color_viridis_d(name = "Rooting Medium") +
  labs(
    title = "B. Number of Roots (Neg. Bin. GLMM)",
    x = "Medium",
    y = "Root Count (Raw) + Model Mean ± 95% CI"
  ) +
  coord_cartesian(ylim = c(0, max(dat2$Root_Count, na.rm = TRUE) + 4)) +
  theme_classic()


# ----------------------------
# C) Time-to-rooting panel (Kaplan–Meier curves)
# This is descriptive and aligns with your survival modeling.
# ----------------------------
km_fit <- survfit(Surv(time, event) ~ Medium, data = dat_surv)

pal <- palette(hcl.colors(6, "viridis"))


pC_obj <- ggsurvplot(
  km_fit,
  data = dat_surv,
  fun = "event",                
  censor = TRUE,
  palette = pal,
  conf.int = FALSE,
  risk.table = FALSE,
  xlim = c(0, 28),
  xlab = "Days Post Plating",
  ylab = "Probability Rooted",
  legend.title = "Medium",
  ggtheme = theme_classic()
)

pC <- pC_obj$plot + labs(title = "C. Time to rooting (Kaplan–Meier)")

# ----------------------------
# Combine into Figure 1
# ----------------------------
pA <- pA +
  theme(legend.position = "none")

pB <- pB +
  theme(legend.position = "right")

pC <- pC +
  theme(legend.position = "none")


fig1 <- (pA | pB) / pC

fig1

ggsave("Rooting_multiplot_fuggle_nodes.pdf", fig1, width = 11, height = 8)


# ============================================================
# Table 1 (model-based summaries)
# Columns:
# - n
# - Rooting prob (95% CI)
# - Mean root count (95% CI)
# - Hazard ratio vs Medium 1 (from coxme fixed effects)
# ============================================================

# n per medium
n_df <- dat2 %>%
  count(Medium, name = "n_explants")

# Rooting prob strings
succ_tbl <- succ_df %>%
  transmute(
    Medium,
    rooting_prob = prob,
    rooting_prob_CI = sprintf("%.3f (%.3f–%.3f)", prob, asymp.LCL, asymp.UCL)
  )

# Root count strings
roots_tbl <- roots_df %>%
  transmute(
    Medium,
    mean_roots = mean_count,
    mean_roots_CI = sprintf("%.2f (%.2f–%.2f)", mean_count, LCL, UCL)
  )

# Hazard ratios vs Medium 1 from coxme:
# Medium 1 is the reference; coxme returns coefficients for Medium2..6
cox_coef <- fixef(m_time)
cox_se   <- sqrt(diag(vcov(m_time)))

hr_df <- tibble(
  term = names(cox_coef),
  coef = as.numeric(cox_coef),
  se = as.numeric(cox_se)
) %>%
  mutate(
    Medium = factor(gsub("Medium", "", term), levels = levels(dat2$Medium)),
    HR = exp(coef),
    HR_LCL = exp(coef - 1.96 * se),
    HR_UCL = exp(coef + 1.96 * se),
    HR_CI = sprintf("%.3f (%.3f–%.3f)", HR, HR_LCL, HR_UCL)
  ) %>%
  dplyr::select(Medium, HR_CI)

# Add Medium 1 row (reference)
hr_df <- bind_rows(
  tibble(Medium = factor("1", levels = levels(dat2$Medium)),
         HR_CI = "1.000 (reference)"),
  hr_df
) %>%
  arrange(Medium)

# Assemble Table 1
table1 <- n_df %>%
  left_join(succ_tbl, by = "Medium") %>%
  left_join(roots_tbl, by = "Medium") %>%
  left_join(hr_df, by = "Medium") %>%
  mutate(
    Medium = as.character(Medium)
  ) %>%
  dplyr::select(
    Medium,
    n_explants,
    rooting_prob_CI,
    mean_roots_CI,
    HR_CI
  ) %>%
  rename(
    `Rooting probability (95% CI)` = rooting_prob_CI,
    `Mean root count (95% CI)` = mean_roots_CI,
    `Hazard ratio for rooting vs Medium 1 (95% CI)` = HR_CI
  )
