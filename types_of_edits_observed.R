
# Create a data frame of observed gene edits
edits_df <- data.frame(
  edit = c( 1, -1, 0, 0, -9, -8, -5, -5, -2, -2, -2, -2, -1,  1,  1,  1,  1,  1,  1,  1, 1, 0, -3, -3, -32),
  gene = c( "48594", "68957", "48594", "68957", "68957", "68957", "68957", "48594", "48594", "48594", "48594", "48594", "48594", "68957", "68957", "68957", "68957", "48594", "48594", "48594", "68957","68957", "50448", "50357", "68957"),
  line = c("Fuggle 2-1", "Fuggle 2-1", "Fuggle 2-1", "Fuggle 2-1", "Fuggle 6-1", "Fuggle 5-1", "Fuggle 5-1", "Fuggle 6-1", "Fuggle 5-1", "Fuggle 4-1", "Fuggle 3-1", "Tettnanger 1-1", "Fuggle 6-1", "Fuggle 6-1", "Fuggle 4-1", "Fuggle 1-4", "Fuggle 1-1", "Fuggle 5-1", "Fuggle 3-1", "Fuggle 1-1", "Tettnanger 1-1", "Fuggle 3-1","Fuggle G48-1", "Fuggle G48-1", "Fuggle 1-2"), stringsAsFactors = FALSE)

edits_df$edit <- factor(
  edits_df$edit,
  levels = sort(unique(edits_df$edit)),
  labels = ifelse(
    sort(unique(edits_df$edit)) > 0,
    paste0("+", sort(unique(edits_df$edit))),
    sort(unique(edits_df$edit))
  )
)

edits_df$edit <-as.factor(edits_df$edit)

edit_levels <- sort(unique(edits_df$edit))

# sum_prop
sum_prop <- edits_df %>%
  count(gene, edit, name = "n") %>%          # counts rows per gene x edit
  group_by(gene) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  complete(gene, edit = edit_levels, fill = list(n = 0, prop = 0)) %>%
  mutate(edit = factor(edit, levels = edit_levels))
  
# total n per gene
n_df <- sum_prop %>%
  group_by(gene) %>%
  summarise(n = sum(n), .groups = "drop")

ggplot(sum_prop, aes(x = gene, y = prop, fill = factor(edit))) +
  geom_col(color = "black", linewidth = 0.2) +
  coord_flip() +
  geom_text(
    data = n_df,
    aes(x = gene, y = 1.1, label = paste0("n = ", n)),
    inherit.aes = FALSE,
    size = 3.5
  ) +
scale_fill_manual(
  values = c(
    "-32" = "#08306B",   # darkest blue
    "-9"  = "#2B59C3",
    "-8"  = "#4F7FD9",
    "-5"  = "#86AEE8",
    "-3"  = "#BBD2F3",
    "-2"  = "#DEE9FA",
    "-1"  = "#EFF3F8",
    "0"   = "#F0F0F0",
    "+1"  = "#F4B6B6"
  ) 
) +
  labs(
    x = NULL,
    y = "Proportion of edits",
    fill = "Edit (bp)"
  ) +
  theme_classic()

ggsave("edit_types_barplot.svg",
       height = 4,
       width = 6,
       dpi=300)
