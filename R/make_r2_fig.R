# create a figure examining cv R^2 for each bnab using ic50 and ic80 outcomes
# similar to figure slapnap_supplemental/single_sensitivity/R/make_auc_fig.R

# --------------------------------------------------------------------------
# Set up
# --------------------------------------------------------------------------
# load required libraries
library("dplyr")
library("forcats")
library("tibble")
library("ggplot2")
library("cowplot")
theme_set(theme_cowplot())
library("here")

# create a tibble with the data
#load(here("py_output", "cvr2_df.cvs"))
rslt_df <- read.csv(file = "./R_output/new_cvr2_df.csv")
#rslt_df <- rslt_df[complete.cases(rslt_df), ] #remove any NA rows


compare_tib <- rslt_df %>%
  as_tibble() %>%
  mutate(
    bnab = toupper(bnab),
    epitope = case_when(
      bnab %in% c("PG16", "PG9", "PGDM1400",
                  "PGT145", "VRC26.08", "VRC26.25") ~ "V1V2",
      bnab %in% c("10-1074", "10-996", "DH270.1",
                  "DH270.5", "DH270.6", "PGT121", "PGT128",
                  "PGT135", "VRC29.03", "VRC38.01", "2G12") ~ "V3",
      bnab %in% c("3BNC117", "B12", "CH01",
                  "HJ16", "NIH45-46", "VRC-CH31", "VRC-PG04",
                  "VRC01", "VRC03", "VRC07") ~ "CD4bs",
      bnab %in% c("PGT151", "VRC34.01") ~ "Fusion peptide",
      bnab %in% c("35O22", "8ANC195") ~ "Subunit interface",
      bnab %in% c("2F5", "4E10") ~ "MPER"
    ))


# --------------------------------------------------------------------------
# Create the plot
# --------------------------------------------------------------------------
epitope_labs <- unique(compare_tib$epitope)[order(unique(compare_tib$epitope), decreasing = TRUE)]
point_size <- .5
font_size <- 5.5
line_size <- 0.25
ci_line_size <- .3
ci_line_width <- .8
continuous_scale <- c(.2, 6/5)

cd4bs_plot <- compare_tib %>%
  filter(epitope == "CD4bs") %>%
  ggplot(aes(x = forcats::fct_reorder(as.factor(bnab), desc(epitope)),
             y = r2, ymin = cil, ymax = ciu, color = Method,
             group = paste0(epitope, "_", Method))) +
  geom_point(position = position_dodge(width = 0.75, preserve = "total"),
             aes(size = n), show.legend = FALSE) +
  scale_size_continuous(range = continuous_scale) + 
  #geom_text(aes(label = n, hjust = 1.5)) + 
  geom_errorbar(aes(ymin = cil, ymax = ciu), size = ci_line_size, width = .5,
                position = position_dodge(width = 0.75, preserve = "total")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = line_size) +
  ylim(c(-.5, 1)) + #may need to update
  ylab(expression(paste("CV-", R^2))) + 
  ggtitle("CD4bs") +
  xlab("") +
  #labs(y = NULL) +
  guides(x = guide_axis(n.dodge = 2), shape = FALSE, color = FALSE) +
  theme(plot.title = element_text(hjust = 0.5, size = font_size),
        axis.title.x = element_blank(),
        text = element_text(size = font_size),
        axis.text = element_text(size = font_size * .75),
        axis.line = element_line(size = line_size),
        axis.ticks = element_line(size = line_size),
        plot.margin = margin(0, 0.1, 0.2, 0.1, "cm")) + 
  background_grid(size.major = line_size, size.minor = line_size) +
  scale_color_manual(values = c("#808080", "#A9A9A9"))

fusion_plot <- compare_tib %>%
  filter(epitope == "Fusion peptide") %>%
  ggplot(aes(x = forcats::fct_reorder(as.factor(bnab), desc(epitope)),
             y = r2, ymin = cil, ymax = ciu, color = Method,
             group = paste0(epitope, "_", Method))) +
  geom_point(position = position_dodge(width = 0.75, preserve = "total"),
             aes(size = n), show.legend = FALSE) +
  scale_size_continuous(range = continuous_scale) + 
  geom_errorbar(aes(ymin = cil, ymax = ciu), size = ci_line_size, width = .3,
                position = position_dodge(width = 0.75, preserve = "total")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = line_size) +
  ylim(c(-.5, 1)) + #may need to update
  labs(y = NULL) +
  ggtitle("Fusion peptide") +
  xlab("") +
  guides(x = guide_axis(n.dodge = 2), shape = FALSE, y = "none", color = FALSE) +
  theme(plot.title = element_text(hjust = 0.5, size = font_size),
        axis.title.x = element_blank(),
        text = element_text(size = font_size),
        axis.text = element_text(size = font_size * .75),
        axis.line = element_line(size = line_size),
        axis.ticks = element_line(size = line_size),
        plot.margin = margin(0, 0.1, 0.2, 0, "cm")) + 
  background_grid(size.major = line_size, size.minor = line_size) +
  scale_color_manual(values = c("#808080", "#A9A9A9"))

mper_plot <- compare_tib %>%
  filter(epitope == "MPER") %>%
  ggplot(aes(x = forcats::fct_reorder(as.factor(bnab), desc(epitope)),
             y = r2, ymin = cil, ymax = ciu, color = Method,
             group = paste0(epitope, "_", Method))) +
  geom_point(position = position_dodge(width = 0.75, preserve = "total"),
             aes(size = n)) +
  scale_size_continuous(range = continuous_scale) + 
  geom_errorbar(aes(ymin = cil, ymax = ciu), size = ci_line_size, width = .3,
                position = position_dodge(width = 0.75, preserve = "total")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = line_size) +
  ylim(c(-.5, 1)) + #may need to update
  labs(y = NULL) +
  ggtitle("MPER") +
  xlab("") +
  guides(x = guide_axis(n.dodge = 2), y = "none") +
  theme(plot.title = element_text(hjust = 0.5, size = font_size),
        axis.title.x = element_blank(),
        text = element_text(size = font_size),
        legend.title = element_text(size = font_size),
        legend.text = element_text(size = font_size * .75),
        legend.margin = margin(0, 0, 0.1, -0.15, "cm"),
        legend.box.margin = margin(0, 0, 0.1, -0.15, "cm"),
        legend.spacing.y = unit(.001, 'cm'),
        legend.key.size = unit(.5, "lines"),
        axis.text = element_text(size = font_size * .75),
        axis.line = element_line(size = line_size),
        axis.ticks = element_line(size = line_size),
        plot.margin = margin(0, 0.15, 0.2, 0, "cm")) + 
  #guides(shape = guide_legend(override.aes = list(size = 0.5))) + 
  #guides(color = guide_legend(override.aes = list(size = 0.5))) + 
  #theme(legend.title = element_text(size = 3), 
   #     legend.text = element_text(size = 3)) + 
  background_grid(size.major = line_size, size.minor = line_size) +
  scale_color_manual(values = c("#808080", "#A9A9A9"))

subunit_plot <- compare_tib %>%
  filter(epitope == "Subunit interface") %>%
  ggplot(aes(x = forcats::fct_reorder(as.factor(bnab), desc(epitope)),
             y = r2, ymin = cil, ymax = ciu, color = Method,
             group = paste0(epitope, "_", Method))) +
  geom_point(position = position_dodge(width = 0.75, preserve = "total"),
             aes(size = n), show.legend = FALSE) +
  scale_size_continuous(range = continuous_scale) +
  geom_errorbar(aes(ymin = cil, ymax = ciu), size = ci_line_size, width = .25,
                position = position_dodge(width = 0.75, preserve = "total")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = line_size) +
  ylim(c(-.5, 1)) +
  labs(y = NULL) +
  ggtitle("Subunit interface") +
  xlab("") +
  guides(x = guide_axis(n.dodge = 2), shape = FALSE, y = "none", color = FALSE) +
  theme(plot.title = element_text(hjust = 0.5, size = font_size),
        axis.title.x = element_blank(),
        text = element_text(size = font_size),
        axis.text = element_text(size = font_size * .75),
        axis.line = element_line(size = line_size),
        axis.ticks = element_line(size = line_size),
        plot.margin = margin(0, 0.1, 0.2, 0, "cm")) + 
  background_grid(size.major = line_size, size.minor = line_size) +
  scale_color_manual(values = c("#808080", "#A9A9A9"))

v1v2_plot <- compare_tib %>%
  filter(epitope == "V1V2") %>%
  ggplot(aes(x = forcats::fct_reorder(as.factor(bnab), desc(epitope)),
             y = r2, ymin = cil, ymax = ciu, color = Method,
             group = paste0(epitope, "_", Method))) +
  geom_point(position = position_dodge(width = 0.75, preserve = "total"),
             aes(size = n), show.legend = FALSE) +
  scale_size_continuous(range = continuous_scale) +
  geom_errorbar(aes(ymin = cil, ymax = ciu), size = ci_line_size, width = .4,
                position = position_dodge(width = 0.75, preserve = "total")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = line_size) +
  ylim(c(-.5, 1)) +
  ylab(expression(paste("CV-", R^2))) +
  ggtitle("V1V2") +
  xlab("") +
  guides(x = guide_axis(n.dodge = 2), shape = FALSE, color = FALSE) +
  theme(plot.title = element_text(hjust = 0.5, size = font_size),
        axis.title.x = element_blank(),
        text = element_text(size = font_size),
        axis.text = element_text(size = font_size * .75),
        axis.line = element_line(size = line_size),
        axis.ticks = element_line(size = line_size),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm")) + 
  background_grid(size.major = line_size, size.minor = line_size) +
  scale_color_manual(values = c("#808080", "#A9A9A9"))

v3_plot <- compare_tib %>%
  filter(epitope == "V3") %>%
  ggplot(aes(x = forcats::fct_reorder(as.factor(bnab), desc(epitope)),
             y = r2, ymin = cil, ymax = ciu, color = Method,
             group = paste0(epitope, "_", Method))) +
  geom_point(position = position_dodge(width = 0.75, preserve = "total"),
             aes(size = n), show.legend = FALSE) + 
  scale_size_continuous(range = continuous_scale) +
  geom_errorbar(aes(ymin = cil, ymax = ciu), size = ci_line_size, width = .5,
                position = position_dodge(width = 0.75, preserve = "total")) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "red", size = line_size) +
  ylim(c(-.5, 1)) +
  labs(y = NULL) +
  ggtitle("V3") +
  xlab("") +
  guides(x = guide_axis(n.dodge = 2), y = "none", shape = FALSE, color = FALSE) +
  theme(plot.title = element_text(hjust = 0.5, size = font_size),
        axis.title.x = element_blank(),
        text = element_text(size = font_size),
        axis.text = element_text(size = font_size * .75),
        axis.line = element_line(size = line_size),
        axis.ticks = element_line(size = line_size),
        plot.margin = margin(0.1, 0.1, 0.1, 0, "cm")) + 
  background_grid(size.major = line_size, size.minor = line_size) +
  scale_color_manual(values = c("#808080", "#A9A9A9"))


compare_plot <- plot_grid(
  plot_grid(
    plot_grid(v1v2_plot, v3_plot, nrow = 1),
    plot_grid(cd4bs_plot, fusion_plot, subunit_plot, mper_plot, nrow = 1,
              rel_widths = c(3, 1, 1, 1.5)), #may need to change
              #rel_widths = c(1.5, 0.5, 0.5, 1)),
    nrow = 2
  ),
  ggplot() + ggtitle("bnAb") + guides(x = "none", y = "none") +
    theme(plot.title = element_text(hjust = 0.5, face = "plain", size = font_size)),
  nrow = 2, rel_heights = c(1, 0.02)
)


fig_width <- 5
fig_height <- 2.25
ggsave(filename = "./R_output/new_r2_fig.png",
       plot = compare_plot, 
       width = fig_width, height = fig_height, 
       units = "in", dpi = 300)
