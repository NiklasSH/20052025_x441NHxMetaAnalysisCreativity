# --------- Load Libraries ---------
library(metafor)
library(tidyverse)
library(ggtext)
library(gridExtra)
library(broom)

# --------- Ensure Results directory exists ---------
output_dir <- "Results"
if (!dir.exists(output_dir)) dir.create(output_dir)

# --------- Load Datasets ---------
df_performance <- read_csv2("Human-AI_Creative_Performance.csv")
df_diversity <- read_csv2("Human-AI_Diversity.csv")
df_versus      <- read_csv2("Human_vs_AI.csv")

# --------- Mutate Date ---------
df_performance <- df_performance %>% 
  mutate(Date = as.Date(Date, format = ifelse(nchar(Date)==4, "%Y-%m-%d")))
df_diversity <- df_diversity %>% 
  mutate(Date = as.Date(Date, format = ifelse(nchar(Date)==4, "%Y-%m-%d")))
df_versus <- df_versus %>% 
  mutate(Date = as.Date(Date, format = ifelse(nchar(Date)==4, "%Y-%m-%d")))

# --------- Function to Compute Sampling Variance of d ---------
compute_vi <- function(cohens_d, n_treatment, n_control) {
  vi <- ((n_treatment + n_control) / (n_treatment * n_control)) + (cohens_d^2 / (2 * (n_treatment + n_control)))
  return(vi)
}

# --------- Wrapper to Compute vi If Not Present ---------
safe_calc_effects <- function(df) {
  if (!"vi" %in% names(df)) {
    df$vi <- compute_vi(df$cohens_d, df$n_treatment, df$n_control)
  }
  return(df)
}

# --------- Apply to All Datasets ---------
df_performance <- safe_calc_effects(df_performance)
df_diversity <- safe_calc_effects(df_diversity)
df_versus      <- safe_calc_effects(df_versus)

# --------- Meta-Analysis Function ---------
run_meta <- function(df, label, plot_filename, allow_moderator = TRUE, save_width, save_height) {
  cat("\n\n========== Meta-Analysis:", label, "==========\n\n")
  
  # 1) Fit random-effects model
  res         <- rma(yi = cohens_d, vi = vi, data = df, method = "REML")
  summary_res <- summary(res)
  print(summary_res)
  
  # 3) Prepare data frame for plotting
  #    -- compute per-study weight = 1/(vi + tau2)
  tau2   <- res$tau2
  weights <- 1 / (df$vi + tau2)
  
  plot_data <- df %>%
    mutate(
      effect = cohens_d,
      se     = sqrt(vi),
      ci.lb  = effect - 1.96 * se,
      ci.ub  = effect + 1.96 * se,
      weight = weights / sum(weights) * 100,  # percent weight
      study  = ID
    ) %>%
    arrange(desc(Date)) %>% 
    mutate(order = row_number())
  
  # 3b) classify each study’s significance & choose a colour
  plot_data <- plot_data %>%
    mutate(
      sign = case_when(
        ci.lb >  0          ~ "positive",
        ci.ub <  0          ~ "negative",
        TRUE                ~ "ns"
      ),
      col = case_when(
        sign == "positive"  ~ "green",
        sign == "negative"  ~ "red",
        sign == "ns"        ~ "grey40"
      )
    )
  
  # combine data for plotting
  plot_df <- plot_data   # no overall row
  
  # compute a little horizontal padding so text sits just outside the widest CI
  x_max <- max(plot_df$ci.ub, na.rm = TRUE)
  x_pad <- x_max + 0.05 * diff(range(plot_df$ci.lb, plot_df$ci.ub))
  
  # 4) Build ggplot forest plot
  dark_blue <- "#1f4e79"
  p <- ggplot(plot_df, aes(x = effect, y = order)) +
    
    # A) overall CI band (grey, 30% opacity)
    annotate("rect",
             xmin = summary_res$ci.lb,
             xmax = summary_res$ci.ub,
             ymin = -Inf, ymax = Inf,
             fill = "orange", alpha = 0.3) +
    # **new** black zero‐line
    geom_vline(xintercept = 0,
               color      = "black",
               size       = 0.7,
               alpha = 0.6) +
  
    # B) study‐level CI blocks (light‐blue, 30% opacity)
    geom_segment(aes(x    = ci.lb,
                     xend = ci.ub,
                     y    = order,
                     yend = order,
                     color = col),
                 size  = 7,
                 alpha = 0.3) +
    
    # 2) small dashed vertical at each studys own d
    geom_segment(
      aes(x    = effect,
          xend = effect,
          y    = order - 0.3,
          yend = order + 0.3,
          color = col),
      size     = 1
    ) +
    # 3) Weight squares
    geom_point(aes(size = weight, color = col),
               shape = 15) +
    scale_color_identity(
      name   = "Significance",
      guide  = "legend",
      labels = c(
        "green"  = "Significant positive",
        "red"    = "Significant negative",
        "grey40" = "Non‑significant"
      ),
      breaks = c("green", "red", "grey40")
    ) +
    # 4) Red line at summary d
    geom_vline(xintercept = summary_res$b[1],
               linetype = "dashed",
               color      = "orange",
               size       = 1) +
    # F) overall d label in the same boxed style as studies
    geom_label(
      data = tibble(
        effect = summary_res$b[1],
        order  = max(plot_df$order) + 0.3   # place just below the bottom study
      ),
      aes(
        x     = effect,
        y     = order,
        label = sprintf("d = %.3f", effect)
      ),
      fill       = "white",
      label.size = 0.2,
      color      = "orange",
      size       = 3.5,
      hjust      = 0.5
    ) +
    # 7) axes & theme tweaks
    scale_y_reverse(
      breaks = plot_df$order,
      labels = as.character(plot_df$Date),
      expand = expansion(add = c(0.5, 0.5)),     # gives 0.5 “row” of padding at top & bottom
      sec.axis = sec_axis(
        trans = ~ .,
        breaks = plot_df$order,
        labels = plot_df$study,
        name   = NULL    # no extra title on that side
      )
    ) +
    scale_size_continuous(
      range = c(3, 8),
      name  = "Study weight (%)"
    ) +
    # stack the two guides one on top of the other
    guides(
      color = guide_legend(nrow = 1, order = 1),  # importance: color first in the layout
      size  = guide_legend(nrow = 1, order = 2)
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position      = "bottom",
      legend.box           = "vertical",        # <-- stack legends
      plot.title.position  = "plot",            # <-- allow centering across full width
      plot.title           = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.y.left   = element_markdown(hjust = 0),
      axis.text.y.right  = element_markdown(hjust = 0),
      axis.ticks.y.right = element_blank()
    ) +
    labs(
      x     = bquote("Effect size (" * italic(d) * ")"),
      y     = NULL,
      title = paste("Forest Plot:", label)
    )
  
  # 5) Save & display
  ggsave(
    filename = file.path(output_dir, paste0(plot_filename, ".pdf")),
    plot     = p,
    width    = save_width,
    height   = save_height,
    dpi      = 600
    )
  
  # 5b) Violin plot of effect‐size heterogeneity
  p_violin <- ggplot(plot_df, aes(x = "", y = effect)) +
    geom_violin(fill    = "skyblue",
                alpha   = 0.6,
                size    = 0.5) +
    geom_hline(aes(yintercept = summary_res$b[1], linetype = "Meta‑analytic mean"),
               color     = "orange",
               size      = 0.8) +
    # ▶ annotate the median of the observed effects
    stat_summary(aes(shape = "Observed median"),
                 fun    = median,
                 geom   = "point",
                 size   = 3,
                 fill   = "white",
                 color  = "black") +
    geom_jitter(aes(color = col),
                width  = 0.1,
                size   = 4) +
    scale_color_identity() +
    scale_y_continuous() +
    scale_linetype_manual(name   = NULL,
                          values = c("Meta‑analytic mean" = "dashed")) +
    scale_shape_manual(name   = NULL,
                       values = c("Observed median"     = 23)) +
    labs(x     = NULL,
         y     = bquote("Effect size (" * italic(d) * ")"),
         title = paste("Effect‐Size Distribution:", label)) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.text.x        = element_blank(),
      axis.ticks.x       = element_blank(),
      plot.title         = element_text(size = 14, face = "bold"),
      legend.position    = "bottom",
      # stack the two keys horizontally
      legend.box         = "horizontal",
      legend.spacing.x   = unit(0.5, "cm")
    )
  
  ggsave(
    filename = file.path(output_dir, paste0(plot_filename, "_violin.pdf")),
    plot     = p_violin,
    width    = 6.8,
    height   = 4.5,
    dpi      = 600
  )
  
  print(p_violin)
  
  # 5c) Export summary table (and moderator tables) as PDF
  
  #  – main summary
  sum_df <- tibble::tibble(
    k      = res$k,
    Effect = round(summary_res$b, 3),
    SE     = round(summary_res$se,3),
    ci95   = sprintf("[%.3f,%.3f]", summary_res$ci.lb, summary_res$ci.ub),
    pval   = round(summary_res$pval, 3),
    Q      = round(res$QE,          2),
    df     = res$k - 1,
    pQ     = round(res$QEp,         3),
    i2     = round(summary_res$I2, 1),
    tau2   = round(res$tau2,       3)
  )
  
  res_grob <- tableGrob(
    sum_df, 
    rows = NULL,
    theme = ttheme_minimal(
      core = list(fg_params = list(hjust = 0, x = 0.1)),
      colhead = list(fg_params = list(fontface = "bold"))
    )
  )
  
  # ▶ funnel plot — observed
  pdf(file.path(output_dir, paste0(plot_filename, "_funnel_observed.pdf")), width = save_width, height = save_height)
  funnel(res, main = paste("Funnel plot —", label))
  dev.off()

  # ▶ funnel plot — trim-and-fill
  tf <- trimfill(res)
  pdf(file.path(output_dir, paste0(plot_filename, "_funnel_trimfill.pdf")), width = save_width, height = save_height)
  funnel(tf, main = paste("Trim‑and‑fill —", label), col = "red")
  dev.off()

  # 3) Moderators: GenAI_Type & Participants
  if(allow_moderator == TRUE){
    mods <- list()
    # GenAI_Type
    if(n_distinct(df$GenAI_Type)>1){
      mg <- rma(yi=cohens_d, vi=vi, mods=~factor(GenAI_Type), data=df, method="REML")
      mg_df <- tidy(mg, conf.int=TRUE) %>% filter(term!="(Intercept)")
      mods[["GenAI_Type"]] <- mg_df
      # table
      g_grob <- tableGrob(mg_df %>% select(term, estimate, std.error, conf.low, conf.high, p.value),
                          rows=NULL, theme=ttheme_minimal(core=list(fg_params=list(hjust=0))))
      ggsave(file.path(output_dir, paste0(plot_filename, "_mod_GenAI_Type_table.pdf")), g_grob, width=6.8, height=3)
    }
    # Participants
    if(n_distinct(df$Participants)>1){
      mp <- rma(yi=cohens_d, vi=vi, mods=~factor(Participants), data=df, method="REML")
      mp_df <- tidy(mp, conf.int=TRUE) %>% filter(term!="(Intercept)")
      mods[["Participants"]] <- mp_df
      # table
      p_grob <- tableGrob(mp_df %>% select(term, estimate, std.error, conf.low, conf.high, p.value),
                          rows=NULL, theme=ttheme_minimal(core=list(fg_params=list(hjust=0))))
      ggsave(file.path(output_dir, paste0(plot_filename, "_mod_Participants_table.pdf")), p_grob, width=6.8, height=3)
    }
    # combined moderator plot
    if(length(mods)>0){
      mod_all <- bind_rows(mods, .id="Moderator")
      p_mod <- ggplot(mod_all, aes(x=term, y=estimate, color=Moderator)) +
        geom_pointrange(aes(ymin=conf.low, ymax=conf.high), position=position_dodge(width=0.7)) +
        geom_hline(yintercept=0, linetype="dashed") + coord_flip() +
        labs(title=paste("Moderator effects —", label), x=NULL,
             y=bquote("Effect size (" * italic(d) * ")")) +
        theme_minimal()
      ggsave(file.path(output_dir, paste0(plot_filename, "_mod_plot.pdf")), p_mod, width=6.8, height=4.5)
    }
  }
  
  return(sum_df)
}
# --------- Run All Meta-Analyses ---------
enh_res <- run_meta(df_performance, "Human‑AI Creative Performance", "0_plot_performance", allow_moderator = TRUE,save_width  = 6.8, save_height = 5)
nov_res <- run_meta(df_diversity,  "AI Effect on Creative Diversity",  "0_plot_diversity", allow_moderator = FALSE,save_width  = 6.8, save_height = 3)
vs_res  <- run_meta(df_versus,       "Human versus AI",                   "0_plot_versus", allow_moderator = TRUE,save_width  = 6.8, save_height = 5)

# combined table
all_results <- bind_rows(
  enh_res  %>% mutate(Analysis="Performance"),
  nov_res  %>% mutate(Analysis="Diversity"),
  vs_res   %>% mutate(Analysis="Human_vs_AI")
) %>%
  mutate(
    Signif   = case_when(
      pval <  .001 ~ "***",
      pval <  .01  ~ "**",
      pval <  .05  ~ "*",
      TRUE         ~ ""
    ),
    `p (sig)` = paste0(format(pval, nsmall=3), Signif)
  )

all_grob <- tableGrob(
  all_results %>% 
    select(
      Analysis,
      Effect,
      SE,
      ci95,        
      `p (sig)`,
      Q,
      df,
      pQ,          
      i2,          
      tau2
    ),
  rows  = NULL,
  theme = ttheme_minimal(
    core    = list(fg_params = list(hjust=0, x=0.1, fontsize=10)),
    colhead  = list(fg_params = list(fontface="bold", fontsize=11))
  )
)

ggsave(
  filename = file.path(output_dir, "all_meta_analyses.pdf"),
  plot     = all_grob,
  width    = 6.8,
  height   = 2.5,
  dpi      = 600
  )