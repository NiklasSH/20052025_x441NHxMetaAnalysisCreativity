# --------- Load Libraries ---------
library(metafor)
library(tidyverse)
library(ggtext)
library(gridExtra)
library(broom)
library(here)

# --------- Ensure Results directory exists ---------
output_dir <- "Plots"
if (!dir.exists(output_dir)) dir.create(output_dir)

# --------- Load Datasets ---------
df_performance <- read_csv2( here("Data","Human-AI_Creative_Performance.csv") )
df_diversity   <- read_csv2( here("Data","Human-AI_Diversity.csv") )
df_versus      <- read_csv2( here("Data","Human_vs_AI.csv") )

# --------- Mutate Date ---------
df_performance <- df_performance %>% 
  mutate(Date = as.Date(Date, format = ifelse(nchar(Date)==4, "%Y-%m-%d")))
df_diversity <- df_diversity %>% 
  mutate(Date = as.Date(Date, format = ifelse(nchar(Date)==4, "%Y-%m-%d")))
df_versus <- df_versus %>% 
  mutate(Date = as.Date(Date, format = ifelse(nchar(Date)==4, "%Y-%m-%d")))

# --------- Computation of d ---------
compute_cohens_d <- function(df) {
  # 1. Means & SDs
  if (all(c("M_treatment", "M_control", "SD_treatment", "SD_control", "n_treatment", "n_control") %in% names(df))) {
    pooled_sd <- sqrt( ((df$n_treatment - 1) * df$SD_treatment^2 + (df$n_control - 1) * df$SD_control^2) / (df$n_treatment + df$n_control - 2) )
    df$cohens_d <- (df$M_treatment - df$M_control) / pooled_sd
    return(df)
  }
  # 2. F-value
  if (all(c("F_value", "n_treatment", "n_control") %in% names(df))) {
    df$cohens_d <- sqrt(df$F_value * (df$n_treatment + df$n_control) / (df$n_treatment * df$n_control))
    return(df)
  }
  # 3. Chi-square
  if (all(c("Chi_square", "n_total") %in% names(df))) {
    df$cohens_d <- sqrt(df$Chi_square / df$n_total)
    return(df)
  }
  # 4. Unstandardised beta
  if (all(c("unstandardised_beta", "SD_beta", "n_treatment", "n_control") %in% names(df))) {
    df$cohens_d <- df$unstandardised_beta / df$SD_beta
    return(df)
  }
  # 5. Mann–Whitney U (Z value)
  if (all(c("Z_value", "n_total") %in% names(df))) {
    r <- df$Z_value / sqrt(df$n_total)
    r <- pmin(pmax(r, -0.99), 0.99)  # prevent division by zero
    df$cohens_d <- (2 * r) / sqrt(1 - r^2)
    return(df)
  }
  # If none of the above, keep cohens_d as is (or NA)
  if (!"cohens_d" %in% names(df)) {
    df$cohens_d <- NA
  }
  return(df)
}

# --------- Function to Compute Sampling Variance of d ---------
compute_vi <- function(cohens_d, n_treatment, n_control) {
  vi <- ((n_treatment + n_control) / (n_treatment * n_control)) + (cohens_d^2 / (2 * (n_treatment + n_control)))
  return(vi)
}

# --------- Wrapper to Compute vi If Not Present ---------
safe_calc_effects <- function(df) {
  # If cohens_d is missing, try to compute it
  if (!"cohens_d" %in% names(df) || any(is.na(df$cohens_d))) {
    df <- compute_cohens_d(df)
  }
  # If vi is missing, compute it
  if (!"vi" %in% names(df) || any(is.na(df$vi))) {
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

  # 2) Prepare data frame for plotting
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
    arrange(desc(effect)) %>% 
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
        sign == "positive"  ~ "#006400",    # darkgreen
        sign == "negative"  ~ "#8B0000",    # darkred
        sign == "ns"        ~ "grey40"
      )
    )
  
  # combine data for plotting
  plot_df <- plot_data   # no overall row
  
  # compute a little horizontal padding so text sits just outside the widest CI
  x_max <- max(plot_df$ci.ub, na.rm = TRUE)
  x_pad <- x_max + 0.05 * diff(range(plot_df$ci.lb, plot_df$ci.ub))
  
  # 3) Build ggplot forest plot
  dark_blue <- "#1f4e79"
  p <- ggplot(plot_df, aes(x = effect, y = order)) +
    
    # 1) Overall CI band (orange, 30% opacity)
    annotate("rect",
             xmin = summary_res$ci.lb,
             xmax = summary_res$ci.ub,
             ymin = -Inf, ymax = Inf,
             fill = "orange", alpha = 0.3) +
    # 2) black zero‐line
    geom_vline(xintercept = 0,
               color      = "black",
               linewidth  = 0.7,
               alpha = 0.6) +
  
    # 3) Study‐level CI segments, coloured by significance (alpha = 0.3)
    geom_segment(aes(x    = ci.lb,
                     xend = ci.ub,
                     y    = order,
                     yend = order,
                     color = col),
                 size  = 4,
                 alpha = 0.3) +
    
    # 4) Small vertical lines at each study’s observed effect size
    geom_segment(
      aes(x    = effect,
          xend = effect,
          y    = order - 0.3,
          yend = order + 0.3,
          color = col),
      size     = 1
    ) +
    # 5) Weight squares (sized by study weight, coloured by significance)
    geom_point(aes(size = weight, color = col),
               shape = 15) +
    scale_color_identity(
      guide  = "legend",
      labels = c(
        "#006400" = "95% CI > 0",
        "#8B0000" = "95% CI < 0",
        "grey40"  = "95% CI = 0"
      ),
      breaks = c("#006400", "#8B0000", "grey40")
    ) +
    # 6) Summary effect (Cohen’s d) dashed line (orange)
    geom_vline(xintercept = summary_res$b[1],
               linetype   = "dashed",
               color      = "orange",
               linewidth  = 0.75) +
    # 7) Overall d label in the same boxed style as studies
    geom_label(
      data = tibble(
        effect = summary_res$b[1],
        order  = max(plot_df$order) + 0.3   # place just below the bottom study
      ),
      aes(
        x     = effect,
        y     = order,
        label = sprintf("Cohen's d = %.3f", effect)
      ),
      fill       = "white",
      label.size = 0.2,
      color      = "orange",
      size       = 3.5,
      hjust      = 0.5
    ) +
    # 8) axes & theme tweaks
    scale_y_reverse(
      breaks = plot_df$order,
      labels = as.character(plot_df$study),
      expand = expansion(add = c(0.5, 0.5)),     # gives 0.5 “row” of padding at top & bottom
      sec.axis = sec_axis(
        transform = ~ .,
        breaks    = plot_df$order,
        labels    = plot_df$study,
        name      = NULL    # no extra title on that side
      )
    ) +
    geom_text(aes(x = x_pad, y = order, label = sprintf("%.1f%%", weight)),
              hjust = 0, size = 4) +      # show weight in right column
    scale_size_continuous(
      range = c(3, 8),
      guide = "none"
    ) +
    # stack the two guides one on top of the other
    guides(
      color = guide_legend(nrow = 1, order = 1),  # importance: color first in the layout
    ) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position      = "bottom",
      legend.text          = element_text(size = 12),
      legend.title         = element_text(size = 12),
      legend.box           = "vertical",        # <-- stack legends
      axis.text            = element_text(size = 12),
      axis.text.y.left     = element_markdown(size = 12, hjust = 0),
      axis.text.y.right  = element_blank(),
      axis.ticks.y.right = element_blank(),
      plot.margin          = margin(5, 20, 5, 5, "pt")
    ) +
    coord_cartesian(clip = "off") +
    labs(
      x     = bquote("Effect size (Cohen's " * italic(d) * ")"),
      y     = NULL,
    )
  
  # 4) Save & display
  ggsave(
    filename = file.path(output_dir, paste0(plot_filename, "_forest.pdf")),
    plot     = p,
    width    = save_width,
    height   = save_height,
    dpi      = 600
    )
  
  # 5a_GenAI) Violin plot of effect‐size by GenAI_Type
  p_violin_genai <- ggplot(plot_df, aes(x = GenAI_Type, y = effect)) +
    geom_violin(fill    = "skyblue",
                alpha   = 0.6,
                size    = 0.5) +
    geom_hline(yintercept = summary_res$b[1],
               linetype    = "dashed",
               color       = "orange",
               linewidth   = 0.75,
               show.legend = FALSE) +
    # ▶ annotate the median of the observed effects
    stat_summary(fun     = median,
                 geom    = "point",
                 shape   = 21,        # circle with fill+border
                 size    = 4,         # outer diameter
                 stroke  = 1,         # border thickness
                 fill    = "white",   # inner color
                 color   = "black",
                 show.legend = FALSE) +
    geom_jitter(aes(color = col),
                width  = 0.1,
                size   = 4) +
    scale_color_identity() +
    scale_y_continuous() +
    labs(x     = NULL,
         y     = bquote("Effect size (Cohen's " * italic(d) * ")")) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title           = element_text(size = 12),
      axis.text            = element_text(size = 12),
      axis.text.x          = element_text(size = 12, angle = 30, hjust = 1),
      panel.grid.major     = element_line(color = "grey80"),
      panel.grid.minor     = element_blank(),
      plot.title           = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.margin          = margin(5, 5, 5, 5)
    )
  
  ggsave(
    filename = file.path(output_dir, paste0(plot_filename, "_violin_GenAI.pdf")),
    plot     = p_violin_genai,
    width    = 6.8,
    height   = 4.5,
    dpi      = 600
  )
  
  print(p_violin_genai)
  
  # 5b_Part) Violin plot of effect‐size by Participants
  p_violin_part <- ggplot(plot_df, aes(x = Participants, y = effect)) +
    geom_violin(fill    = "skyblue",
                alpha   = 0.6,
                size    = 0.5) +
    geom_hline(yintercept = summary_res$b[1],
               linetype    = "dashed",
               color       = "orange",
               linewidth   = 0.75,
               show.legend = FALSE) +
    # ▶ annotate the median of the observed effects
    stat_summary(fun     = median,
                 geom    = "point",
                 shape   = 21,        # circle with fill+border
                 size    = 4,         # outer diameter
                 stroke  = 1,         # border thickness
                 fill    = "white",   # inner color
                 color   = "black",
                 show.legend = FALSE) +
    geom_jitter(aes(color = col),
                width  = 0.1,
                size   = 4) +
    scale_color_identity() +
    scale_y_continuous() +
    labs(x     = NULL,
         y     = bquote("Effect size (Cohen's " * italic(d) * ")")) +
    theme_minimal(base_size = 12) +
    theme(
      axis.title           = element_text(size = 12),
      axis.text            = element_text(size = 12),
      axis.text.x          = element_text(size = 12, angle = 30, hjust = 1),
      panel.grid.major     = element_line(color = "grey80"),
      panel.grid.minor     = element_blank(),
      plot.title           = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.margin          = margin(5, 5, 5, 5)
    )
  
  # Save & show
  ggsave(
    filename = file.path(output_dir, paste0(plot_filename, "_violin_participants.pdf")),
    plot     = p_violin_part,
    width    = 6.8,
    height   = 4.5,
    dpi      = 600
  )
  print(p_violin_part)
  
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
  
  # 5d) Funnel plot — observed
  pdf(file.path(output_dir, paste0(plot_filename, "_funnel_observed.pdf")), width = 6.8, height = 5)
  funnel(res, main = NULL)
  dev.off()

  # 5e) Funnel plot — trim-and-fill
  tf <- trimfill(res)
  pdf(file.path(output_dir, paste0(plot_filename, "_funnel_trimfill.pdf")), width = 6.8, height = 5)
  funnel(tf, main = NULL, col = "#8B0000")
  dev.off()

  # 6) Moderators: GenAI_Type & Participants
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
        labs(x=NULL,
             y=bquote("Effect size (" * italic(d) * ")")) +
        theme_minimal(base_size = 12) +
        theme(
          legend.position      = "right",
          legend.text          = element_text(size = 12),
          legend.title         = element_text(size = 12),
          legend.box           = "vertical",
          axis.title           = element_text(size = 12),
          axis.text            = element_text(size = 12),
          panel.grid.major     = element_line(color = "grey80"),
          panel.grid.minor     = element_blank(),
          plot.title           = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.margin          = margin(5, 5, 5, 5)
        )
      ggsave(file.path(output_dir, paste0(plot_filename, "_mod_plot.pdf")), p_mod, width=6.8, height=4.5)
    }
  }
  
  return(sum_df)
}
# --------- Run All Meta-Analyses ---------
enh_res <- run_meta(df_performance, "Human‑AI Creative Performance", "plot_performance", allow_moderator = TRUE,save_width  = 6.8, save_height = 5)
nov_res <- run_meta(df_diversity,  "AI Effect on Creative Diversity",  "plot_diversity", allow_moderator = FALSE,save_width  = 6.8, save_height = 3)
vs_res  <- run_meta(df_versus,       "Human versus AI",                   "plot_versus", allow_moderator = TRUE,save_width  = 6.8, save_height = 5)

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