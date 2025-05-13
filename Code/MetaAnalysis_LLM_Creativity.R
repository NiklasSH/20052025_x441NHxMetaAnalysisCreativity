# --------- Load Libraries ---------
library(metafor)
library(tidyverse)
library(ggtext)
library(gridExtra)
library(broom)
library(here)
library(knitr)
library(xtable)

# --------- Ensure Results directory exists ---------
output_dir <- "Plots"
if (!dir.exists(output_dir)) dir.create(output_dir)

# --------- Load Datasets ---------
df_performance <- read_csv2( here("Data","Human-AI_Creative_Performance.csv") )
df_diversity   <- read_csv2( here("Data","Human-AI_Diversity.csv") )
df_versus      <- read_csv2( here("Data","Human_vs_AI.csv") )

# --------- Computation of d ---------
compute_cohens_d <- function(df) {
  # make sure every "possible" column exists
  all.possible <- c(
    "M_treatment","SD_treatment","SE_treatment","M_control","SD_control","SE_control","n_treatment","n_control",
    "F_value","Chi_square","unstandardised_beta","SE_beta","Z_value","n_total",
    "cohens_d", "t_value"  # in case the author already supplied it
  )
  missing.cols <- setdiff(all.possible, names(df))
  for(col in missing.cols) df[[col]] <- NA_real_
  
  df %>%
    mutate(
      cohens_d = if_else(
        !is.na(cohens_d),   # 1) keep author‐reported d
        cohens_d,
        case_when(
          # 2) Means & SDs
          !is.na(M_treatment) & !is.na(M_control) &
            !is.na(SD_treatment) & !is.na(SD_control) ~
            (M_treatment - M_control) /
            sqrt(((n_treatment - 1)*SD_treatment^2 +
                    (n_control  - 1)*SD_control^2) /
                   (n_treatment + n_control - 2)),
          
          # 3) Means & SEs
          !is.na(M_treatment) & !is.na(M_control) &
            !is.na(SE_treatment) & !is.na(SE_control) ~
            (M_treatment - M_control) /
            sqrt(((n_treatment - 1)*(SE_treatment*sqrt(n_treatment))^2 +
                    (n_control  - 1)*(SE_control*sqrt(n_control))^2) /
                   (n_treatment + n_control - 2)),
          
          # 4) F‐value
          !is.na(F_value) ~
            sqrt(F_value) * sqrt(1/n_treatment + 1/n_control),
          
          # 5) Unstandardised beta
          !is.na(unstandardised_beta) & !is.na(SE_beta) ~
            (unstandardised_beta / SE_beta)*sqrt(1/n_control + 1/n_treatment),
          
          # 5) t_value
          !is.na(t_value) ~
            t_value*sqrt(1/n_control + 1/n_treatment),
          
          TRUE ~ NA_real_
        )
      )
    )
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
  # — Compute Hedges’ g correction factor J and its variance
  df <- df %>%
    mutate(
      # degrees of freedom for correction
      df_total  = n_treatment + n_control - 2,
      J         = 1 - (3 / (4 * df_total - 1)),
      hedges_g  = cohens_d * J,
      vi_g      = vi * J^2
      )
  return(df)
}

# --------- Apply to All Datasets ---------
df_performance <- safe_calc_effects(df_performance)
df_diversity <- safe_calc_effects(df_diversity)
df_versus      <- safe_calc_effects(df_versus)

# --------- Inspect all computed effect‐sizes ---------
inspect_effects <- bind_rows(
  df_performance %>%
    select(ID, Report_ID, cohens_d, hedges_g) %>%
    mutate(dataset = "performance"),
  df_diversity %>%
    select(ID, Report_ID, cohens_d, hedges_g) %>%
    mutate(dataset = "diversity"),
  df_versus %>%
    select(ID, Report_ID, cohens_d, hedges_g) %>%
    mutate(dataset = "versus")
  )

# optionally write to CSV for closer inspection
  write_csv(inspect_effects, file.path(output_dir, "all_effect_sizes.csv"))

# --------- Meta-Analysis Function ---------
run_meta <- function(df, label, plot_filename, allow_moderator = TRUE) {
  cat("\n\n========== Meta-Analysis:", label, "==========\n\n")
  
  # 1) Fit random-effects model
  res         <- rma(yi = hedges_g, vi = vi_g, data = df, method = "REML")
  summary_res <- summary(res)
  print(summary_res)
  
  # — Change #1: Influence diagnostics on raw data only
  if (allow_moderator) {
    inf <- influence(res)
    print(inf)
    # — Change #4a: save influence plot to PDF
    pdf(
      file.path(output_dir,
                paste0(plot_filename, "_influence_diagnostics.pdf")),
      width  = 6.5,
      height = 6.5
      )
    plot(inf)
    dev.off()
  }
  
  # — Change #2: Leave-one-out sensitivity analysis
  sens <- leave1out(res)
  print(sens)
  est <- sens$estimate  
  pdf(
    file.path(output_dir,
                  paste0(plot_filename, "_leave1out_sensitivity.pdf")),
    width  = 6.5,
    height = 6.5
    )  
  plot(est,  
       type = "b",  
       xlab = "Omitted study",  
       ylab = "Pooled estimate")  
  dev.off()# basic leave-one-out plot

  # — Change #3: Publication-bias checks only if k ≥ 10
  if (res$k >= 10) {
    egger_test <- regtest(res, model = "rma")
    print(egger_test)
    
    # raw funnel
    pdf(file.path(output_dir, paste0(plot_filename, "_funnel_plot.pdf")),
        width = 6.5, height = 6.5)
    funnel(res)
    dev.off()
    
    # trimmed-and-filled funnel
    tf <- trimfill(res)
    print(summary(tf))
    orig_k <- res$k
    pdf(file.path(output_dir, paste0(plot_filename, "_funnel_trimfill.pdf")),
        width = 6.5, height = 6.5)
    funnel(res, pch = 19, col = "black")
    if (tf$k0 > 0) {
      pts <- (orig_k + 1):(orig_k + tf$k0)
      points(
        x   = tf$yi[pts],
        y   = sqrt(tf$vi[pts]),
        pch = 21, bg = "red", col = "red"
      )
    }
    dev.off()
  } else {
    message("Skipping Egger’s test and funnel/trim‐and‐fill (k = ", res$k, " < 10).")
  }
  
  # 2) Prepare data frame for plotting
  #    -- compute per-study weight = 1/(vi + tau2)
  tau2   <- res$tau2
  weights <- 1 / (df$vi_g + tau2)
  
  plot_data <- df %>%
    mutate(
      effect = hedges_g,
      se     = sqrt(vi_g),
      ci.lb  = effect - 1.96 * se,
      ci.ub  = effect + 1.96 * se,
      weight = weights / sum(weights) * 100,  # percent weight
      study  = paste0(ID, "_", Report_ID)) %>%
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
  
  # Compute dynamic height
  n_rows        <- nrow(plot_df)
  row_height    <- 0.25    # cm per row
  extra_space   <- 3      # cm for margins
  plot_height_cm <- n_rows * row_height + extra_space
  
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
      linewidth = 1
    ) +
    scale_color_identity(
      guide  = "legend",
      labels = c(
        "#006400" = "95% CI > 0",
        "#8B0000" = "95% CI < 0",
        "grey40"  = "95% CI = 0"
      ),
      breaks = c("#006400", "#8B0000", "grey40")
    ) +
    # 6) Summary effect (Hedges’ g) dashed line (orange)
    geom_vline(xintercept = summary_res$b[1],
               linetype   = "dashed",
               color      = "orange",
               linewidth  = 0.75) +
    # 7) Overall g label in the same boxed style as studies
    geom_label(
      data = tibble(
        effect = summary_res$b[1],
        order  = max(plot_df$order) + 0.3   # place just below the bottom study
      ),
      aes(
        x     = effect,
        y     = order,
        label = sprintf("Hedges' g = %.3f", effect)
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
      expand = expansion(add = c(0.3, 0.3)),     # gives 0.5 “row” of padding at top & bottom
      sec.axis = sec_axis(
        transform = ~ .,
        breaks    = plot_df$order,
        labels    = plot_df$study,
        name      = NULL    # no extra title on that side
      )
    ) +
    geom_text(aes(x = x_pad, y = order, label = sprintf("%.1f%%", weight)),
              hjust = 0, size = 4) + # show weight in right column
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
      axis.text.y.left     = element_markdown(size = 12),
      axis.text.y.right    = element_blank(),
      axis.ticks.y.right   = element_blank(),
      plot.margin          = margin(5,5,5,5, "mm")
    ) +
    coord_cartesian(clip = "off") +
    labs(
      x     = bquote("Effect size (Hedges' " * italic(g) * ")"),
      y     = NULL,
    )
  
  # 4a) Save & display
  ggsave(
    filename = file.path(output_dir, paste0(plot_filename, "_forest.pdf")),
    plot     = p,
    dpi      = 600,
    width    = 6.5, 
    height   = plot_height_cm,
    units    = "in"
  )
  # 0) If this is the diversity‐dataset, bail out immediately
  if (identical(deparse(substitute(df)), "df_diversity")) {
    message("Skipping violin‐plot generation for df_diversity.")
  } else {
    
    # 1) Define for each moderator which values to drop
    exclude_levels <- list(
      GenAI_Model = c("2 Models", "3 models", "5 models", "> 5 models", "n.d."),
      Participants = "Not disclosed",
      Task_Type   = "ideation other"
    )
    
    # 2) Identify which of your static mods are actually present & non‐NA
    static_mods   <- c("GenAI_Model", "Participants", "Task_Type")
    present_static <- intersect(static_mods, names(df))
    non_na_static <- present_static[
      vapply(df[present_static], function(col) any(!is.na(col)), logical(1))
    ]
    moderators <- non_na_static
    
    # 3) Loop over each moderator
    for (mod in moderators) {
      # start from the full plot_df each time
      plot_data_mod <- plot_df
      
      # 3a) drop the unwanted levels for *this* moderator, if any
      if (mod %in% names(exclude_levels)) {
        plot_data_mod <- plot_data_mod %>%
          filter(!.data[[mod]] %in% exclude_levels[[mod]])
      }
      
      # 3b) now drop any levels with fewer than 2 obs
      lvl_counts  <- table(plot_data_mod[[mod]], useNA="no")
      keep_levels <- names(lvl_counts)[lvl_counts >= 2]
      plot_data_mod <- plot_data_mod %>%
        filter(.data[[mod]] %in% keep_levels) %>%
        droplevels()
      
      # --- decide whether this plot should get a y‑axis label ----
      y_lab <- if (mod == "GenAI_Model") {
        bquote("Effect size (Hedges' " * italic(g) * ")")
      } else {
        NULL        # no y‑label for the other moderators
      }
      
      # 3c) build & save your violin
      p_raw <- ggplot(plot_data_mod, aes_string(x = mod, y = "effect")) +
        geom_violin(
          fill      = "skyblue",
          width     = 0.8,
          scale     = "width",
          alpha     = 0.6,
          linewidth = 0.5,
          trim      = FALSE,
          adjust    = 1.5
        ) +
        geom_hline(
          yintercept = summary_res$b[1],
          linetype   = "dashed",
          color      = "orange",
          linewidth  = 0.75
        ) +
        geom_jitter(
          data   = plot_data_mod,
          aes(color = col),
          width  = 0.1,
          size   = 2
        ) +
        stat_summary(
          fun    = median,
          geom   = "point",
          shape  = 21,
          size   = 2,
          stroke = 1,
          fill   = "white"
        ) +
        scale_color_identity() +
        labs(
          x = NULL,
          y = y_lab
        ) +
        theme_minimal(base_size = 12) +
        theme(
          axis.text.x     = element_text(angle = 30, hjust = 1),
          panel.grid.minor= element_blank()
        )
      
      ggsave(
        filename = file.path(output_dir, paste0(plot_filename, "_violin_", mod, ".pdf")),
        plot     = p_raw,
        dpi      = 600,
        width    = 6.5
      )
    }
  }
  
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
  
  # make a reusable table‐theme with larger text & padding
  tt <- ttheme_minimal(
    base_size = 12,                       # default font size
    padding   = unit(c(4, 4), "mm"),      # cell padding goes here, not inside fg_params
    core      = list(
      fg_params = list(
        fontsize = 12,                    # core text size
        hjust    = 0                      # left‐align
      )
    ),
    colhead   = list(
      fg_params = list(
        fontsize = 14,                    # header text size
        fontface = "bold"
      )
    )
  )
  
  # 6) Moderators
  if (allow_moderator) {
    # reusable big-font table theme
    tt <- ttheme_minimal(
      base_size = 12,
      padding   = unit(c(4,4), "mm"),
      core      = list(fg_params = list(fontsize = 12, hjust = 0)),
      colhead   = list(fg_params = list(fontsize = 14, fontface = "bold"))
    )
    
    # which moderator variables to test
    mod_vars <- c(
      "GenAI_Type", "Participants", "Recruitment_Source", "Task_Type",
      "GenAI_Model", "Creativity_Measurement", "Measurement_Evaluator"
    )
    
    for (mod in mod_vars) {
      # only if it exists in df and has more than one level:
      if (mod %in% names(df) && n_distinct(df[[mod]], na.rm = TRUE) > 1) {
        # — drop any levels with <2 observations (use same data as violin)
        lvl_counts  <- table(plot_df[[mod]], useNA="no")
        keep_levels <- names(lvl_counts)[lvl_counts >= 2]
        # if fewer than 2 levels remain, skip this moderator altogether
        if (length(keep_levels) < 2) {
          message("Skipping moderator ", mod, ": only ", length(keep_levels),
                  " level(s) after filtering.")
          next
        }
        df_mod      <- df %>% filter(.data[[mod]] %in% keep_levels)
        # 1) fit meta-regression
        form    <- as.formula(paste0("~ 0 + factor(", mod, ")"))
        fit     <- rma(yi = hedges_g, vi = vi_g, mods = form, data = df_mod, method = "REML")
        fit_df  <- tidy(fit, conf.int = TRUE) %>% filter(term != "(Intercept)")
        # 3) make & save plot
        p_mod <- ggplot(fit_df, aes(x = term, y = estimate)) +
          geom_pointrange(aes(ymin = conf.low, ymax = conf.high), size = 0.6) +
          geom_hline(yintercept = 0, linetype = "dashed") +
          coord_flip() +
          labs(
            x = NULL,
            y = bquote("Hedges’" ~ italic(g) ~ "effect")
          ) +
          theme_minimal(base_size = 12) +
          theme(
            axis.title   = element_text(size = 13),
            axis.text    = element_text(size = 11),
            plot.margin  = margin(5,5,5,5,"mm")
          )
        # dynamic plot sizing:
        n_rows     <- nrow(fit_df)
        row_h_cm   <- 0.5             # cm per row
        extra_cm   <- 3               # cm top+bottom
        plot_h_cm  <- n_rows * row_h_cm + extra_cm
        ggsave(
          file.path(output_dir, paste0(plot_filename, "_mod_", mod, "_plot.pdf")),
          p_mod,
          width  = 6.5, 
          height = plot_h_cm / 2.54,
          units  = "in",
          dpi    = 600,
        )
      }
    }
  }
  return(sum_df)
}
# --------- Run All Meta-Analyses ---------
enh_res <- run_meta(df_performance, "Human‑AI Creative Performance", "plot_performance_raw", allow_moderator = TRUE)
nov_res <- run_meta(df_diversity,  "AI Effect on Creative Diversity",  "plot_diversity_raw", allow_moderator = TRUE)
vs_res  <- run_meta(df_versus,       "Human versus AI",                   "plot_versus_raw", allow_moderator = TRUE)

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

raw_tab <- all_results %>% 
  select(Analysis, Effect, SE, ci95, `p (sig)`, Q, df, pQ, i2, tau2)

raw_xt <- xtable(
  raw_tab,
  caption = "Raw‐level meta‐analysis results for each outcome. Hedges’ $g$, standard errors, 95\\% CIs, p‐values, and heterogeneity statistics (Q, df, pQ, I\\textsuperscript{2}, $\\tau^2$).",
  label   = "tab:meta_raw",
  align   = c("l", "l", rep("r", 9))
)

print(
  raw_xt,
  type           = "latex",
  file           = file.path(output_dir, "all_meta_analyses_raw.tex"),
  include.rownames = FALSE,
  booktabs       = TRUE,
  sanitize.text.function = identity  # allows your LaTeX in ci95 (with $…$) to pass through
)

# — dataset-Labels hinzufügen
df_performance <- df_performance %>% mutate(dataset = "Creative Performance")
df_diversity   <- df_diversity   %>% mutate(dataset = "Creative Diversity")
df_versus      <- df_versus      %>% mutate(dataset = "Human_vs_AI")

# --------- 1 große Tabelle aller Moderatoren ---------
moderators <- c(
  "GenAI_Type",
  "GenAI_Model",
  "Participants",
  "Recruitment_Source",
  "Task_Type",
  "Creativity_Measurement",
  "Measurement_Evaluator"
)

df_all <- bind_rows(df_performance, df_diversity, df_versus)

combined_table <- bind_rows(
  lapply(moderators, function(mod) {
    df_all %>%
      filter(!is.na(.data[[mod]])) %>%
      group_by(
        Moderator      = mod,
        Characteristic = .data[[mod]]
      ) %>%
      summarise(
        Creative_Performance = sum(dataset == "Creative Performance"),
        Creative_Diversity   = sum(dataset == "Creative Diversity"),
        Human_vs_AI          = sum(dataset == "Human_vs_AI"),
        Total                = n(),
        .groups = "drop"
      )
  })
)

# — für jeden Moderator eine eigene LaTeX-Tabelle —
for(mod in moderators) {
  sub <- combined_table %>% 
    filter(Moderator == mod) %>% 
    select(-Moderator) %>% 
    arrange(desc(Total))
    
  
  xt <- xtable(
    sub,
    label   = paste0("tab:", mod),
    align   = c("l","l","r","r","r","r")
  )
  print(
    xt,
    type             = "latex",
    file             = file.path(output_dir, paste0("table_", mod, ".tex")),
    include.rownames = FALSE,
    booktabs         = TRUE,
    caption.placement= "top",
    sanitize.text.function = identity
  )
}