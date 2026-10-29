############################################################
# PUBLICATION-QUALITY KAPLAN-MEIER PANEL - EFS
#
# Genes:
# TACSTD2, EGFR, MET
#
# FIXED CUTOFFS:
# TACSTD2 = 7.73
# EGFR    = 5.74
# MET     = 4.65
#
# COLORS:
# High = BLUE
# Low  = RED
#
# P-VALUE:
# Exact log-rank p-value shown in the upper-right corner
#
# RISK TABLE:
# Clean custom "Number at risk" table
# First and last values are kept completely inside the box
#
# INPUT:
# ACC292_genes_clinics_v3.txt
#
# OUTPUT:
# 1. EFS_KM_panel_TACSTD2_EGFR_MET_publication.pdf
# 2. EFS_KM_panel_TACSTD2_EGFR_MET_publication.png
# 3. EFS_KM_TACSTD2_EGFR_MET_summary.csv
############################################################


# ==========================================================
# 1. INSTALL AND LOAD PACKAGES
# ==========================================================

packages <- c(
  "survival",
  "ggplot2",
  "cowplot",
  "ggtext"
)


new.packages <- packages[
  !(packages %in% rownames(installed.packages()))
]


if(length(new.packages) > 0){

  install.packages(
    new.packages,
    dependencies = TRUE
  )

}


library(survival)
library(ggplot2)
library(cowplot)
library(ggtext)


# ==========================================================
# 2. READ DATA
# ==========================================================

file.name <- "ACC266_genes_clinics.txt"


dat <- read.delim(

  file.name,

  header = TRUE,

  sep = "\t",

  check.names = FALSE,

  stringsAsFactors = FALSE,

  na.strings = c(
    "",
    "NA",
    "NaN"
  )

)


# ==========================================================
# 3. CHECK SURVIVAL VARIABLES
# ==========================================================

required.variables <- c(
  "EFS_time",
  "EFS_status"
)


if(
  !all(
    required.variables %in% colnames(dat)
  )
){

  stop(
    "EFS_time and/or EFS_status are missing from the dataset."
  )

}


# Convert survival variables to numeric

dat$EFS_time <- suppressWarnings(
  as.numeric(
    dat$EFS_time
  )
)


dat$EFS_status <- suppressWarnings(
  as.numeric(
    dat$EFS_status
  )
)


# ==========================================================
# 4. DEFINE GENES
# ==========================================================

genes <- c(
  "TACSTD2",
  "EGFR",
  "MET"
)


# ==========================================================
# 5. FIXED CUTOFFS
#
# These values are imposed and NOT recalculated.
# ==========================================================

cutoffs <- c(

  TACSTD2 = 7.73,

  EGFR = 5.74,

  MET = 4.65

)


# ==========================================================
# 6. CHECK GENES
# ==========================================================

missing.genes <- setdiff(
  genes,
  colnames(dat)
)


if(length(missing.genes) > 0){

  stop(

    paste(

      "Missing genes:",

      paste(
        missing.genes,
        collapse = ", "
      )

    )

  )

}


# Convert expression values to numeric

for(gene in genes){

  dat[[gene]] <- suppressWarnings(
    as.numeric(
      dat[[gene]]
    )
  )

}


# ==========================================================
# 7. COLORS
# ==========================================================

group.colors <- c(

  "High" = "blue",

  "Low" = "red"

)


# ==========================================================
# 8. TIME AXIS SETTINGS
# ==========================================================

BREAK.TIME <- 50


valid.times <- dat$EFS_time[

  is.finite(dat$EFS_time) &

    dat$EFS_time >= 0

]


MAX.TIME <- max(
  valid.times,
  na.rm = TRUE
)


# Round maximum follow-up upward

X.MAX <- ceiling(
  MAX.TIME / 10
) * 10


# Time points used in the Number-at-Risk table

RISK.TIMES <- seq(

  0,

  floor(
    X.MAX / BREAK.TIME
  ) * BREAK.TIME,

  by = BREAK.TIME

)


# ==========================================================
# 9. INTERNAL POSITIONING OF EDGE RISK NUMBERS
#
# Keeps the first and last values completely inside
# the Number-at-Risk rectangle.
# ==========================================================

INNER.NUDGE <- BREAK.TIME * 0.04


cat(
  "\nMaximum follow-up:",
  MAX.TIME,
  "\n"
)


cat(
  "X-axis maximum:",
  X.MAX,
  "\n"
)


cat(
  "Risk-table times:",
  RISK.TIMES,
  "\n"
)


# ==========================================================
# 10. FUNCTION TO FORMAT EXACT P-VALUE
#
# IMPORTANT:
# No thresholds such as p < 0.0001 are used.
#
# Small values:
# p = 1.62e-05
#
# Larger values:
# p = 0.0123
# ==========================================================

format.pvalue <- function(p){


  if(is.na(p)){

    return(
      "p = NA"
    )

  }


  if(p < 0.001){

    return(

      paste0(

        "p = ",

        formatC(
          p,
          format = "e",
          digits = 2
        )

      )

    )

  }


  return(

    paste0(

      "p = ",

      formatC(
        p,
        format = "f",
        digits = 4
      )

    )

  )

}


# ==========================================================
# 11. FUNCTION TO CREATE KM DATA
# ==========================================================

create.km.data <- function(fit){


  ss <- summary(

    fit,

    censored = TRUE

  )


  km.data <- data.frame(

    time =
      ss$time,

    survival =
      ss$surv,

    n.risk =
      ss$n.risk,

    n.event =
      ss$n.event,

    n.censor =
      ss$n.censor,

    strata =
      as.character(
        ss$strata
      )

  )


  # --------------------------------------------------------
  # Clean strata names
  # --------------------------------------------------------

  km.data$Group <- sub(

    "^Group=",

    "",

    km.data$strata

  )


  km.data$Group <- factor(

    km.data$Group,

    levels = c(
      "High",
      "Low"
    )

  )


  # --------------------------------------------------------
  # Add time zero
  # --------------------------------------------------------

  initial <- data.frame(

    time = c(
      0,
      0
    ),

    survival = c(
      1,
      1
    ),

    n.risk = c(
      NA,
      NA
    ),

    n.event = c(
      0,
      0
    ),

    n.censor = c(
      0,
      0
    ),

    strata = c(
      "Group=High",
      "Group=Low"
    ),

    Group = factor(

      c(
        "High",
        "Low"
      ),

      levels = c(
        "High",
        "Low"
      )

    )

  )


  km.data <- rbind(
    initial,
    km.data
  )


  km.data <- km.data[

    order(
      km.data$Group,
      km.data$time
    ),

    ,

    drop = FALSE

  ]


  rownames(
    km.data
  ) <- NULL


  return(
    km.data
  )

}


# ==========================================================
# 12. FUNCTION TO CREATE NUMBER-AT-RISK DATA
# ==========================================================

create.risk.data <- function(
    fit,
    times
){


  rs <- summary(

    fit,

    times = times,

    extend = TRUE

  )


  risk.data <- data.frame(

    time =
      rs$time,

    n.risk =
      rs$n.risk,

    strata =
      as.character(
        rs$strata
      )

  )


  # --------------------------------------------------------
  # Clean group names
  # --------------------------------------------------------

  risk.data$Group <- sub(

    "^Group=",

    "",

    risk.data$strata

  )


  risk.data$Group <- factor(

    risk.data$Group,

    levels = c(
      "High",
      "Low"
    )

  )


  return(
    risk.data
  )

}


# ==========================================================
# 13. FUNCTION TO CREATE COMPLETE PANEL FOR ONE GENE
# ==========================================================

create.gene.panel <- function(

    data,

    gene,

    cutoff,

    xmax,

    risk.times,

    inner.nudge

){


  # ========================================================
  # 13.1 PREPARE DATA
  # ========================================================

  d <- data.frame(

    expression = suppressWarnings(
      as.numeric(
        data[[gene]]
      )
    ),

    time = suppressWarnings(
      as.numeric(
        data$EFS_time
      )
    ),

    status = suppressWarnings(
      as.numeric(
        data$EFS_status
      )
    )

  )


  # --------------------------------------------------------
  # Keep evaluable patients
  # --------------------------------------------------------

  d <- d[

    complete.cases(d) &

      d$status %in% c(
        0,
        1
      ) &

      d$time > 0,

    ,

    drop = FALSE

  ]


  # ========================================================
  # 13.2 DEFINE GROUPS
  #
  # High > cutoff
  # Low <= cutoff
  # ========================================================

  d$Group <- ifelse(

    d$expression > cutoff,

    "High",

    "Low"

  )


  d$Group <- factor(

    d$Group,

    levels = c(
      "High",
      "Low"
    )

  )


  # ========================================================
  # 13.3 SAMPLE COUNTS
  # ========================================================

  n.high <- sum(
    d$Group == "High"
  )


  n.low <- sum(
    d$Group == "Low"
  )


  events.high <- sum(

    d$status[
      d$Group == "High"
    ] == 1

  )


  events.low <- sum(

    d$status[
      d$Group == "Low"
    ] == 1

  )


  # ========================================================
  # 14. KAPLAN-MEIER FIT
  # ========================================================

  fit <- survfit(

    Surv(
      time,
      status
    ) ~ Group,

    data = d

  )


  # ========================================================
  # 15. LOG-RANK TEST
  # ========================================================

  lr <- survdiff(

    Surv(
      time,
      status
    ) ~ Group,

    data = d

  )


  p.value <- 1 - pchisq(

    lr$chisq,

    df =
      length(
        lr$n
      ) - 1

  )


  # Exact p-value text

  p.label <- format.pvalue(
    p.value
  )


  # ========================================================
  # 16. COX MODEL
  #
  # High is the reference group.
  # HR therefore corresponds to Low versus High.
  # ========================================================

  cox.fit <- coxph(

    Surv(
      time,
      status
    ) ~ Group,

    data = d

  )


  cox.summary <- summary(
    cox.fit
  )


  HR.low.vs.high <- cox.summary$conf.int[
    1,
    "exp(coef)"
  ]


  CI.low <- cox.summary$conf.int[
    1,
    "lower .95"
  ]


  CI.high <- cox.summary$conf.int[
    1,
    "upper .95"
  ]


  # ========================================================
  # 17. KM DATA
  # ========================================================

  km.data <- create.km.data(
    fit
  )


  # ========================================================
  # 18. CENSOR DATA
  # ========================================================

  censor.data <- km.data[

    !is.na(
      km.data$n.censor
    ) &

      km.data$n.censor > 0,

    ,

    drop = FALSE

  ]


  # ========================================================
  # 19. MAIN KM PLOT
  # ========================================================

  km.plot <- ggplot(

    km.data,

    aes(

      x = time,

      y = survival,

      color = Group,

      group = Group

    )

  ) +


    # ------------------------------------------------------
    # KM curves
    # ------------------------------------------------------

    geom_step(

      linewidth = 1.10,

      direction = "hv",

      lineend = "butt"

    ) +


    # ------------------------------------------------------
    # Censor marks
    # ------------------------------------------------------

    geom_point(

      data = censor.data,

      aes(

        x = time,

        y = survival,

        color = Group

      ),

      inherit.aes = FALSE,

      shape = 3,

      size = 1.65,

      stroke = 0.7

    ) +


    # ------------------------------------------------------
    # Colors
    # ------------------------------------------------------

    scale_color_manual(

      values =
        group.colors,

      breaks = c(
        "High",
        "Low"
      ),

      labels = c(
        "High",
        "Low"
      ),

      name = NULL

    ) +


    # ------------------------------------------------------
    # X axis
    # ------------------------------------------------------

    scale_x_continuous(

      limits = c(
        0,
        xmax
      ),

      breaks =
        risk.times,

      expand = expansion(
        mult = c(
          0,
          0
        )
      )

    ) +


    # ------------------------------------------------------
    # Y axis
    # ------------------------------------------------------

    scale_y_continuous(

      limits = c(
        0,
        1
      ),

      breaks = c(
        0,
        0.25,
        0.50,
        0.75,
        1.00
      ),

      labels = c(
        "0.00",
        "0.25",
        "0.50",
        "0.75",
        "1.00"
      ),

      expand = expansion(
        mult = c(
          0,
          0.02
        )
      )

    ) +


    # ======================================================
    # EXACT P-VALUE
    #
    # UPPER-RIGHT CORNER
    #
    # x = 96% of maximum X
    # y = 95% survival
    #
    # hjust = 1 means text extends to the LEFT,
    # remaining completely inside the graph.
    # ======================================================

    annotate(

      "text",

      x =
        xmax * 0.96,

      y =
        0.95,

      label =
        p.label,

      hjust =
        1,

      vjust =
        1,

      size =
        4.0,

      fontface =
        "plain",

      color =
        "black"

    ) +


    # ------------------------------------------------------
    # Labels
    # ------------------------------------------------------

    labs(

      title =
        gene,

      subtitle = paste0(

        "Cutoff = ",

        sprintf(
          "%.2f",
          cutoff
        )

      ),

      x =
        NULL,

      y =
        "Event-Free Survival"

    ) +


    # ------------------------------------------------------
    # Theme
    # ------------------------------------------------------

    theme_classic(
      base_size = 11
    ) +


    theme(


      plot.title = element_text(

        hjust = 0.5,

        face = "bold",

        size = 15,

        margin = margin(
          b = 2
        )

      ),


      plot.subtitle = element_text(

        hjust = 0.5,

        size = 9,

        margin = margin(
          b = 8
        )

      ),


      axis.title.y = element_text(

        size = 11,

        face = "bold",

        margin = margin(
          r = 8
        )

      ),


      axis.text.x = element_text(
        size = 9
      ),


      axis.text.y = element_text(
        size = 9
      ),


      axis.line = element_line(
        linewidth = 0.45
      ),


      axis.ticks = element_line(
        linewidth = 0.4
      ),


      # Subtle horizontal grid

      panel.grid.major.y = element_line(

        color = "grey90",

        linewidth = 0.35

      ),


      panel.grid.major.x =
        element_blank(),


      panel.grid.minor =
        element_blank(),


      legend.position =
        "top",


      legend.text = element_text(
        size = 10
      ),


      legend.key.width =
        grid::unit(
          1.3,
          "cm"
        ),


      plot.margin = margin(

        t = 4,

        r = 8,

        b = 3,

        l = 8

      )

    ) +


    guides(

      color = guide_legend(

        override.aes = list(

          linewidth = 1.4,

          shape = NA

        )

      )

    )


  # ========================================================
  # 20. CREATE NUMBER-AT-RISK DATA
  # ========================================================

  risk.data <- create.risk.data(

    fit =
      fit,

    times =
      risk.times

  )


  # ========================================================
  # 21. DIVIDE RISK TABLE INTO THREE PARTS
  #
  # FIRST COLUMN:
  # left justified inside the box
  #
  # MIDDLE:
  # centered
  #
  # LAST:
  # right justified inside the box
  # ========================================================

  first.time <- min(
    risk.times
  )


  last.time <- max(
    risk.times
  )


  risk.first <- risk.data[

    risk.data$time == first.time,

    ,

    drop = FALSE

  ]


  risk.middle <- risk.data[

    risk.data$time > first.time &

      risk.data$time < last.time,

    ,

    drop = FALSE

  ]


  risk.last <- risk.data[

    risk.data$time == last.time,

    ,

    drop = FALSE

  ]


  # ========================================================
  # 22. NUMBER-AT-RISK TABLE
  # ========================================================

  risk.plot <- ggplot() +


    # ------------------------------------------------------
    # FIRST COLUMN
    #
    # All digits stay inside the rectangle
    # ------------------------------------------------------

    geom_text(

      data =
        risk.first,

      aes(

        x =
          first.time + inner.nudge,

        y =
          Group,

        label =
          n.risk,

        color =
          Group

      ),

      hjust =
        0,

      size =
        4.25,

      fontface =
        "bold",

      show.legend =
        FALSE

    ) +


    # ------------------------------------------------------
    # MIDDLE COLUMNS
    # ------------------------------------------------------

    geom_text(

      data =
        risk.middle,

      aes(

        x =
          time,

        y =
          Group,

        label =
          n.risk,

        color =
          Group

      ),

      hjust =
        0.5,

      size =
        4.25,

      fontface =
        "bold",

      show.legend =
        FALSE

    ) +


    # ------------------------------------------------------
    # LAST COLUMN
    # ------------------------------------------------------

    geom_text(

      data =
        risk.last,

      aes(

        x =
          last.time - inner.nudge,

        y =
          Group,

        label =
          n.risk,

        color =
          Group

      ),

      hjust =
        1,

      size =
        4.25,

      fontface =
        "bold",

      show.legend =
        FALSE

    ) +


    # ------------------------------------------------------
    # Colors
    # ------------------------------------------------------

    scale_color_manual(

      values =
        group.colors,

      guide =
        "none"

    ) +


    # ------------------------------------------------------
    # X axis
    #
    # True time locations remain:
    # 0, 50, 100, 150...
    # ------------------------------------------------------

    scale_x_continuous(

      limits = c(
        0,
        xmax
      ),

      breaks =
        risk.times,

      labels =
        risk.times,

      expand = expansion(
        mult = c(
          0,
          0
        )
      )

    ) +


    # ------------------------------------------------------
    # Y axis
    #
    # Low on top
    # High below
    # ------------------------------------------------------

    scale_y_discrete(

      limits = c(
        "High",
        "Low"
      ),

      labels = c(

        "High" =
          "<span style='color:blue'><b>High</b></span>",

        "Low" =
          "<span style='color:red'><b>Low</b></span>"

      ),

      drop =
        FALSE

    ) +


    # ------------------------------------------------------
    # Labels
    # ------------------------------------------------------

    labs(

      title =
        "Number at risk",

      x =
        "EFS time",

      y =
        NULL

    ) +


    # ------------------------------------------------------
    # Risk-table theme
    # ------------------------------------------------------

    theme_bw(
      base_size = 10
    ) +


    theme(


      plot.title = element_text(

        hjust = 0,

        face = "bold",

        size = 10.5,

        margin = margin(
          b = 5
        )

      ),


      # ----------------------------------------------------
      # Colored High / Low labels
      # ----------------------------------------------------

      axis.text.y =
        ggtext::element_markdown(

          size = 10,

          hjust = 1,

          margin = margin(
            r = 8
          )

        ),


      # ----------------------------------------------------
      # Time labels
      # ----------------------------------------------------

      axis.text.x = element_text(

        size = 8.5,

        color = "black"

      ),


      axis.title.x = element_text(

        size = 9,

        face = "bold",

        margin = margin(
          t = 5
        )

      ),


      axis.title.y =
        element_blank(),


      axis.ticks.y =
        element_blank(),


      # ----------------------------------------------------
      # Vertical grid
      # ----------------------------------------------------

      panel.grid.major.x = element_line(

        color = "grey86",

        linewidth = 0.40

      ),


      # ----------------------------------------------------
      # Horizontal separator
      # ----------------------------------------------------

      panel.grid.major.y = element_line(

        color = "grey91",

        linewidth = 0.40

      ),


      panel.grid.minor =
        element_blank(),


      # ----------------------------------------------------
      # Table rectangle
      # ----------------------------------------------------

      panel.border = element_rect(

        color = "grey30",

        fill = NA,

        linewidth = 0.55

      ),


      axis.ticks.x = element_line(
        linewidth = 0.35
      ),


      plot.margin = margin(

        t = 2,

        r = 10,

        b = 4,

        l = 10

      )

    )


  # ========================================================
  # 23. REMOVE INDIVIDUAL LEGEND
  # ========================================================

  km.no.legend <- km.plot +

    theme(
      legend.position = "none"
    )


  # ========================================================
  # 24. COMBINE KM + NUMBER-AT-RISK TABLE
  # ========================================================

  combined <- plot_grid(

    km.no.legend,

    risk.plot,

    ncol = 1,

    align = "v",

    axis = "lr",

    rel_heights = c(

      4.0,

      1.30

    )

  )


  # ========================================================
  # 25. RETURN
  # ========================================================

  return(

    list(

      combined =
        combined,

      km.plot =
        km.plot,

      risk.plot =
        risk.plot,

      gene =
        gene,

      cutoff =
        cutoff,

      n_high =
        n.high,

      n_low =
        n.low,

      events_high =
        events.high,

      events_low =
        events.low,

      logrank_p =
        p.value,

      HR_low_vs_high =
        HR.low.vs.high,

      CI_low =
        CI.low,

      CI_high =
        CI.high

    )

  )

}


# ==========================================================
# 26. CREATE ALL THREE GENE PANELS
# ==========================================================

results <- list()


for(gene in genes){


  cat(

    "\nCreating ",

    gene,

    " | cutoff = ",

    cutoffs[gene],

    "\n",

    sep = ""

  )


  results[[gene]] <- create.gene.panel(

    data =
      dat,

    gene =
      gene,

    cutoff =
      cutoffs[gene],

    xmax =
      X.MAX,

    risk.times =
      RISK.TIMES,

    inner.nudge =
      INNER.NUDGE

  )

}


# ==========================================================
# 27. CREATE SINGLE COMMON LEGEND
# ==========================================================

legend.source <- results[["TACSTD2"]]$km.plot +

  theme(

    legend.position =
      "top",

    legend.justification =
      "center",

    legend.box.just =
      "center"

  )


common.legend <- cowplot::get_legend(
  legend.source
)


# ==========================================================
# 28. EXTRACT PANELS
# ==========================================================

panel.TACSTD2 <-
  results[["TACSTD2"]]$combined


panel.EGFR <-
  results[["EGFR"]]$combined


panel.MET <-
  results[["MET"]]$combined


# ==========================================================
# 29. COMBINE THREE GENE PANELS
# ==========================================================

three.genes <- plot_grid(

  panel.TACSTD2,

  panel.EGFR,

  panel.MET,

  ncol = 3,

  align = "hv",

  axis = "tblr",

  rel_widths = c(
    1,
    1,
    1
  )

)


# ==========================================================
# 30. GLOBAL TITLE
# ==========================================================

global.title <- ggdraw() +

  draw_label(

    "Event-Free Survival according to biomarker expression",

    x =
      0.5,

    hjust =
      0.5,

    fontface =
      "bold",

    size =
      17

  )


# ==========================================================
# 31. FINAL PANEL
# ==========================================================

final.panel <- plot_grid(

  global.title,

  common.legend,

  three.genes,

  ncol = 1,

  rel_heights = c(

    0.07,

    0.07,

    0.86

  )

)


# ==========================================================
# 32. SAVE PDF
# ==========================================================

pdf.file <-
  "EFS_KM_panel_TACSTD2_EGFR_MET_publication.pdf"


if(
  capabilities("cairo")
){

  ggsave(

    filename =
      pdf.file,

    plot =
      final.panel,

    device =
      cairo_pdf,

    width =
      16.5,

    height =
      7.5,

    units =
      "in",

    limitsize =
      FALSE

  )

} else {


  grDevices::pdf(

    file =
      pdf.file,

    width =
      16.5,

    height =
      7.5,

    useDingbats =
      FALSE

  )


  print(
    final.panel
  )


  dev.off()

}


# ==========================================================
# 33. SAVE HIGH-RESOLUTION PNG
# ==========================================================

png.file <-
  "EFS_KM_panel_TACSTD2_EGFR_MET_publication.png"


ggsave(

  filename =
    png.file,

  plot =
    final.panel,

  width =
    16.5,

  height =
    7.5,

  units =
    "in",

  dpi =
    400,

  limitsize =
    FALSE

)


# ==========================================================
# 34. CREATE SUMMARY TABLE
# ==========================================================

summary.results <- do.call(

  rbind,

  lapply(

    results,

    function(x){

      data.frame(

        Gene =
          x$gene,

        Fixed_cutoff =
          x$cutoff,

        High_n =
          x$n_high,

        High_events =
          x$events_high,

        Low_n =
          x$n_low,

        Low_events =
          x$events_low,

        LogRank_p =
          x$logrank_p,

        HR_Low_vs_High =
          x$HR_low_vs_high,

        CI95_low =
          x$CI_low,

        CI95_high =
          x$CI_high

      )

    }

  )

)


rownames(
  summary.results
) <- NULL


# ==========================================================
# 35. SAVE SUMMARY
# ==========================================================

write.csv(

  summary.results,

  file =
    "EFS_KM_TACSTD2_EGFR_MET_summary.csv",

  row.names =
    FALSE

)


# ==========================================================
# 36. PRINT RESULTS
# ==========================================================

cat(
  "\n\n========================================\n"
)


cat(
  "KAPLAN-MEIER RESULTS\n"
)


cat(
  "========================================\n"
)


print(
  summary.results
)


cat(
  "\nPDF saved as:\n",
  pdf.file,
  "\n"
)


cat(
  "\nPNG saved as:\n",
  png.file,
  "\n"
)


# ==========================================================
# 37. DISPLAY FINAL PANEL
# ==========================================================

print(
  final.panel
)


############################################################
# END OF SCRIPT
############################################################