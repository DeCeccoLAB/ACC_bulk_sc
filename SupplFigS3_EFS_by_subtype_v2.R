#!/usr/bin/env Rscript
# =============================================================================
#  Supplementary Figure S3 (v2)
#  Event-free survival by ACC molecular subtype and expression of
#  therapeutically relevant genes
#
#  Usage:
#    Rscript SupplFigS3_EFS_by_subtype_v2.R ACC292_genes_clinics_v3.txt [out_dir]
#
#  Input columns: GeneSymbol (sample ID), 12 genes, ClustersB (ACC1/ACC2),
#                 EFS_status (1 = event, 0 = censored), EFS_time (months)
#  Output:        Supplementary_FigureS3_EFS_by_subtype.pdf (one page)
#                 SupplFigS3_cox_continuous.tsv, SupplFigS3_cox_cutoff.tsv
#  Packages:      survival, ggplot2 (>= 3.4), patchwork, dplyr
# =============================================================================

## ---- 0. Setup ---------------------------------------------------------------
for (loc in c("en_US.UTF-8", "C.UTF-8", "English_United States.utf8"))
  if (nzchar(suppressWarnings(Sys.setlocale("LC_CTYPE", loc)))) break

suppressPackageStartupMessages({
  library(survival); library(ggplot2); library(patchwork); library(dplyr)
})

args    <- commandArgs(trailingOnly = TRUE)
infile  <- if (length(args) >= 1) args[1] else "ACC292_genes_clinics_v3.txt"
out_dir <- if (length(args) >= 2) args[2] else "."
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

GENES   <- c("VTCN1", "MYB", "NECTIN4", "TACSTD2", "STING1", "KIT",
             "FAP", "RBPJ", "CDK4", "CDK6", "EGFR", "MET")
CUTOFFS <- c(TACSTD2 = 7.73, EGFR = 5.74, MET = 4.65)     # Figure 3B cutoffs

## Page geometry (mm). Panel B gets a fixed vertical pitch per estimate so that
## symbols can never touch: PITCH_MM must exceed the largest symbol (~2 mm).
PAGE_W   <- 183
H_A      <- 66      # row A
H_B      <- 205     # row B (forest + table + legend)
H_C      <- 80      # row C
PAGE_H   <- H_A + H_B + H_C

## Colours / text
COL <- c(I = "#E0612B", II = "#2A6FC9", crude = "#8A9099", adj = "#1F2328")
INK <- "#1F2328"; MUTED <- "#60666E"; BAND <- "#F0F1F3"
FONT <- "sans"      # "Arial"/"Helvetica" if available
FS   <- 7           # base font size (pt)
txt_mm <- function(pt) pt / .pt

theme_s3 <- function() {
  theme_classic(base_size = FS, base_family = FONT) +
    theme(axis.line  = element_line(linewidth = 0.3, colour = MUTED),
          axis.ticks = element_line(linewidth = 0.3, colour = MUTED),
          axis.text  = element_text(size = FS - 0.5, colour = INK),
          axis.title = element_text(size = FS, colour = INK),
          plot.title = element_text(size = FS + 2, face = "bold", colour = INK,
                                    margin = margin(0, 0, 2, 0)),
          plot.subtitle = element_text(size = FS - 0.8, colour = INK, lineheight = 1.15),
          legend.text = element_text(size = FS - 0.5, colour = INK),
          legend.background = element_blank(),
          legend.key = element_blank(),
          plot.margin = margin(3, 5, 3, 3))
}

## p-value formatting with real superscripts (needs UTF-8)
SUP <- c("-" = "⁻", "0" = "⁰", "1" = "¹", "2" = "²", "3" = "³",
         "4" = "⁴", "5" = "⁵", "6" = "⁶", "7" = "⁷", "8" = "⁸", "9" = "⁹")
fp <- function(p) vapply(p, function(x) {
  if (is.na(x)) return("NA")
  if (x >= 0.1)   return(sprintf("%.2f", x))
  if (x >= 0.001) return(sprintf("%.3f", x))
  e <- floor(log10(x)); m <- x / 10^e
  paste0(sprintf("%.1f", m), " × 10", paste(SUP[strsplit(as.character(e), "")[[1]]], collapse = ""))
}, character(1))
fci <- function(h, l, u) sprintf("%.2f (%.2f–%.2f)", h, l, u)

## ---- 1. Data ------------------------------------------------------------------
raw <- read.delim(infile, check.names = FALSE, stringsAsFactors = FALSE)
raw$Subtype <- factor(ifelse(raw$ClustersB == "ACC1", "ACC-I", "ACC-II"), c("ACC-I", "ACC-II"))
delta_med <- sapply(GENES, function(g)
  median(raw[[g]][raw$Subtype == "ACC-II"]) - median(raw[[g]][raw$Subtype == "ACC-I"]))

has_efs <- !is.na(raw$EFS_time) & !is.na(raw$EFS_status)
d <- raw[has_efs & raw$EFS_time > 0, ]       # time = 0 excluded (as in CutoffFinder)
d$cohort <- dplyr::case_when(grepl("^ACC_", d$GeneSymbol) ~ "Ferrarotto",
                             grepl("^GSM4632", d$GeneSymbol) ~ "Persson",
                             grepl("^P[0-9]", d$GeneSymbol) ~ "Dou", TRUE ~ "other")
for (g in GENES) d[[paste0(g, "_z")]] <- as.numeric(scale(d[[g]]))
N  <- table(d$Subtype); EV <- tapply(d$EFS_status, d$Subtype, sum)
message(sprintf("EFS available %d; time = 0 excluded %d; analysed %d (%d events)",
                sum(has_efs), sum(has_efs & raw$EFS_time <= 0), nrow(d), sum(d$EFS_status)))

## ---- 2. Statistics ------------------------------------------------------------
S <- function(rhs, data = d) coxph(as.formula(paste("Surv(EFS_time, EFS_status) ~", rhs)), data = data)
hr <- function(fit, term) {
  b <- coef(fit)[term]; ci <- exp(confint(fit)[term, ])
  c(HR = unname(exp(b)), lo = unname(ci[1]), hi = unname(ci[2]),
    p = unname(summary(fit)$coefficients[term, "Pr(>|z|)"]))
}
lrt <- function(f0, f1) anova(f0, f1)[2, "Pr(>|Chi|)"]

# subtype alone
f_sub  <- S("Subtype")
hr_sub <- hr(f_sub, "SubtypeACC-II")
p_lr   <- pchisq(survdiff(Surv(EFS_time, EFS_status) ~ Subtype, d)$chisq, 1, lower.tail = FALSE)
km_sub <- survfit(Surv(EFS_time, EFS_status) ~ Subtype, d)
med    <- summary(km_sub)$table[, "median"]

# continuous models, HR per 1 SD
cont <- bind_rows(lapply(GENES, function(g) {
  z <- paste0(g, "_z")
  fa <- S(paste(z, "+ Subtype"))
  m <- list(crude = hr(S(z), z), adj = hr(fa, z),
            adj_cohort = hr(S(paste(z, "+ Subtype + strata(cohort)")), z),
            ACC1 = hr(S(z, d[d$Subtype == "ACC-I", ]), z),
            ACC2 = hr(S(z, d[d$Subtype == "ACC-II", ]), z))
  out <- data.frame(gene = g, delta_median = unname(delta_med[g]))
  for (k in names(m)) out[paste0(k, "_", names(m[[k]]))] <- as.list(m[[k]])
  out$p_interaction <- lrt(fa, S(paste(z, "* Subtype")))
  out
}))
for (k in c("crude", "adj", "ACC1", "ACC2")) cont[[paste0(k, "_FDR")]] <- p.adjust(cont[[paste0(k, "_p")]], "BH")
write.table(cont, file.path(out_dir, "SupplFigS3_cox_continuous.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# dichotomised at Figure 3B cutoffs
cutt <- bind_rows(lapply(names(CUTOFFS), function(g) {
  d$hi <- factor(ifelse(d[[g]] > CUTOFFS[[g]], "High", "Low"), c("Low", "High"))
  tb <- table(d$Subtype, d$hi)
  fa <- coxph(Surv(EFS_time, EFS_status) ~ hi + Subtype, d)
  a  <- hr(fa, "hiHigh")
  data.frame(gene = g, cutoff = CUTOFFS[[g]],
             hi_I = tb["ACC-I", "High"],  lo_I = tb["ACC-I", "Low"],
             hi_II = tb["ACC-II", "High"], lo_II = tb["ACC-II", "Low"],
             fisher_p = fisher.test(tb)$p.value,
             adj_HR = a["HR"], adj_lo = a["lo"], adj_hi = a["hi"], adj_p = a["p"],
             strat_logrank_p = pchisq(survdiff(Surv(EFS_time, EFS_status) ~ hi + strata(Subtype), d)$chisq,
                                      1, lower.tail = FALSE),
             p_interaction = lrt(fa, coxph(Surv(EFS_time, EFS_status) ~ hi * Subtype, d)))
}))
write.table(cutt, file.path(out_dir, "SupplFigS3_cox_cutoff.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

## ---- 3. Helpers for Kaplan-Meier ---------------------------------------------
km_steps <- function(fit) {
  st <- if (is.null(fit$strata)) rep("all", length(fit$time)) else
    rep(sub("^[^=]*=", "", names(fit$strata)), fit$strata)
  df <- data.frame(time = fit$time, surv = fit$surv, cens = fit$n.censor > 0, strata = st)
  bind_rows(data.frame(time = 0, surv = 1, cens = FALSE, strata = unique(st)), df)
}
XMAX <- 180

## ---- 4. Panel A ---------------------------------------------------------------
kmA <- km_steps(km_sub)
kmA$strata <- factor(kmA$strata, c("ACC-I", "ACC-II"))
labA <- sprintf("%s (n = %d; %d events)", levels(d$Subtype), N, EV)

pA1 <- ggplot(kmA, aes(time, surv, colour = strata)) +
  geom_step(linewidth = 0.6) +
  geom_point(data = subset(kmA, cens), shape = 3, size = 1.1, stroke = 0.4, show.legend = FALSE) +
  scale_colour_manual(values = unname(COL[c("I", "II")]), labels = labA, name = NULL) +
  scale_x_continuous(breaks = seq(0, XMAX, 30)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25)) +
  coord_cartesian(xlim = c(-14, XMAX + 3), ylim = c(0, 1.02), expand = FALSE) +
  annotate("text", x = XMAX, y = 1.0, hjust = 1, vjust = 1, size = txt_mm(FS - 0.5),
           family = FONT, colour = INK, lineheight = 1.15,
           label = sprintf("ACC-II vs ACC-I: HR %s\nlog-rank p = %s\nmedian EFS %.1f vs %.1f months",
                           fci(hr_sub[1], hr_sub[2], hr_sub[3]), fp(p_lr), med[2], med[1])) +
  labs(title = "A", x = NULL, y = "EFS probability") +
  theme_s3() +
  theme(legend.position = c(1, 0.48), legend.justification = c(1, 0.5),
        legend.key.width = unit(5, "mm"), legend.key.height = unit(3, "mm"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.margin = margin(3, 5, 0, 3))

risk <- expand.grid(time = seq(0, XMAX, 30), Subtype = levels(d$Subtype))
risk$n <- mapply(function(t, s) sum(d$EFS_time >= t & d$Subtype == s), risk$time, risk$Subtype)
risk$Subtype <- factor(risk$Subtype, rev(levels(d$Subtype)))
pA1r <- ggplot(risk, aes(time, Subtype, label = n, colour = Subtype)) +
  geom_text(size = txt_mm(FS - 0.5), family = FONT, show.legend = FALSE) +
  scale_colour_manual(values = c(`ACC-I` = unname(COL["I"]), `ACC-II` = unname(COL["II"]))) +
  scale_x_continuous(breaks = seq(0, XMAX, 30)) +
  coord_cartesian(xlim = c(-14, XMAX + 3), expand = TRUE) +
  labs(x = "Months from diagnosis", y = NULL, subtitle = "Number at risk") +
  theme_s3() +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_text(colour = unname(COL[c("II", "I")]), face = "bold"),
        plot.subtitle = element_text(colour = MUTED), plot.margin = margin(0, 5, 3, 3))

shr <- bind_rows(
  data.frame(gene = cutt$gene, Subtype = "ACC-I",  pct = 100 * cutt$hi_I  / (cutt$hi_I + cutt$lo_I)),
  data.frame(gene = cutt$gene, Subtype = "ACC-II", pct = 100 * cutt$hi_II / (cutt$hi_II + cutt$lo_II)))
shr$gene <- factor(shr$gene, rev(names(CUTOFFS)))
shr$Subtype <- factor(shr$Subtype, c("ACC-II", "ACC-I"))
pA2 <- ggplot(shr, aes(pct, gene, fill = Subtype)) +
  geom_col(position = position_dodge(0.8), width = 0.72) +
  geom_text(aes(label = sprintf("%.0f%%", pct)), position = position_dodge(0.8),
            hjust = -0.2, size = txt_mm(FS - 1), family = FONT, colour = INK) +
  scale_fill_manual(values = c(`ACC-I` = unname(COL["I"]), `ACC-II` = unname(COL["II"])),
                    breaks = c("ACC-I", "ACC-II"), name = NULL) +
  scale_x_continuous(limits = c(0, 115), breaks = c(0, 50, 100), expand = c(0, 0)) +
  labs(title = " ", x = "% high-expressing tumors", y = NULL,
       subtitle = sprintf("Figure 3B cutoffs\nFisher p ≤ %s for each gene", fp(max(cutt$fisher_p)))) +
  theme_s3() +
  theme(legend.position = "bottom", legend.key.size = unit(3, "mm"),
        axis.text.y = element_text(face = "bold"), plot.subtitle = element_text(colour = MUTED))

rowA <- ((pA1 / pA1r) + plot_layout(heights = c(1, 0.24))) | pA2
rowA <- rowA + plot_layout(widths = c(2.15, 1))

## ---- 5. Panel B (forest) ------------------------------------------------------
## One row per estimate. Each gene block = 4 rows (pitch 1) + a gap of GAP rows.
EST   <- c("crude", "adj", "ACC1", "ACC2")
ELAB  <- c(crude = sprintf("Crude, all patients (n = %d)", nrow(d)),
           adj   = "Adjusted for ACC subtype",
           ACC1  = sprintf("Within ACC-I (n = %d)", N[1]),
           ACC2  = sprintf("Within ACC-II (n = %d)", N[2]))
ECOL  <- c(crude = unname(COL["crude"]), adj = unname(COL["adj"]), ACC1 = unname(COL["I"]), ACC2 = unname(COL["II"]))
ESHP  <- c(crude = 21, adj = 23, ACC1 = 22, ACC2 = 24)    # circle, diamond, square, triangle
GAP   <- 1.6
ord   <- cont$gene[order(cont$delta_median)]            # ACC-I-enriched (top) -> ACC-II-enriched
nG    <- length(ord); BLOCK <- length(EST) + GAP
top_of <- function(i) (nG - i) * BLOCK + length(EST) - 1  # y of first row of block i

fr <- bind_rows(lapply(EST, function(k) data.frame(
  gene = cont$gene, est = k,
  HR = cont[[paste0(k, "_HR")]], lo = cont[[paste0(k, "_lo")]],
  hi = cont[[paste0(k, "_hi")]], p = cont[[paste0(k, "_p")]])))
fr$y   <- top_of(match(fr$gene, ord)) - (match(fr$est, EST) - 1)
fr$est <- factor(fr$est, EST)
fr$sig <- fr$p < 0.05
XL <- c(0.4, 2.2)
fr$loc <- pmax(fr$lo, XL[1]); fr$hic <- pmin(fr$hi, XL[2])

blk <- data.frame(gene = ord, ymid = top_of(seq_len(nG)) - (length(EST) - 1) / 2)
blk$ymin <- blk$ymid - BLOCK / 2; blk$ymax <- blk$ymid + BLOCK / 2
YL <- c(min(blk$ymin), max(blk$ymax))
bands <- blk[seq(1, nG, 2), ]

pB <- ggplot(fr) +
  geom_rect(data = bands, aes(xmin = XL[1], xmax = XL[2], ymin = ymin, ymax = ymax),
            fill = BAND, inherit.aes = FALSE) +
  geom_vline(xintercept = 1, linetype = "22", linewidth = 0.3, colour = MUTED) +
  geom_segment(aes(x = loc, xend = hic, y = y, yend = y, colour = est), linewidth = 0.45) +
  geom_point(aes(HR, y, colour = est, shape = est, fill = ifelse(sig, as.character(est), "open")),
             size = 1.7, stroke = 0.55) +
  scale_colour_manual(values = ECOL, labels = ELAB, name = NULL) +
  scale_shape_manual(values = ESHP, labels = ELAB, name = NULL) +
  scale_fill_manual(values = c(ECOL, open = "white"), guide = "none") +
  scale_x_log10(breaks = c(0.5, 0.75, 1, 1.5, 2), labels = c("0.5", "0.75", "1", "1.5", "2")) +
  scale_y_continuous(breaks = blk$ymid, labels = blk$gene) +
  coord_cartesian(xlim = XL, ylim = YL, expand = FALSE) +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE,
                               override.aes = list(fill = ECOL, size = 2, linewidth = 0.45)),
         shape  = guide_legend(nrow = 2, byrow = TRUE)) +
  labs(title = "B", y = NULL,
       x = "HR for an EFS event per 1-SD increase in expression (log scale)",
       subtitle = "← favourable          unfavourable →",
       caption = "Filled symbols: p < 0.05; open symbols: p ≥ 0.05") +
  theme_s3() +
  theme(axis.text.y = element_text(face = "bold", size = FS),
        axis.ticks.y = element_blank(),
        plot.subtitle = element_text(hjust = 0.5, colour = MUTED),
        plot.caption = element_text(hjust = 0, colour = MUTED, size = FS - 1),
        legend.position = "bottom", legend.justification = "left",
        legend.key.width = unit(6, "mm"), legend.key.height = unit(3.2, "mm"),
        legend.margin = margin(0, 0, 0, -45),
        plot.margin = margin(3, 0, 3, 3))

## Table aligned on the same y axis: one line per gene (block centre)
tcols <- data.frame(
  key = c("dm", "hr", "pc", "pa", "pi"),
  head = c("Δ median\n(II − I)", "Adjusted HR\n(95% CI)", "p\ncrude", "p\nadjusted", "p\ninteraction"),
  x = c(0.5, 1.95, 3.35, 4.45, 5.55))
tb <- data.frame(gene = cont$gene, y = blk$ymid[match(cont$gene, blk$gene)],
                 dm = sprintf("%+.2f", cont$delta_median),
                 hr = fci(cont$adj_HR, cont$adj_lo, cont$adj_hi),
                 pc = fp(cont$crude_p), pa = fp(cont$adj_p), pi = fp(cont$p_interaction),
                 bc = cont$crude_p < 0.05, ba = cont$adj_p < 0.05, bi = cont$p_interaction < 0.05)
tbl <- bind_rows(lapply(seq_len(nrow(tcols)), function(j) data.frame(
  x = tcols$x[j], y = tb$y, lab = tb[[tcols$key[j]]],
  bold = switch(tcols$key[j], pc = tb$bc, pa = tb$ba, pi = tb$bi, rep(FALSE, nrow(tb))))))
pBt <- ggplot(tbl) +
  geom_rect(data = bands, aes(xmin = 0, xmax = 6.15, ymin = ymin, ymax = ymax),
            fill = BAND, inherit.aes = FALSE) +
  geom_text(aes(x, y, label = lab, fontface = ifelse(bold, "bold", "plain")),
            size = txt_mm(FS - 0.6), family = FONT, colour = INK) +
  annotate("text", x = tcols$x, y = YL[2] + 1.0, label = tcols$head, vjust = 0, lineheight = 0.95,
           size = txt_mm(FS - 0.8), family = FONT, fontface = "bold", colour = MUTED) +
  coord_cartesian(xlim = c(0, 6.15), ylim = YL, expand = FALSE, clip = "off") +
  theme_void(base_family = FONT) +
  theme(plot.margin = margin(3, 3, 3, 2))

rowB <- pB + pBt + plot_layout(widths = c(1.3, 1))

## ---- 6. Panel C ---------------------------------------------------------------
GRP  <- c("ACC-I high", "ACC-I low", "ACC-II high", "ACC-II low")
GCOL <- setNames(unname(COL[c("I", "I", "II", "II")]), GRP)
GLTY <- setNames(c("solid", "22", "solid", "22"), GRP)
pC <- lapply(seq_len(nrow(cutt)), function(i) {
  r <- cutt[i, ]; g <- r$gene
  d$grp <- factor(paste(d$Subtype, ifelse(d[[g]] > r$cutoff, "high", "low")), GRP)
  keep  <- names(which(table(d$grp) >= 3))
  km <- km_steps(survfit(Surv(EFS_time, EFS_status) ~ grp, d[d$grp %in% keep, ]))
  km$strata <- factor(km$strata, GRP)
  ggplot(km, aes(time, surv, colour = strata, linetype = strata)) +
    geom_step(linewidth = 0.55) +
    geom_point(data = subset(km, cens), shape = 3, size = 0.9, stroke = 0.35, show.legend = FALSE) +
    scale_colour_manual(values = GCOL, drop = FALSE, name = NULL) +
    scale_linetype_manual(values = GLTY, drop = FALSE, name = NULL) +
    scale_x_continuous(breaks = seq(0, XMAX, 60)) +
    scale_y_continuous(breaks = seq(0, 1, 0.25)) +
    coord_cartesian(xlim = c(0, XMAX + 3), ylim = c(0, 1.02), expand = FALSE) +
    labs(title = if (i == 1) "C" else " ", x = "Months from diagnosis",
         y = if (i == 1) "EFS probability" else NULL,
         subtitle = sprintf(paste0("%s (high: > %.2f)\nSubtype-adjusted HR %s\n",
                                   "Stratified log-rank p = %s\nInteraction p = %s\n",
                                   "High/low: ACC-I %d/%d; ACC-II %d/%d"),
                            g, r$cutoff, fci(r$adj_HR, r$adj_lo, r$adj_hi),
                            fp(r$strat_logrank_p), fp(r$p_interaction),
                            r$hi_I, r$lo_I, r$hi_II, r$lo_II)) +
    theme_s3() + theme(legend.key.width = unit(7, "mm"))
})
rowC <- wrap_plots(pC, nrow = 1) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

## ---- 7. Assemble and save -----------------------------------------------------
fig <- wrap_elements(full = rowA) / wrap_elements(full = rowB) / wrap_elements(full = rowC) +
  plot_layout(heights = c(H_A, H_B, H_C))
pdf_file <- file.path(out_dir, "Supplementary_FigureS3_EFS_by_subtype.pdf")
ggsave(pdf_file, fig, width = PAGE_W, height = PAGE_H, units = "mm", device = cairo_pdf)

## Panel B geometry: vertical distance between adjacent estimates (mm)
pitch_mm <- (H_B - 45) / diff(YL)   # ~45 mm of row B are title, axis, legend
message(sprintf("Panel B: %d estimate rows, approx. %.1f mm between adjacent rows", length(EST) * nG, pitch_mm))
message("Saved: ", pdf_file)
