############################################################
# ACC - ROC ANALYSIS: TUMOR vs NORMAL
#
# Genes:
# VTCN1, MYB, NECTIN4, TACSTD2, STING1, KIT,
# FAP, RBPJ, CDK4, CDK6, EGFR, MET
#
# OUTPUT:
# 1. ACC_ROC_panel.pdf
# 2. ACC_ROC_AUC_results.csv
############################################################


############################################################
# 1. PACKAGES
############################################################

if (!requireNamespace("pROC", quietly = TRUE)) {
  install.packages("pROC")
}

library(pROC)


############################################################
# 2. IMPORT DATA
############################################################
data <- read.delim(
  "ACC292_genes.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)


############################################################
# 3. GENES TO ANALYSE
############################################################

genes <- c(
  "VTCN1",
  "MYB",
  "NECTIN4",
  "TACSTD2",
  "STING1",
  "KIT",
  "FAP",
  "RBPJ",
  "CDK4",
  "CDK6",
  "EGFR",
  "MET"
)


############################################################
# 4. DEFINE TUMOR / NORMAL GROUP
############################################################

data$Group <- factor(
  data[["Covariate 1"]],
  levels = c("norm", "tum")
)

# Check sample numbers
print(table(data$Group))


############################################################
# 5. ROC ANALYSIS
############################################################

results <- data.frame()

roc_list <- list()


for (gene in genes) {

  # Select samples without missing values
  tmp <- data[
    !is.na(data[[gene]]) &
    !is.na(data$Group),
    c(gene, "Group")
  ]

  ##########################################################
  # ROC curve
  ##########################################################

  roc_obj <- roc(
    response = tmp$Group,
    predictor = tmp[[gene]],
    levels = c("norm", "tum"),
    direction = "auto",
    quiet = TRUE
  )

  roc_list[[gene]] <- roc_obj


  ##########################################################
  # AUC
  ##########################################################

  auc_value <- as.numeric(
    auc(roc_obj)
  )


  ##########################################################
  # 95% CI OF AUC - DeLong method
  ##########################################################

  auc_ci <- ci.auc(
    roc_obj,
    method = "delong"
  )


  ##########################################################
  # P-VALUE: Wilcoxon rank-sum test
  ##########################################################

  wtest <- wilcox.test(
    tmp[[gene]] ~ tmp$Group,
    exact = FALSE
  )


  ##########################################################
  # MEDIANS
  ##########################################################

  median_norm <- median(
    tmp[tmp$Group == "norm", gene],
    na.rm = TRUE
  )

  median_tum <- median(
    tmp[tmp$Group == "tum", gene],
    na.rm = TRUE
  )


  ##########################################################
  # EXPRESSION DIRECTION
  ##########################################################

  expression_direction <- ifelse(
    median_tum > median_norm,
    "Higher in tumor",
    "Lower in tumor"
  )


  ##########################################################
  # SAVE RESULTS
  ##########################################################

  results <- rbind(
    results,
    data.frame(
      Gene = gene,
      AUC = auc_value,
      CI_low = as.numeric(auc_ci[1]),
      CI_high = as.numeric(auc_ci[3]),
      P_value = wtest$p.value,
      Median_normal = median_norm,
      Median_tumor = median_tum,
      Direction = expression_direction
    )
  )
}


############################################################
# 6. FDR CORRECTION
############################################################

results$FDR_BH <- p.adjust(
  results$P_value,
  method = "BH"
)


############################################################
# 7. FORMAT RESULTS
############################################################

results$AUC <- round(
  results$AUC,
  3
)

results$CI_low <- round(
  results$CI_low,
  3
)

results$CI_high <- round(
  results$CI_high,
  3
)

results$Median_normal <- round(
  results$Median_normal,
  3
)

results$Median_tumor <- round(
  results$Median_tumor,
  3
)


############################################################
# CREATE 95% CI COLUMN
############################################################

results$CI_95 <- paste0(
  sprintf("%.3f", results$CI_low),
  " - ",
  sprintf("%.3f", results$CI_high)
)


############################################################
# FINAL TABLE
############################################################

final_table <- results[, c(
  "Gene",
  "AUC",
  "CI_95",
  "P_value",
  "FDR_BH",
  "Median_normal",
  "Median_tumor",
  "Direction"
)]


############################################################
# ORDER BY AUC
############################################################

final_table <- final_table[
  order(
    final_table$AUC,
    decreasing = TRUE
  ),
]


############################################################
# PRINT TABLE
############################################################

print(final_table)


############################################################
# 8. EXPORT TABLE
############################################################

write.csv(
  final_table,
  file = "ACC_ROC_AUC_results.csv",
  row.names = FALSE
)

write.table(
  final_table,
  file = "ACC_ROC_AUC_results.txt",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


############################################################
# 9. ROC PANEL IN PDF
############################################################

pdf(
  file = "ACC_ROC_panel.pdf",
  width = 12,
  height = 9,
  onefile = TRUE
)


############################################################
# PANEL: 3 ROWS x 4 COLUMNS
############################################################

par(
  mfrow = c(3, 4),
  mar = c(4.5, 4.5, 3.5, 1.5),
  mgp = c(2.5, 0.8, 0),
  las = 1
)


############################################################
# PLOT EACH ROC
############################################################

for (gene in genes) {

  roc_obj <- roc_list[[gene]]

  row_gene <- results[
    results$Gene == gene,
  ]


  ##########################################################
  # ROC PLOT - BLUE CURVE
  ##########################################################

  plot(
    roc_obj,
    legacy.axes = TRUE,

    col = "blue",
    lwd = 2.5,

    main = gene,

    xlab = "1 - Specificity",
    ylab = "Sensitivity",

    cex.main = 1.2,
    cex.lab = 1.0,
    cex.axis = 0.9,

    print.auc = FALSE
  )


  ##########################################################
  # REFERENCE DIAGONAL
  ##########################################################

  abline(
    a = 0,
    b = 1,
    lty = 2,
    lwd = 1,
    col = "grey60"
  )


  ##########################################################
  # P-VALUE LABEL
  ##########################################################

  if (row_gene$P_value < 0.001) {

    p_label <- "p < 0.001"

  } else {

    p_label <- paste0(
      "p = ",
      sprintf(
        "%.3f",
        row_gene$P_value
      )
    )
  }


  ##########################################################
  # TEXT: AUC + CI + P-VALUE
  ##########################################################

  legend(
    "bottomright",

    legend = c(

      paste0(
        "AUC = ",
        sprintf(
          "%.2f",
          row_gene$AUC
        )
      ),

      paste0(
        "95% CI: ",
        sprintf(
          "%.2f",
          row_gene$CI_low
        ),
        "-",
        sprintf(
          "%.2f",
          row_gene$CI_high
        )
      ),

      p_label
    ),

    bty = "n",
    cex = 0.82
  )
}


############################################################
# CLOSE PDF
############################################################

dev.off()


############################################################
# 10. MESSAGE
############################################################

cat(
  "\nAnalysis completed.\n",
  "ROC panel saved as: ACC_ROC_panel.pdf\n",
  "Results saved as: ACC_ROC_AUC_results.csv\n"
)
