suppressPackageStartupMessages({
  library("DESeq2")
  library("data.table")
  library("ggplot2")
})

args <- commandArgs(trailingOnly = TRUE)
getArg <- function(flag) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) stop(paste("Missing", flag))
  args[[idx + 1]]
}

counts_file <- getArg("--counts")
samples_file <- getArg("--samples")
outdir <- getArg("--outdir")

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

condition_col <- "{{CONDITION_COL}}"
alpha <- as.numeric("{{ALPHA}}")
lfc_threshold <- as.numeric("{{LFC_THRESHOLD}}")
contrast_enabled <- "{{CONTRAST_ENABLED}}"
contrast_ref <- "{{CONTRAST_REF}}"
contrast_test <- "{{CONTRAST_TEST}}"

samples <- fread(samples_file)
if (!("sample" %in% names(samples))) stop("samples.csv must contain column: sample")
if (!(condition_col %in% names(samples))) stop(paste("samples.csv must contain column:", condition_col))

samples$sample <- as.character(samples$sample)
samples[[condition_col]] <- as.factor(samples[[condition_col]])

if (contrast_enabled == "1") {
  samples[[condition_col]] <- relevel(samples[[condition_col]], ref=contrast_ref)
}

fc <- fread(counts_file)
if (!("Geneid" %in% names(fc))) stop("Counts TSV must contain 'Geneid' column.")

common <- intersect(names(fc), samples$sample)
if (length(common) < 2) stop("Need >=2 sample columns in counts that match samples.csv sample values.")

counts <- as.data.frame(fc[, ..common])
rownames(counts) <- fc$Geneid

for (nm in colnames(counts)) {
  counts[[nm]] <- as.integer(round(as.numeric(counts[[nm]])))
}
counts[] <- lapply(counts, function(x) { x[is.na(x)] <- 0L; x })

coldata <- as.data.frame(samples[match(colnames(counts), samples$sample), ])
rownames(coldata) <- coldata$sample

dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(counts),
  colData = coldata,
  design = as.formula(paste0("~ ", condition_col))
)
dds <- dds[rowSums(counts(dds)) > 0,]
dds <- DESeq(dds)

if (contrast_enabled == "1") {
  contrast <- c(condition_col, contrast_test, contrast_ref)
  res <- results(dds, contrast=contrast, alpha=alpha, lfcThreshold=lfc_threshold)
} else {
  res <- results(dds, alpha=alpha, lfcThreshold=lfc_threshold)
}

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[order(is.na(res_df$padj), res_df$padj, is.na(res_df$pvalue), res_df$pvalue), ]
fwrite(res_df, file.path(outdir, "deseq_results.tsv"), sep="\t", na="NA")

sink(file.path(outdir, "session_info.txt"))
sessionInfo()
sink()

png(file.path(outdir, "MA_plot.png"), width=1200, height=900, res=150)
plotMA(res, ylim=c(-5, 5))
dev.off()

vsd <- vst(dds, blind=TRUE)
p <- plotPCA(vsd, intgroup=c(condition_col), returnData=TRUE)
percentVar <- round(100 * attr(p, "percentVar"))

p_pca <- ggplot(p, aes(x=PC1, y=PC2, color=.data[[condition_col]])) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "%")) +
  ylab(paste0("PC2: ", percentVar[2], "%")) +
  theme_minimal()

ggsave(file.path(outdir, "pca_plot.png"), plot=p_pca, width=7, height=5, dpi=150)

vol <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), ]
vol$neglog10padj <- -log10(vol$padj)

p_vol <- ggplot(vol, aes(x=log2FoldChange, y=neglog10padj)) +
  geom_point(alpha=0.6) +
  theme_minimal() +
  xlab("log2 fold change") +
  ylab("-log10(adjusted p-value)")

ggsave(file.path(outdir, "volcano_plot.png"), plot=p_vol, width=7, height=5, dpi=150)

cat("Done. Wrote results to:", outdir, "\n")
