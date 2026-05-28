# install.packages("BiocManager")
# BiocManager::install("biomaRt")
library(biomaRt)
library(BayesMRclus)

data("whradjbmi_insulin", package = "BayesMRclus")
dat <- data.frame(SNP = whradjbmi_insulin[, "SNP"],
                  beta_exposure = whradjbmi_insulin[, "beta.exposure"],
                  beta_outcome = whradjbmi_insulin[, "beta.outcome"],
                  se_exposure = whradjbmi_insulin[, "se.exposure"],
                  se_outcome = whradjbmi_insulin[, "se.outcome"])

snp_names <- dat$SNP

mart <- useMart("ENSEMBL_MART_SNP",
                dataset  = "hsapiens_snp",
                # host     = "https://grch37.ensembl.org")  # GRCh37/hg19
                host     = "https://www.ensembl.org")
                # host     = "https://useast.ensembl.org")

batch_getBM <- function(snps, mart, batch_size = 30, pause = 2) {
  batches <- split(snps, ceiling(seq_along(snps) / batch_size))
  results <- vector("list", length(batches))
  for (i in seq_along(batches)) {
    cat(sprintf("Batch %d / %d ...\n", i, length(batches)))
    ok <- FALSE
    for (attempt in 1:3) {           # retry up to 3 times
      tryCatch({
        results[[i]] <- getBM(
          attributes = c("refsnp_id", "chr_name",
                         "chrom_start", "associated_gene"),
          filters    = "snp_filter",
          values     = batches[[i]],
          mart       = mart
        )
        ok <- TRUE
      }, error = function(e) {
        cat(sprintf("  Attempt %d failed: %s\n", attempt, e$message))
        Sys.sleep(pause * attempt)   # back off before retrying
      })
      if (ok) break
    }
    if (!ok) {
      warning(sprintf("Batch %d failed after 3 attempts; skipping.", i))
      results[[i]] <- data.frame()
    }
    Sys.sleep(pause)                 # be polite between batches
  }
  do.call(rbind, results)
}

httr::set_config(httr::timeout(120))
anno <- batch_getBM(snp_names, mart)
