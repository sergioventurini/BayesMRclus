# =============================================================================
# 0.  Internal helpers (not exported)
# =============================================================================

# ── Hungarian maximisation ────────────────────────────────────────────────────
# Returns assign_row_to_col: integer vector length n,
#   assign_row_to_col[i] = j  <=>  row i matched to column j.
.hungarian_max <- function(profit) {
  n <- nrow(profit)

  if (requireNamespace("clue", quietly = TRUE))
    return(as.integer(clue::solve_LSAP(profit, maximum = TRUE)))

  if (n <= 8 && requireNamespace("combinat", quietly = TRUE)) {
    perms    <- combinat::permn(seq_len(n))
    best_val <- -Inf
    best_p   <- seq_len(n)
    for (perm in perms) {
      v <- sum(profit[cbind(seq_len(n), unlist(perm))])
      if (v > best_val) { best_val <- v; best_p <- unlist(perm) }
    }
    return(best_p)
  }

  # Greedy fallback
  atrc <- integer(n); used <- logical(n)
  for (i in seq_len(n)) {
    avail    <- which(!used)
    j        <- avail[which.max(profit[i, avail])]
    atrc[i]  <- j
    used[j]  <- TRUE
  }
  atrc
}


# ── Display-label builders ────────────────────────────────────────────────────

# Build A-side display labels.
# Returns named character vector: A-integer (as string) -> display label.
.build_a_labels <- function(part_A, cluster_labels_A) {
  a_ints <- sort(unique(part_A))
  if (!is.null(cluster_labels_A)) {
    labs <- cluster_labels_A[as.character(a_ints)]
    if (anyNA(labs)) {
      warning("Some A integers not in cluster_labels_A; using integer labels.")
      labs[is.na(labs)] <- paste0("C", a_ints[is.na(labs)])
    }
  } else {
    labs <- paste0("C", a_ints)
  }
  setNames(as.character(labs), as.character(a_ints))
}

# Build B-side display labels from the matched partition.
# Recovers original B codes via cm$orig_B_of_matched, then applies
# cluster_labels_B if supplied.
# Returns named character vector: matched-integer (as string) -> display label.
.build_b_labels <- function(cm, cluster_labels_B) {
  b_ints_matched <- sort(unique(cm$part_B_matched))
  orig_codes     <- cm$orig_B_of_matched[as.character(b_ints_matched)]
  if (!is.null(cluster_labels_B)) {
    labs <- cluster_labels_B[as.character(orig_codes)]
    labs[is.na(labs)] <- paste0("C", orig_codes[is.na(labs)])
  } else {
    labs <- paste0("C", orig_codes)
  }
  setNames(as.character(labs), as.character(b_ints_matched))
}


# =============================================================================
# 1.  Optimal label matching
# =============================================================================

#' Build a contingency table and optimally match cluster labels
#'
#' Finds the permutation of B labels that maximises total overlap with A
#' (Hungarian algorithm on the contingency table used as a profit matrix).
#'
#' @param part_A  Integer vector length p
#' @param part_B  Integer vector length p
#' @return List:
#'   $table              contingency table, original labels (rows=A, cols=B)
#'   $table_matched      contingency table after relabelling B
#'   $part_B_matched     integer vector: B replaced by matched A-row indices
#'   $mapping            named integer vector:
#'                         names  = original B labels (character)
#'                         values = matched A-row index
#'   $orig_B_of_matched  named character vector:
#'                         names  = matched integer (character)
#'                         values = original B label (character)
#' @export
contingency_match <- function(part_A, part_B) {
  tab <- table(A = part_A, B = part_B)
  rA  <- nrow(tab);  rB <- ncol(tab);  n <- max(rA, rB)

  profit <- matrix(0L, n, n)
  profit[seq_len(rA), seq_len(rB)] <- tab

  atrc           <- .hungarian_max(profit)   # row i -> col atrc[i]
  actr           <- integer(n)
  actr[atrc]     <- seq_len(n)              # col j -> row actr[j]
  actr_B         <- actr[seq_len(rB)]       # length rB

  b_orig         <- colnames(tab)           # original B labels (character)
  mapping        <- setNames(actr_B, b_orig)
  stopifnot(!anyNA(mapping))

  orig_B_of_matched <- setNames(b_orig, as.character(actr_B))
  part_B_matched    <- as.integer(unname(mapping[as.character(part_B)]))

  list(
    table             = tab,
    table_matched     = table(A = part_A, B = part_B_matched),
    part_B_matched    = part_B_matched,
    mapping           = mapping,
    orig_B_of_matched = orig_B_of_matched
  )
}


# =============================================================================
# 2.  Numerical agreement measures
# =============================================================================

#' Adjusted Rand Index (ARI)
#'
#' Computed from scratch (no external package needed).
#' ARI = 1: perfect agreement; ARI = 0: chance; ARI < 0: worse than chance.
#'
#' @param part_A  Integer vector length p
#' @param part_B  Integer vector length p
#' @return  Scalar ARI
#' @export
ari <- function(part_A, part_B) {
  n      <- length(part_A)
  tab    <- table(part_A, part_B)
  sum_ab <- sum(choose(tab, 2))
  sum_a  <- sum(choose(rowSums(tab), 2))
  sum_b  <- sum(choose(colSums(tab), 2))
  expct  <- sum_a * sum_b / choose(n, 2)
  denom  <- (sum_a + sum_b) / 2 - expct
  if (denom == 0) return(1)
  (sum_ab - expct) / denom
}

#' Variation of Information (VI)
#'
#' \eqn{VI(A, B) = H(A|B) + H(B|A)}, where H denotes conditional entropy.
#' Ranges from 0 (perfect agreement) to \eqn{log(p)} (maximum disagreement).
#' Optionally normalised to \eqn{[0,1]} by dividing by \eqn{log(p)}.
#'
#' @param part_A    Integer vector length p
#' @param part_B    Integer vector length p
#' @param normalise Logical; if TRUE divide by log(p) (default TRUE)
#' @return  Scalar VI (or normalised VI)
#' @export
vi <- function(part_A, part_B, normalise = TRUE) {
  n    <- length(part_A)
  tab  <- table(part_A, part_B) / n
  p_A  <- rowSums(tab);  p_B <- colSums(tab)
  H_A  <- -sum(p_A[p_A > 0] * log(p_A[p_A > 0]))
  H_B  <- -sum(p_B[p_B > 0] * log(p_B[p_B > 0]))
  nz   <- tab > 0
  I_AB <- sum(tab[nz] * log(tab[nz] / outer(p_A, p_B)[nz]))
  vi_val <- H_A + H_B - 2 * I_AB
  if (normalise) {
    H_AB   <- -sum(tab[nz] * log(tab[nz]))
    # vi_val <- vi_val / log(n)  # alternative (less standard) normalization
    vi_val <- if (H_AB > 0) vi_val / H_AB else 0
  }
  vi_val
}

#' Fowlkes-Mallows Index (FMI)
#'
#' Geometric mean of pairwise precision and recall.
#' FMI = 1: perfect agreement; FMI = 0: no agreement.
#'
#' @param part_A  Integer vector length p
#' @param part_B  Integer vector length p
#' @return  Scalar FMI
#' @export
fmi <- function(part_A, part_B) {
  tab <- table(part_A, part_B)
  TP  <- sum(choose(tab, 2))
  if (TP == 0) return(0)
  FP  <- sum(choose(colSums(tab), 2)) - TP
  FN  <- sum(choose(rowSums(tab), 2)) - TP
  sqrt(TP / (TP + FP) * TP / (TP + FN))
}

#' Print all three agreement measures and the contingency tables
#'
#' @param part_A   Integer vector length p
#' @param part_B   Integer vector length p
#' @param name_A   Method A label (default "MRClust")
#' @param name_B   Method B label (default "BNPM-MR")
#' @param verbose  Print to console (default FALSE)
#' @return  Invisibly: list($measures, $contingency)
#' @export
cluster_agreement_summary <- function(part_A, part_B,
                                      name_A  = "MRClust",
                                      name_B  = "BNPM-MR",
                                      verbose = FALSE) {
  cm <- contingency_match(part_A, part_B)
  df <- data.frame(
    measure = c("ARI", "NVI", "FMI", "K_A", "K_B"),
    value   = c(round(ari(part_A, part_B), 4),
                round(vi( part_A, part_B), 4),
                round(fmi(part_A, part_B), 4),
                length(unique(part_A)),
                length(unique(part_B))),
    note    = c("Adjusted Rand Index    (1=perfect, 0=chance)",
                "Normalised VI          (0=perfect)",
                "Fowlkes-Mallows Index  (1=perfect)",
                paste("Clusters in", name_A),
                paste("Clusters in", name_B)),
    stringsAsFactors = FALSE
  )
  if (verbose) {
    cat(sprintf("\n=== Cluster agreement: %s vs. %s ===\n", name_A, name_B))
    print(df, row.names = FALSE)
    cat(sprintf("\nContingency table (rows=%s, cols=%s, original labels):\n",
                name_A, name_B))
    print(cm$table)
    cat(sprintf("\nContingency table (%s columns optimally relabelled):\n", name_B))
    print(cm$table_matched)
  }
  invisible(list(measures = df, contingency = cm))
}


# =============================================================================
# 3.  Alluvial diagram
# =============================================================================

#' Alluvial (Sankey) diagram: SNP flow between two partitions
#'
#' Each SNP is a flow. Left column = method A clusters; right = method B.
#' Width of each ribbon is proportional to the number of SNPs in that
#' (A-cluster, B-cluster) combination. SNPs are coloured by their A-cluster
#' assignment, so you can immediately see which A-clusters are "split" or
#' "merged" by method B.
#'
#' Requires: ggplot2, ggalluvial, dplyr
#'
#' @param part_A           Integer vector length p
#' @param part_B           Integer vector length p
#' @param snp_names        Optional character vector length p
#' @param name_A           Method A label
#' @param name_B           Method B label
#' @param cluster_labels_A Named character: A integer -> description (optional)
#' @param cluster_labels_B Named character: B integer -> description (optional)
#' @param file             PNG path or NULL (returns ggplot object)
plot_alluvial <- function(part_A, part_B,
                          snp_names        = NULL,
                          name_A           = "MRClust",
                          name_B           = "BNPM-MR",
                          cluster_labels_A = NULL,
                          cluster_labels_B = NULL,
                          file             = NULL) {
  for (pkg in c("ggplot2", "ggalluvial", "dplyr"))
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(sprintf("Package '%s' required.", pkg))
  library(ggplot2); library(ggalluvial); library(dplyr)

  p  <- length(part_A)
  cm <- contingency_match(part_A, part_B)

  a_map <- .build_a_labels(part_A, cluster_labels_A)
  b_map <- .build_b_labels(cm, cluster_labels_B)

  lab_A <- factor(a_map[as.character(part_A)],          levels = a_map)
  lab_B <- factor(b_map[as.character(cm$part_B_matched)], levels = b_map)

  df     <- data.frame(
    snp = if (is.null(snp_names)) seq_len(p) else snp_names,
    A   = lab_A,
    B   = lab_B,
    stringsAsFactors = FALSE
  )
  df_agg <- df %>% dplyr::count(A, B, name = "freq")

  pal_A <- setNames(
    rep(c("#D62728","#1F77B4","#2CA02C","#FF7F0E",
          "#9467BD","#8C564B","#E377C2","#7F7F7F"),
        length.out = length(levels(lab_A))),
    levels(lab_A)
  )

  g <- ggplot(df_agg, aes(axis1 = A, axis2 = B, y = freq)) +
    geom_alluvium(aes(fill = A), width = 1/5, alpha = 0.75, knot.pos = 0.4) +
    geom_stratum(width = 1/5, fill = "grey90", color = "grey40") +
    geom_text(stat = "stratum",
              aes(label = after_stat(stratum)),
              size = 3.2, fontface = "bold") +
    scale_x_discrete(limits = c(name_A, name_B), expand = c(0.12, 0.12)) +
    scale_fill_manual(values = pal_A) +
    labs(title = sprintf("SNP cluster assignment: %s vs. %s", name_A, name_B),
         y = "Number of SNPs", fill = name_A) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right", panel.grid.major.x = element_blank())

  if (!is.null(file)) { ggsave(file, g, width = 7, height = 6, dpi = 130); invisible(g) }
  else g
}


# =============================================================================
# 4.  Wald-ratio scatter, two panels
# =============================================================================

#' Scatter plot: Wald ratios coloured by cluster, two panels
#'
#' Cluster colours are harmonised across panels via optimal matching.
#' Discordant SNPs (different matched cluster) are circled in black.
#'
#' The discord flag is computed by comparing matched integer indices
#' (part_A vs. part_B_matched), NOT by comparing display strings.
#' This eliminates the NA-panel bug caused by string comparison failures.
#'
#' @param part_A              Integer vector length p
#' @param part_B              Integer vector length p
#' @param wald                Numeric vector length p
#' @param se_wald             Numeric vector length p (|sy / bx|)
#' @param name_A              Method A label
#' @param name_B              Method B label
#' @param cluster_labels_A    Named character: A integer -> description (optional)
#' @param cluster_labels_B    Named character: B integer -> description (optional)
#' @param highlight_discordant  Circle discordant SNPs (default TRUE)
#' @param file                PNG path or NULL
plot_wald_clusters <- function(part_A, part_B,
                               wald, se_wald,
                               name_A               = "MRClust",
                               name_B               = "BNPM-MR",
                               cluster_labels_A     = NULL,
                               cluster_labels_B     = NULL,
                               highlight_discordant = TRUE,
                               file                 = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required.")
  library(ggplot2)

  p  <- length(part_A)
  cm <- contingency_match(part_A, part_B)

  a_map <- .build_a_labels(part_A, cluster_labels_A)
  b_map <- .build_b_labels(cm, cluster_labels_B)

  lab_A_str <- a_map[as.character(part_A)]
  lab_B_str <- b_map[as.character(cm$part_B_matched)]

  # Unified colour palette (union of both label sets)
  all_levels <- union(as.character(a_map), as.character(b_map))
  pal <- setNames(
    rep(c("#D62728","#1F77B4","#2CA02C","#FF7F0E",
          "#9467BD","#8C564B","#E377C2","#7F7F7F"),
        length.out = length(all_levels)),
    all_levels
  )

  # Discord: compare matched integer indices — never produces NA
  discordant <- (part_A != cm$part_B_matched)

  df_A <- data.frame(wald    = wald,
                     prec    = 1 / se_wald,
                     cluster = factor(lab_A_str, levels = all_levels),
                     method  = name_A,
                     discord = discordant,
                     stringsAsFactors = FALSE)
  df_B <- data.frame(wald    = wald,
                     prec    = 1 / se_wald,
                     cluster = factor(lab_B_str, levels = all_levels),
                     method  = name_B,
                     discord = discordant,
                     stringsAsFactors = FALSE)

  df        <- rbind(df_A, df_B)
  df$method <- factor(df$method, levels = c(name_A, name_B))

  g <- ggplot(df, aes(x = wald, y = prec, colour = cluster)) +
    geom_point(size = 1.2, alpha = 0.75) +
    facet_wrap(~ method, ncol = 2) +
    scale_colour_manual(values = pal, drop = TRUE, na.value = "grey50") +
    labs(x      = "Wald ratio",
         y      = "Precision (1 / SE)",
         colour = "Cluster",
         # title  = sprintf("Cluster assignments: %s vs. %s", name_A, name_B)) +
         title  = "") +
    theme_bw(base_size = 11) +
    theme(strip.text = element_text(face = "bold"))

  if (highlight_discordant) {
    df_disc <- df[df$discord, ]          # discord is always TRUE/FALSE, never NA
    g <- g +
      geom_point(data = df_disc,
                 aes(x = wald, y = prec),
                 shape = 1, size = 3.2, colour = "black",
                 stroke = 0.6, inherit.aes = FALSE) +
      labs(caption = "Open circle: SNP assigned to different cluster by the two methods")
  }

  if (!is.null(file)) { ggsave(file, g, width = 9, height = 5, dpi = 130); invisible(g) }
  else g
}


# =============================================================================
# 5.  Contingency heatmap
# =============================================================================

#' Heatmap of the (optimally matched) contingency table
#'
#' Row axis = method A (own labels).
#' Column axis = method B (own labels, recovered from cm$orig_B_of_matched).
#' A perfect match appears as a filled diagonal.
#'
#' @param part_A           Integer vector length p
#' @param part_B           Integer vector length p
#' @param name_A           Method A label
#' @param name_B           Method B label
#' @param cluster_labels_A Named character: A integer -> description (optional)
#' @param cluster_labels_B Named character: B integer -> description (optional)
#' @param file             PNG path or NULL
plot_contingency_heatmap <- function(part_A, part_B,
                                     name_A           = "MRClust",
                                     name_B           = "BNPM-MR",
                                     cluster_labels_A = NULL,
                                     cluster_labels_B = NULL,
                                     file             = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required.")
  library(ggplot2)

  cm    <- contingency_match(part_A, part_B)
  a_map <- .build_a_labels(part_A, cluster_labels_A)
  b_map <- .build_b_labels(cm, cluster_labels_B)

  df        <- as.data.frame(cm$table_matched)
  names(df) <- c("A", "B", "n")

  df$A <- factor(a_map[as.character(df$A)], levels = rev(as.character(a_map)))
  df$B <- factor(b_map[as.character(df$B)], levels = as.character(b_map))

  g <- ggplot(df, aes(x = B, y = A, fill = n)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(n > 0, as.character(n), "")),
              size = 3.5, fontface = "bold", colour = "white") +
    scale_fill_gradient(low = "#deebf7", high = "#08519c", name = "SNPs") +
    labs(x        = name_B,
         y        = name_A,
         title    = sprintf("Contingency table: %s vs. %s", name_A, name_B),
         subtitle = sprintf("(%s columns relabelled to maximize diagonal)",
                            name_B)) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          axis.text.y = element_text(hjust = 1),
          panel.grid  = element_blank())

  if (!is.null(file)) { ggsave(file, g, width = 6, height = 5, dpi = 130); invisible(g) }
  else g
}


# =============================================================================
# 6.  Per-SNP concordance table
# =============================================================================

#' Per-SNP cluster assignment and concordance flag
#'
#' A SNP is concordant when its matched B index equals its A index.
#' Concordance is always a clean TRUE/FALSE (no NAs).
#'
#' @param part_A           Integer vector length p
#' @param part_B           Integer vector length p
#' @param snp_names        Optional character vector length p
#' @param wald             Optional numeric vector length p
#' @param name_A           Method A label
#' @param name_B           Method B label
#' @param cluster_labels_A Named character: A integer -> description (optional)
#' @param cluster_labels_B Named character: B integer -> description (optional)
#' @return  Data frame (sorted by Wald ratio if supplied)
snp_agreement_table <- function(part_A, part_B,
                                snp_names        = NULL,
                                wald             = NULL,
                                name_A           = "MRClust",
                                name_B           = "BNPM-MR",
                                cluster_labels_A = NULL,
                                cluster_labels_B = NULL) {
  p  <- length(part_A)
  cm <- contingency_match(part_A, part_B)

  a_map <- .build_a_labels(part_A, cluster_labels_A)
  b_map <- .build_b_labels(cm, cluster_labels_B)

  lab_A  <- a_map[as.character(part_A)]
  lab_B  <- b_map[as.character(cm$part_B_matched)]
  agree  <- (part_A == cm$part_B_matched)   # always TRUE/FALSE

  df <- data.frame(
    snp   = if (is.null(snp_names)) seq_len(p) else snp_names,
    wald  = if (is.null(wald))      NA_real_   else round(wald, 4),
    A     = as.character(lab_A),
    B     = as.character(lab_B),
    agree = agree,
    stringsAsFactors = FALSE
  )
  names(df)[3:4] <- c(name_A, name_B)
  if (!is.null(wald)) df <- df[order(wald), ]
  df
}


# =============================================================================
# 7.  Full pipeline
# =============================================================================

#' Run the complete cluster comparison pipeline
#'
#' @param part_mrclust           Integer vector length p
#' @param part_bnpmr             Integer vector length p
#' @param snp_names              Optional character vector length p
#' @param wald                   Numeric vector length p (Wald ratios)
#' @param se_wald                Numeric vector length p (SE of Wald ratios)
#' @param cluster_labels_mrclust Named character: MRClust integer -> description
#'   e.g. c("1"="Positive","2"="Negative","3"="Null","4"="Junk")
#' @param cluster_labels_bnpmr   Named character: BNPM-MR integer -> description
#'   (optional; clusters labelled "C1", "C2", … if NULL)
#' @param fig_dir   Directory for PNG output (NULL = display only)
#' @param csv_dir   Directory for CSV output (NULL = skip)
#' @param verbose   Print numerical results (default FALSE)
#' @return  Named list (invisibly)
#' @export
compare_partitions <- function(part_mrclust,
                               part_bnpmr,
                               snp_names              = NULL,
                               wald                   = NULL,
                               se_wald                = NULL,
                               cluster_labels_mrclust = NULL,
                               cluster_labels_bnpmr   = NULL,
                               fig_dir                = NULL,
                               csv_dir                = NULL,
                               verbose                = FALSE) {
  if (!is.null(fig_dir)) dir.create(fig_dir, showWarnings = FALSE)
  if (!is.null(csv_dir)) dir.create(csv_dir, showWarnings = FALSE)

  agree <- cluster_agreement_summary(
    part_mrclust, part_bnpmr,
    name_A = "MRClust", name_B = "BNPM-MR", verbose = verbose
  )

  snp_tab <- snp_agreement_table(
    part_mrclust, part_bnpmr,
    snp_names        = snp_names,
    wald             = wald,
    cluster_labels_A = cluster_labels_mrclust,
    cluster_labels_B = cluster_labels_bnpmr
  )
  if (verbose) {
    cat(sprintf("\nSNPs concordant: %d / %d (%.1f%%)\n",
                sum(snp_tab$agree), nrow(snp_tab),
                100 * mean(snp_tab$agree)))
    cat("\nDiscordant SNPs:\n")
    print(snp_tab[!snp_tab$agree, ], row.names = FALSE)
  }
  if (!is.null(csv_dir))
    write.csv(snp_tab, file.path(csv_dir, "snp_cluster_agreement.csv"),
              row.names = FALSE)

  f_all <- if (!is.null(fig_dir)) file.path(fig_dir, "fig_alluvial.png") else NULL
  fig_alluvial <- plot_alluvial(
    part_mrclust, part_bnpmr,
    snp_names        = snp_names,
    cluster_labels_A = cluster_labels_mrclust,
    cluster_labels_B = cluster_labels_bnpmr,
    file             = f_all
  )
  if (is.null(f_all)) print(fig_alluvial)

  fig_scatter <- NULL
  if (!is.null(wald) && !is.null(se_wald)) {
    f_sc <- if (!is.null(fig_dir)) file.path(fig_dir, "fig_cluster_scatter.png") else NULL
    fig_scatter <- plot_wald_clusters(
      part_mrclust, part_bnpmr,
      wald                 = wald,
      se_wald              = se_wald,
      cluster_labels_A     = cluster_labels_mrclust,
      cluster_labels_B     = cluster_labels_bnpmr,
      highlight_discordant = TRUE,
      file                 = f_sc
    )
    if (is.null(f_sc)) print(fig_scatter)
  }

  f_hm <- if (!is.null(fig_dir)) file.path(fig_dir, "fig_contingency_heatmap.png") else NULL
  fig_heatmap <- plot_contingency_heatmap(
    part_mrclust, part_bnpmr,
    cluster_labels_A = cluster_labels_mrclust,
    cluster_labels_B = cluster_labels_bnpmr,
    file             = f_hm
  )
  if (is.null(f_hm)) print(fig_heatmap)

  invisible(list(
    agreement    = agree,
    snp_table    = snp_tab,
    fig_alluvial = fig_alluvial,
    fig_scatter  = fig_scatter,
    fig_heatmap  = fig_heatmap
  ))
}
