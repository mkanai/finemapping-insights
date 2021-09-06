library(doParallel)
library(dplyr)
library(ggplot2)
library(patchwork)
library(Matrix)
source("~/src/github.com/mkanai/finemapping-insights/R/const.R")

jaccard <- function(input_m, n_cores) {
  # Non zero values
  A <- tcrossprod(t(input_m))
  # Indices for non-zero values
  im <- which(as.matrix(A) > 0, arr.ind = TRUE)
  # Weighted jaccard similarity
  # wj <- apply(im, 1, function(x) {
  #  # Get min/max of dense matrix (sparse function was slower)
  #  x2 <- rowRanges(as.matrix(input_m[, x]), cols = c(1, 2))
  #  # Get weighted jaccard
  #  s <- colSums2(x2)
  #  s[1] / s[2]
  # })
  # Setup parallelization
  cl <- parallel::makeForkCluster(n_cores)
  doParallel::registerDoParallel(cl)
  # Get block indices for each core
  stop <- dim(im)[1]
  if (stop < n_cores) {
    indices <- rep(1, stop)
    n_cores <- 1
  } else {
    indices <- as.numeric(cut(1:stop, n_cores))
  }
  # Weighted jaccard similarity
  wj <- foreach(i = 1:n_cores, .combine = "c") %dopar% {
    # Subset non-zero indice for core
    im_loop <- im[indices == i, , drop = F]
    X <- apply(im_loop, 1, function(x) {
      # Get min/max of dense matrix (sparse function was slower)
      x2 <- matrixStats::rowRanges(as.matrix(input_m[, x, drop = F]), cols = c(1, 2))
      # Get weighted jaccard
      s <- matrixStats::colSums2(x2)
      s[1] / s[2]
    })
  }
  # Stop parallelization
  parallel::stopCluster(cl)
  # Form sparse matrix of weighted jaccard similarities
  J <- sparseMatrix(
    i = im[, 1],
    j = im[, 2],
    x = wj,
    dims = dim(A),
    dimnames = dimnames(A)
  )
  return(J)
}

find_overlaps <- function(CS.HC.df, n_cores = 4, plot_path = NULL) {
  if (length(unique(CS.HC.df$id)) == 1) {
    return(dplyr::mutate(CS.HC.df, csm_id = 1))
  }

  # Convert to sparse matrix
  CS.HC.mat <- CS.HC.df %>%
    tidytext::cast_sparse(variant, id, pip) %>%
    as(., "dgCMatrix")
  # Compute weighted Jaccard similarity matrix
  system.time(CS.HC.jsm <- jaccard(CS.HC.mat, n_cores))
  # Get dissimilarity matrix
  CS.HC.dist <- as.dist(1 - CS.HC.jsm)
  # Hierarchical clustering using Complete Linkage
  CS.HC.hc1 <- hclust(CS.HC.dist, method = "complete")
  # Plot wj minimum cutoff choices
  x <- seq(0, 1, 0.005)
  CS.HC.hc1.cut <- cutree(CS.HC.hc1, h = x)
  y <- apply(CS.HC.hc1.cut, 2, function(x) {
    length(unique(x))
  })
  if (!is.null(plot_path)) {
    cairo_pdf(plot_path, width = 5, height = 5)
    plot(x, y, pch = 18, ylab = "number of clusters", xlab = "max distance in cluster")
    dev.off()
    write.table(data.frame(cutoff = x, n_clusters = y), stringr::str_replace(plot_path, ".pdf$", ".tsv"), quote = F, row.names = F, sep = "\t")
  }

  # Given cutoff, assign IDs to 95% CSs
  CS.HC.hc1.df <- cutree(CS.HC.hc1, h = 0.9) %>%
    data.frame() %>%
    tibble::rownames_to_column() %>%
    tibble::as_tibble() %>%
    dplyr::rename("id" = "rowname", "csm_id" = ".")
  CS.HC.df <- CS.HC.df %>%
    merge(., CS.HC.hc1.df, by = "id") %>%
    tibble::as_tibble() %>%
    tidyr::separate(id,
      into = c("pop", "trait", "region", "cs_id"),
      sep = ";",
      remove = F
    )
  return(CS.HC.df)
}

df.cs_id <- rgsutil::read_gsfile(fm_insights$get_merged_results_path("fm_only.cs_id", "tsv.bgz")) %>%
  dplyr::rename(pip = susie.pip)

df.csm <-
  dplyr::mutate(df.cs_id, id = paste(cohort, trait, region, susie.cs_id, sep = ";")) %>%
  dplyr::group_split(trait) %>%
  purrr::map_dfr(~ {
    print(.$trait[1])
    find_overlaps(., n_cores = 8, plot_path = sprintf("./figures/cs_hclust/%s.hclust.pdf", .$trait[1]))
  })

rgsutil::write_gsfile(df.csm, fm_insights$get_merged_results_path("fm_only.csm_id", "tsv.bgz"), overwrite = TRUE)


df.csm <- rgsutil::read_gsfile(fm_insights$get_merged_results_path("fm_only.csm_id", "tsv.bgz"))

# unique merged CS-trait pairs
dplyr::select(df.csm, trait, csm_id) %>%
  dplyr::distinct() %>%
  nrow()

# median CS size
dplyr::group_by(df.csm, trait, csm_id) %>%
  dplyr::summarize(count = n()) %>%
  .$count %>%
  median()

# PIP > 0.1 CS
dplyr::group_by(df.csm, trait, csm_id) %>%
  dplyr::summarize(pip01 = max(pip, na.rm = T) > 0.1) %>%
  .$pip01 %>%
  sum()

# X CS
dplyr::filter(df.csm, stringr::str_starts(variant, "X:")) %>%
  dplyr::select(trait, csm_id) %>%
  dplyr::distinct() %>%
  nrow()