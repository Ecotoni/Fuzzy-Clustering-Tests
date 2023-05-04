#' Extension of vegan's cascadeKM to fuzzy clustering.
#'
#' @description This function is an extension of the [vegan::cascadeKM()]
#'     function from the vegan package for the case of fuzzy clustering.
#'     It internally uses the [vegclust::vegclust()] function, and 
#'     creates several partitions forming a cascade from a small 
#'     to a large number of groups.
#'
#' @param data The data matrix. The objects (samples) are the rows.
#' @param inf.gr The number of groups for the partition with the 
#'     smallest number of groups of the cascade (min).
#' @param sup.gr The number of groups for the partition with the 
#'     largest number of groups of the cascade (max).
#' @param iter The number of random starting configurations for each 
#'     value of k.
#' @param criterion The criterion that will be used to select the 
#'     fuzzy best partition. These are taken from the 
#'     [vegclust::vegclustIndex()] function of the vegclust package 
#'     (see Details)
#' @param cluster.method Clustering method to be used. Must be
#'     accepted by vegclust (see Details)
#' @param iter.max Maximum number of iterations for the 
#'     clustering algorithm.
#' @param parallel Number of parallel processes or a predefined 
#'     socket cluster. With parallel = 1 uses ordinary, 
#'     non-parallel processing. The parallel processing is 
#'     done with parallel package.
#' @param ... Additional arguments to be passed to the 
#'     [vegclust::vegclust()] function. Typically includes 
#'     the membership exponent m (defaulted to 2 by the 
#'     vegclust function).
#'
#' @details This function is essentially an extension of the 
#'     [vegan::cascadeKM()] function of the vegan package for fuzzy 
#'     clustering methods. Its functioning is kept identical. 
#'     The range of indices it accepts to select the best
#'     partition is, however, taken from fuzzy partition indices, 
#'     such as the ones coming from Fuzzy C-means (Bezdek 1981).
#'     It currently returns one of our values: partition 
#'     coefficient (PC), normalized partition coefficient (PCN), 
#'     partition entropy (PE) and normalized partition entropy (PEN).
#'     Maximum values of PCN or minimum values of PEN can be 
#'     used as criteria to choose the number of clusters.
#'     
#'     The function defaults to "FCM" for the fuzzy cluster algorithm,
#'     but can accept any of the algorithms implemented in the 
#'     original vegclust function, so long as they are compatible with
#'     the calculation of the fuzzy partition indices. 
#' 
#' @return Function cascadeFCM returns an object of class cascadeFCM 
#'     with items:
#'     * partition	
#'     Table with the partitions found for different numbers of groups K, 
#'     from K = inf.gr to K = sup.gr.
#'     * results	
#'     Values of the criterion to select the best partition.
#'     * criterion	
#'     The name of the criterion used.
#'     * size	
#'     The number of objects found in each group, for all 
#'     partitions (columns).
#'     
#' @export
#'
#' @examples
#'  # Partitioning a (10 x 10) data matrix of random numbers
#'  mat <- matrix(runif(100),10,10)
#'  res <- cascadeFCM(mat, 2, 10, iter = 25, criterion = "PC",
#'                    cluster.method = "FCM", m = 2)
#'  
cascadeFCM <- function (data, inf.gr, sup.gr, iter = 100, criterion = c("PC", "PCN", "PE", "PEN"), 
                        cluster.method = "FCM",
                        iter.max = 100, parallel = getOption("mc.cores"),
                        ...) 
{
  data <- as.matrix(data)
  if (!is.null(nrow(data))) {
    partition <- matrix(NA, nrow(data), sup.gr - inf.gr + 
                          1)
  }
  else {
    partition <- matrix(NA, length(data), sup.gr - inf.gr + 
                          1)
  }
  results <- matrix(NA, 2, sup.gr - inf.gr + 1)
  size <- matrix(NA, sup.gr, sup.gr - inf.gr + 1)
  h <- 1
  if (is.null(parallel)) 
    parallel <- 1
  hasClus <- inherits(parallel, "cluster")
  if (!hasClus && parallel <= 1) {
    tmp <- lapply(inf.gr:sup.gr, function(ii) {
      vegclust::vegclust(data, mobileCenters = ii, method = cluster.method, iter.max = iter.max,
                         ...)
    })
  }
  else {
    if (hasClus || .Platform$OS.type == "windows") {
      if (!hasClus) 
        cl <- makeCluster(parallel)
      tmp <- parLapply(cl, inf.gr:sup.gr, function(ii) vegclust::vegclust(data, mobileCenters = ii, 
                                                                          method = cluster.method, 
                                                                          iter.max = iter.max,
                                                                                ...))
      if (!hasClus) 
        stopCluster(cl)
    }
    else {
      tmp <- mclapply(inf.gr:sup.gr, function(ii) vegclust::vegclust(data, mobileCenters = ii, 
                                                                     method = cluster.method, 
                                                                     iter.max = iter.max,
                                                                     ...), mc.cores = parallel)
    }
  }
  for (ii in inf.gr:sup.gr) { ### need to check vegclust structure:
    idx <- ii - inf.gr + 1
    j <- ii - inf.gr + 1
    size[1:ii, h] <- tmp[[idx]]$size
    h <- h + 1
    partition[, j] <- apply(tmp[[idx]]$memb, MARGIN = 1, which.max)
    results[1, j] <- sum(tmp[[idx]]$withinss)
    results[2, j] <- vegclust::vegclustIndex(tmp[[idx]])[criterion]
  }
  colnames(partition) <- paste(inf.gr:sup.gr, "groups")
  tmp <- rownames(data)
  if (is.null(tmp)) {
    r.name <- c(1:nrow(partition))
  }
  else {
    r.name <- tmp
  }
  rownames(partition) <- r.name
  colnames(results) <- paste(inf.gr:sup.gr, "groups")
  rownames(results) <- c("SSE", criterion)
  colnames(size) <- paste(inf.gr:sup.gr, "groups")
  rownames(size) <- paste("Group", 1:sup.gr)
  tout <- list(partition = partition, results = results, criterion = criterion, 
               size = size)
  class(tout) <- "cascadeFCM"
  tout
}
