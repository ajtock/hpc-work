#!/applications/R/R-4.0.0/bin/Rscript

# Set class for permutation test results object
setClass("permTest",
         representation(alternative = "character",
                        alphaThreshold = "numeric",
                        pval = "numeric",
                        observed = "numeric",
                        expected = "numeric",
                        log2obsexp = "numeric",
                        log2alpha = "numeric",
                        permDist = "numeric",
                        features = "numeric",
                        fam = "character",
                        accession = "character",
                        metric = "character",
                        region = "character"))

# Permutation test function to evaluate if CEN180 metric_name is significantly
# higher or lower at CENATHILA than at perms sets of random centromeric loci
permTest <- function(acc_idx, CENATHILA_CEN180_metrics_list, CENranLoc_CEN180_metrics_list, region_name, metric_name) {
  print(acc[acc_idx])

  fam_names <- sort(unique(CENATHILA_CEN180_metrics_list[[acc_idx]]$Family))

  # Analyse by ATHILA family
  permTestDF_fam <- data.frame()
  for(y in 1:length(fam_names)) {
 
    print(fam_names[y])

    features <- nrow(
                     CENATHILA_CEN180_metrics_list[[acc_idx]][
                       CENATHILA_CEN180_metrics_list[[acc_idx]]$Family == fam_names[y] , ]
                    ) / 2

    if(region_name %in% c("Upstream", "Downstream")) {

      CENATHILA <- mean(
                        CENATHILA_CEN180_metrics_list[[acc_idx]][
                          which(CENATHILA_CEN180_metrics_list[[acc_idx]]$Region == region_name &
                                CENATHILA_CEN180_metrics_list[[acc_idx]]$Family == fam_names[y]) ,
                          which(colnames(CENATHILA_CEN180_metrics_list[[acc_idx]]) == metric_name) ,
                          drop = T],
                        na.rm = T)

      CENranLoc_permsList <- CENranLoc_CEN180_metrics_list[[acc_idx]]
      CENranLoc <- foreach(x = iter(1:perms),
#                           .options.mpi = mpiopts,
                           .combine = "c",
                           .multicombine = T,
                           .maxcombine = perms+1e1,
                           .inorder = F) %dopar% {
        mean(
             CENranLoc_permsList[[x]][
               which(CENranLoc_permsList[[x]]$Region == region_name &
                     CENranLoc_permsList[[x]]$Family == fam_names[y]) ,
               which(colnames(CENranLoc_permsList[[x]]) == metric_name) ,
               drop = T],
             na.rm = T)
      }

    } else if(region_name == "Flanks") {

      CENATHILA <- mean(
                        CENATHILA_CEN180_metrics_list[[acc_idx]][
                          which(CENATHILA_CEN180_metrics_list[[acc_idx]]$Family == fam_names[y]) ,
                          which(colnames(CENATHILA_CEN180_metrics_list[[acc_idx]]) == metric_name) ,
                          drop = T],
                        na.rm = T)

      CENranLoc_permsList <- CENranLoc_CEN180_metrics_list[[acc_idx]]
      CENranLoc <- foreach(x = iter(1:perms),
#                           .options.mpi = mpiopts,
                           .combine = "c",
                           .multicombine = T,
                           .maxcombine = perms+1e1,
                           .inorder = F) %dopar% {
        mean(
             CENranLoc_permsList[[x]][
               which(CENranLoc_permsList[[x]]$Family == fam_names[y]) ,
               which(colnames(CENranLoc_permsList[[x]]) == metric_name) ,
               drop = T],
             na.rm = T)
      }

    }

    if(!is.na(CENATHILA)) {

      permDist <- CENranLoc
      observed <- CENATHILA
      expected <- mean(permDist, na.rm = T)

      # Calculate P-values and significance levels
      if(observed > expected) {
        permPval <- 1 - ( sum(observed > permDist, na.rm = T) / perms )
        MoreOrLessThanRandom <- "MoreThanRandom"
        alphaThreshold <- quantile(permDist, probs = 0.95, na.rm = T)[[1]]
      } else {
        permPval <- 1 - ( sum(observed < permDist, na.rm = T) / perms )
        MoreOrLessThanRandom <- "LessThanRandom"
        alphaThreshold <- quantile(permDist, probs = 0.05, na.rm = T)[[1]]
      }

      if(permPval == 0) { permPval <- minPval }

      permTestResults <- new("permTest",
                             alternative = MoreOrLessThanRandom,
                             alphaThreshold = alphaThreshold,
                             pval = permPval,
                             observed = observed,
                             expected = expected,
                             log2obsexp = log2( (observed+1) / (expected+1) ),
                             log2alpha  = log2( (alphaThreshold+1) / (expected+1) ),
                             permDist = permDist,
                             features = features,
                             fam = fam_names[y],
                             accession = acc[acc_idx],
                             metric = metric_name,
                             region = region_name)

      permTestDF_y <- data.frame(accession = permTestResults@accession,
                                 metric = permTestResults@metric,
                                 region = permTestResults@region,
                                 fam = permTestResults@fam,
                                 features = permTestResults@features,
                                 alternative = permTestResults@alternative,
                                 alphaThreshold = permTestResults@alphaThreshold,
                                 observed = permTestResults@observed,
                                 expected = permTestResults@expected,
                                 log2obsexp = permTestResults@log2obsexp,
                                 log2alpha = permTestResults@log2alpha,
                                 pval = permTestResults@pval)
      permTestDF_fam <- rbind(permTestDF_fam, permTestDF_y)

    }

  }


  # Analyse all ATHILA (not by ATHILA family)
  features <- nrow(
                   CENATHILA_CEN180_metrics_list[[acc_idx]]
                  ) / 2

  if(region_name %in% c("Upstream", "Downstream")) {

    CENATHILA <- mean(
                      CENATHILA_CEN180_metrics_list[[acc_idx]][
                        which(CENATHILA_CEN180_metrics_list[[acc_idx]]$Region == region_name) ,
                        which(colnames(CENATHILA_CEN180_metrics_list[[acc_idx]]) == metric_name) ,
                        drop = T],
                      na.rm = T)

    CENranLoc_permsList <- CENranLoc_CEN180_metrics_list[[acc_idx]]
    CENranLoc <- foreach(x = iter(1:perms),
#                         .options.mpi = mpiopts,
                         .combine = "c",
                         .multicombine = T,
                         .maxcombine = perms+1e1,
                         .inorder = F) %dopar% {
      mean(
           CENranLoc_permsList[[x]][
             which(CENranLoc_permsList[[x]]$Region == region_name) ,
             which(colnames(CENranLoc_permsList[[x]]) == metric_name) ,
             drop = T],
           na.rm = T)
    }

  } else if(region_name == "Flanks") {

    CENATHILA <- mean(
                      CENATHILA_CEN180_metrics_list[[acc_idx]][ ,
                        which(colnames(CENATHILA_CEN180_metrics_list[[acc_idx]]) == metric_name) ,
                        drop = T],
                      na.rm = T)

    CENranLoc_permsList <- CENranLoc_CEN180_metrics_list[[acc_idx]]
    CENranLoc <- foreach(x = iter(1:perms),
#                         .options.mpi = mpiopts,
                         .combine = "c",
                         .multicombine = T,
                         .maxcombine = perms+1e1,
                         .inorder = F) %dopar% {
      mean(
           CENranLoc_permsList[[x]][ ,
             which(colnames(CENranLoc_permsList[[x]]) == metric_name) ,
             drop = T],
           na.rm = T)
    }

  }

  if(!is.na(CENATHILA)) {
 
    permDist <- CENranLoc
    observed <- CENATHILA
    expected <- mean(permDist, na.rm = T)

    # Calculate P-values and significance levels
    if(observed > expected) {
      permPval <- 1 - ( sum(observed > permDist, na.rm = T) / perms )
      MoreOrLessThanRandom <- "MoreThanRandom"
      alphaThreshold <- quantile(permDist, probs = 0.95, na.rm = T)[[1]]
    } else {
      permPval <- 1 - ( sum(observed < permDist, na.rm = T) / perms )
      MoreOrLessThanRandom <- "LessThanRandom"
      alphaThreshold <- quantile(permDist, probs = 0.05, na.rm = T)[[1]]
    }

    if(permPval == 0) { permPval <- minPval }

    permTestResults <- new("permTest",
                           alternative = MoreOrLessThanRandom,
                           alphaThreshold = alphaThreshold,
                           pval = permPval,
                           observed = observed,
                           expected = expected,
                           log2obsexp = log2( (observed+1) / (expected+1) ),
                           log2alpha  = log2( (alphaThreshold+1) / (expected+1) ),
                           permDist = permDist,
                           features = features,
                           fam = "ATHILA",
                           accession = acc[acc_idx],
                           metric = metric_name,
                           region = region_name)

    permTestDF_all <- data.frame(accession = permTestResults@accession,
                                 metric = permTestResults@metric,
                                 region = permTestResults@region,
                                 fam = permTestResults@fam,
                                 features = permTestResults@features,
                                 alternative = permTestResults@alternative,
                                 alphaThreshold = permTestResults@alphaThreshold,
                                 observed = permTestResults@observed,
                                 expected = permTestResults@expected,
                                 log2obsexp = permTestResults@log2obsexp,
                                 log2alpha = permTestResults@log2alpha,
                                 pval = permTestResults@pval)

    permTestDF <- rbind(permTestDF_fam, permTestDF_all)

  } else {

    permTestDF <- permTestDF_fam

  }

  permTestDF[ with(permTestDF, order(fam)) , ]

}
