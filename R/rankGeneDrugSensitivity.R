
#' @importFrom stats complete.cases
#' @importFrom stats p.adjust

##TODO:: Convert this to roxygen2 format
#################################################
## Rank genes based on drug effect in the Connectivity Map
##
## inputs:
##      - data: gene expression data matrix
##            - drugpheno: sensititivity values fo thr drug of interest
##            - type: cell or tissue type for each experiment
##            - duration: experiment duration in hours
##      - batch: experiment batches
##            - single.type: Should the statitsics be computed for each cell/tissue type separately?
##      - nthread: number of parallel threads (bound to the maximum number of cores available)
##
## outputs:
## list of datafraes with the statistics for each gene, for each type
##
## Notes:    duration is not taken into account as only 4 perturbations lasted 12h, the other 6096 lasted 6h
#################################################

rankGeneDrugSensitivity <- function (data, drugpheno, type, batch,
                                     single.type=FALSE, standardize = "SD",
                                     nthread=1, verbose=FALSE) {
  if (nthread != 1) {
    availcore <- parallel::detectCores()
    if (missing(nthread) || nthread < 1 || nthread > availcore) {
      nthread <- availcore
    }
  }
  # Set multicore options
  op <- options()
  options(mc.cores=nthread)
  on.exit(options(op))

  if(is.null(dim(drugpheno))){

    drugpheno <- data.frame(drugpheno)

  } else if(!is(drugpheno, "data.frame")) {
    drugpheno <- as.data.frame(drugpheno)

  }

  if (missing(type) || all(is.na(type))) {
    type <- array("other", dim=nrow(data), dimnames=list(rownames(data)))
  }
  if (missing(batch) || all(is.na(batch))) {
    batch <- array(1, dim=nrow(data), dimnames=list(rownames(data)))
  }
  if (any(c(nrow(drugpheno), length(type), length(batch)) != nrow(data))) {
    stop("length of drugpheno, type, duration, and batch should be equal to the number of rows of data!")
  }
  rownames(drugpheno) <- names(type) <- names(batch) <- rownames(data)

  res <- NULL
  utype <- sort(unique(as.character(type)))
  ltype <- list("all"=utype)
  if (single.type) {
    ltype <- c(ltype, as.list(utype))
    names(ltype)[-1] <- utype
  }
  res <- NULL
  ccix <- complete.cases(data, type, batch, drugpheno)
  nn <- sum(ccix)
  if(!any(unlist(lapply(drugpheno,is.factor)))){
     if(ncol(drugpheno)>1){
      ##FIXME:: NAMES!
      nc <- lapply(seq_len(ncol(drugpheno)), function(i){

        est <- paste("estimate", i, sep=".")
        se <-  paste("se", i, sep=".")
        tstat <- paste("tstat", i, sep=".")

        nc <- c(est, se, tstat)
        return(nc)

      })
      nc  <- c(nc, n=nn, "fstat"=NA, "pvalue"=NA, "fdr")
    } else {
      nc  <- c("estimate", "se", "n", "tstat", "fstat", "pvalue", "df", "fdr")
    }
  } else {
    nc  <- c("estimate", "se", "n", "pvalue", "fdr")
  }



  for (ll in seq_along(ltype)) {
    iix <- !is.na(type) & is.element(type, ltype[[ll]])

    data.not.all.na <- apply(data[iix,,drop=FALSE], 1, function(x) {
      any(!is.na(x))
    })
    drugpheno.not.all.na <- apply(drugpheno[iix,,drop=FALSE], 1, function(x) {
      any(!is.na(x))
    })
    type.not.all.na <- vapply(type[iix], function(x) {
      !is.na(x)
    }, FUN.VALUE=logical(1))
    batch.not.all.na <- vapply(batch[iix], function(x) {
      !is.na(x)
    }, FUN.VALUE=logical(1))

    ccix <- data.not.all.na & drugpheno.not.all.na & type.not.all.na & batch.not.all.na



    if (sum(ccix) < 3) {
      ## not enough experiments
      rest <- list(matrix(NA, nrow=ncol(data), ncol=length(nc), dimnames=list(colnames(data), nc)))
      res <- c(res, rest)
    } else {
      splitix <- parallel::splitIndices(nx=ncol(data), ncl=nthread)
      ##TODO:: Can we reimpement this without length?
      splitix <- splitix[vapply(splitix, length, FUN.VALUE=numeric(1)) > 0]
      mcres <- BiocParallel::bplapply(splitix, function(x, data, type, batch, drugpheno, standardize) {
        res <- t(apply(data[ , x, drop=FALSE], 2, geneDrugSensitivity, type=type, batch=batch, drugpheno=drugpheno, verbose=verbose, standardize=standardize))
        return(res)
      }, data=data[iix, , drop=FALSE], type=type[iix], batch=batch[iix], drugpheno=drugpheno[iix,,drop=FALSE], standardize=standardize)
      rest <- do.call(rbind, mcres)
      rest <- cbind(rest, "fdr"=p.adjust(rest[ , "pvalue"], method="fdr"))
      res <- c(res, list(rest))
    }
  }
  names(res) <- names(ltype)
  return(res)
}