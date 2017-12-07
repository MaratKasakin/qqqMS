#' Detect and integrate MRM-MS features
#'
#' Use the \code{matchedFilter} algorithim \code{xcms} for feature detection in MRM-MS which have been converted to LC-MS style \code{.mzML} files.
#'
#' @param files a character vector of \code{.mzML} files which have been convereted from MRM-MS \code{.mzML} to LC-MS style \code{.mzML} files using \code{MRMConverteR}
#' @param phenoData a \code{data.frame}. \code{phenoData} must contain the following columns;
#' \itemize{
#'     \item{\strong{fileName}} file name
#'     \item{\strong{name}} sample name
#' }
#'
#' Providing these two columns are present; \code{phenoData} can contain as many other columns as the user wishes.
#' @param pol a charatcer string of either \code{1} for positive mode or \code{-1} for negative mode. If input data only contains one polarity mode, then \code{polarity} can be \code{NULL}
#' @return a list of two elements
#' \itemize{
#'     \item{values} values for the integrated peak areas of all features detected
#'     \item{definitions} a \code{data.frame} of feature defintions. See \code{\link[xcms]{XCMSnExp-class}} for more details
#' }
#' @export
#' @importFrom methods new
#' @importFrom xcms MatchedFilterParam findChromPeaks PeakDensityParam groupChromPeaks featureValues featureDefinitions ObiwarpParam adjustRtime

#'
integrate_peaks <- function(files, phenoData, pol)
{
  if (length(files) != nrow(phenoData)) {
    stop(
      deparse(substitute(files)),
      ' and ',
      deparse(substitute(phenoData)),
      ' do not have the same dimensions',
      call. = FALSE
    )
  }

  if (!'fileName' %in% names(phenoData)) {
    stop('no `fileName` column found in ',
         deparse(substitute(phenoData)),
         call. = FALSE)
  }

  if (!'name' %in% names(phenoData)) {
    stop('no `name` column found in ', deparse(substitute(phenoData)), call. = FALSE)
  }

  phenoData <-
    data.frame(phenoData, sample_group = rep(1, nrow(phenoData)))

  pheno_ob <- new('NAnnotatedDataFrame', phenoData)

  xcraw <-
    MSnbase::readMSData(files, pdata = pheno_ob, mode = 'onDisk')

  matched_filt_params <-
    MatchedFilterParam(
      fwhm = 40,
      snthresh = 1,
      binSize = 0.01,
      steps = 1,
      mzdiff = -2,
      max = 5
    )

  # polarity check

  plen <- length(unique(xcraw@featureData@data$polarity))

  if (plen == 2) {
    xcraw <- split(xcraw, xcraw@featureData@data$polarity)

    if (pol == '-1') {
      xcraw <- xcraw$`-1`
    }

    if (pol == '1') {
      xcraw <- xcraw$`1`
    }
  }

  xcpeaks <- findChromPeaks(xcraw, matched_filt_params)

  xcpeaks <- adjustRtime(xcpeaks, param = ObiwarpParam(binSize = 0.01))

  grp_parmas <- PeakDensityParam(
    sampleGroups = xcpeaks$sample_group,
    minFraction = 0.1,
    bw = 60,
    binSize = 0.05
  )


  xcgrp <- groupChromPeaks(xcpeaks, param = grp_parmas)

  feature_values <-
    data.frame(featureValues(xcgrp, value = 'into', intensity = 'into'))

  feature_def <- data.frame(featureDefinitions(xcgrp))

  feature_def[, 'mzmed'] <-
    round(feature_def[, 'mzmed'], digits = 3)
  feature_def[, 'rt'] <- round(feature_def[, 'rtmed'], digits = 1)

  feature_def[, 'mzmax'] <-
    round(feature_def[, 'mzmax'], digits = 3)
  feature_def[, 'mzmin'] <-
    round(feature_def[, 'mzmin'], digits = 3)

  return(list(values = feature_values, definitions = feature_def))

}
