#' Detect and integrate MRM-MS features
#'
#'
#' @param files a character vector of \code{.mzML} files which have been convereted from MRM-MS \code{.mzML} to LC-MS style \code{.mzML} files using \code{MRMConverteR}
#' @param phenoData a \code{data.frame}. \code{phenoData} must contain the following columns;
#' \itemize{
#'     \item{\strong{fileName}} file name
#'     \item{\strong{name}} sample name
#' }
#'
#' Providing these two columns are present; \code{phenoData} can contain as many other columns as the user wishes.
#'
#' @return a list of two elements
#' \itemize{
#'     \item{values} values for the integrated peak areas of all features detected
#'     \item{definitions} a \code{data.frame} of feature defintions. See \link[xcms]{featureDefinitions} for more details
#' }
#' @export
#' @importFrom methods new
#' @importFrom xcms readMSData MatchedFilterParam findChromPeaks PeakDensityParam groupChromPeaks featureValues featureDefinitions

integrate_peaks <- function(files, phenoData)
  {

  if(length(files) != nrow(phenoData)) {
    stop(
      deparse(substitute(files)),
      ' and ',
      deparse(substitute(phenoData)),
      ' do not have the same dimensions',
      call. = FALSE
    )
  }

  if (!'fileName' %in% names(phenoData)) {
    stop('no `fileName` column found in ', deparse(substitute(phenoData)), call. = FALSE)
  }

  if (!'name' %in% names(phenoData)) {
    stop('no `name` column found in ', deparse(substitute(phenoData)), call. = FALSE)
  }

  phenoData <- data.frame(phenoData, group = rep(1, nrow(phenoData)))

  pheno_ob <- new('NAnnotatedDataFrame', phenoData)

  xcraw <- readMSData(files, pdata = pheno_ob, mode = 'onDisk')

  matched_filt_params <-
    MatchedFilterParam(
      fwhm = 30,
      snthresh = 1.0,
      binSize = 0.01,
      steps = 1.0,
      mzdif = -2.0
    )

  xcpeaks <- findChromPeaks(xcraw, matched_filt_params)

  grp_parmas <- PeakDensityParam(
    sampleGroups = xcpeaks$group,
    minFraction = 0.1,
    bw = 30,
    binSize = 0.01
  )

  xcgrp <- groupChromPeaks(xcpeaks, param = grp_parmas)

  feature_values <- data.frame(featureValues(xcgrp, value = 'into'))

  feature_def <- data.frame(featureDefinitions(xcgrp))

  feature_def[, 'mzmed'] <- round(feature_def[, 'mzmed'], digits = 3)
  feature_def[, 'rt'] <- round(feature_def[, 'rtmed'], digits = 1)

  feature_def[, 'mzmax'] <- round(feature_def[, 'mzmax'], digits = 3)
  feature_def[, 'mzmin'] <- round(feature_def[, 'mzmin'], digits = 3)

  return(list(values = feature_values), definitions = feature_def)

  }