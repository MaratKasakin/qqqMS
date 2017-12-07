#' Extract Targets
#'
#' Use a user-defined list of MRM-MR targets to create a filtered feature values matrix
#'
#' @param features a list which has been produced using the \link{integrate_peaks} function
#' @param targets a \code{data.frame} of MRM-MS targets to extract the the xcms processed object. \code{targets} must contain the following columns;
#' \itemize{
#'     \item{\strong{name}} the generic name of the target
#'     \item{\strong{mz}} the m/z value of the target (as it appears in the raw file)
#'     \item{\strong{rt}} the approximate retention time of the target
#' }
#'
#' @return a \code{data.frame} of feature values for only those which corresponds to an identified target
#'
#' @export


extract_targets <- function(features, targets)
{
  if (!is.data.frame(targets)) {
    stop(deparse(substitute(targets)), ' is not a `data.frame`', call. = FALSE)
  }

  if (!'name' %in% names(targets)) {
    stop('no `name` column found in ', deparse(substitute(targets)), call. = FALSE)
  }

  if (!'mz' %in% names(targets)) {
    stop('no `mz` column found in ', deparse(substitute(targets)), call. = FALSE)
  }

  if (!'rt' %in% names(targets)) {
    stop('no `rt` column found in ', deparse(substitute(targets)), call. = FALSE)
  }


  if (!is.list(features)) {
    stop(deparse(substitute(features)), ' is not a `list`', call. = FALSE)
  }

  if (length(features) != 2) {
    stop(deparse(substitute(features)), ' is not a list of two elements', call. = FALSE)
  }

  deftmp <- features[['definitions']]
  valtmp <- features[['values']]

  extid <- NULL
  for (i in 1:nrow(targets)) {
    mz <- targets[i, 'mz']
    rt <- targets[i, 'rt']

    idxf <-
      which(deftmp[, 'mzmed'] == mz &
              deftmp[, 'rtmax'] >= (rt - 2.0) &
              deftmp[, 'rtmin'] <= (rt + 2.0))

    if (length(idxf) == 0) {
      idxf <- which(deftmp[, 'mzmax'] >= mz & deftmp[, 'mzmin'] <= mz &
                      deftmp[, 'rtmax'] >= (rt - 1.5) &
                      deftmp[, 'rtmin'] <= (rt + 1.5))
    }

    extid[[i]] <- idxf
  }


  targets_tmp <- data.frame(targets, id = rownames(deftmp)[extid])

  targvals_idx <- match(targets_tmp[, 'id'], rownames(valtmp))

  target_values <- valtmp[targvals_idx, ]
  rownames(target_values) <- targets_tmp[, 'name']

  tval <- t(target_values)

  tval_out <-
    data.frame(
      name = rownames(tval),
      tval,
      row.names = NULL,
      check.names = FALSE
    )

  return(tval_out)
}
