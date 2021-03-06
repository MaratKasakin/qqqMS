context('qqqMS-extract')

test_that('qqqMS-extract',{

  example_files <-
    list.files(system.file(
      'extdata/test_data',package = 'qqqMS'),
      pattern = '.mzML',
      full.names = TRUE)

  example_pheno <-
    read.csv(
      system.file('extdata/test_data/phenoData.csv', package = 'qqqMS'),
      head = TRUE,
      stringsAsFactors = FALSE
    )

  example_targets <-
    read.csv(
      system.file('extdata/test_data/targets.csv', package = 'qqqMS'),
      head = TRUE,
      stringsAsFactors = FALSE
    )

  example_targets <- example_targets[which(example_targets$polarity == '1'),]

  int_res <- integrate_peaks(example_files, example_pheno, pol = '1')

  targ_res <- extract_targets(int_res, example_targets)

  expect_true(ncol(targ_res) == (nrow(example_targets) + 1))

  expect_true(is.data.frame(targ_res))
  expect_true(names(targ_res)[1] == 'name')

  expect_true(all(example_targets[,'name'] %in% names(targ_res)))

  }
)
