context('qqqMS-integrate')

test_that('qqqMS-integrate',{

  example_files <-
    list.files(system.file(
      'extdata/test_data',package = 'qqqMS'),
      pattern = '.mzML',
      full.names = TRUE
    )

  example_pheno <-
    read.csv(
      system.file('extdata/test_data/phenoData.csv', package = 'qqqMS'),
      head = TRUE,
      stringsAsFactors = FALSE
    )

  pheno_a <- example_pheno
  names(pheno_a) <- c('sample', 'id')

  expect_error(integrate_peaks(example_files, pheno_a))
  expect_error(integrate_peaks(example_files, example_pheno[1:2,]))
  expect_error(integrate_peaks(example_files[1:2], example_pheno))

  qqq_ms_peaks <- integrate_peaks(example_files, example_pheno, pol = '1')

  expect_true(is.list(qqq_ms_peaks))
  expect_true(length(qqq_ms_peaks) == 2)

  expect_true(is.data.frame(qqq_ms_peaks[[1]]))
  expect_true(is.data.frame(qqq_ms_peaks[[2]]))

  expect_true(ncol(qqq_ms_peaks[['values']]) == length(example_files))

  }
)
