context('qqqMS')

testthat('qqqMS-integrate',{

  example_files <-
    list.files(system.file(
      'extdata/test_data',
      pattern = '.mzML',
      full.names = TRUE
    ))

  example_pheno <-
    read.csv(
      system.file('extdata/test_data/pheno.csv'),
      head = TRUE,
      stringsAsFactors = FALSE
    )





}
)
