dat_filt <- read.data("Data/counts.csv","Data/expt_design.csv")

annoDat <- annot_data(dat_filt)

test_that("read-in-data is right dimension", {
  expect_equal(ncol(annoDat), 9)
  expect_equal(nrow(annoDat), 14856)

})
