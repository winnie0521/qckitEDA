
dat_filt <- read.data("Data/counts.csv","Data/expt_design.csv")

test_that("read-in-data is right dimension", {
  expect_equal(ncol(dat_filt$wide), 13)
  expect_equal(nrow(dat_filt$wide), 14883)
  expect_equal(ncol(dat_filt$long), 8)
  expect_equal(nrow(dat_filt$long), 178596)
})
