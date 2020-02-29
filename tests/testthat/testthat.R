context("p-value of Q-tilde for strong linear trend")
library(testcorr)

qtilde.pv <- ac.test(1000000 * seq(1, 100, 1), max.lag = 10, plot = FALSE, table = FALSE)$pvqtilde

test_that("p-value of Q-tilde for strong linear trend is 0", {
  expect_equal(qtilde.pv[10], 0)
})
