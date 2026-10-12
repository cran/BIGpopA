test_that(".read_founders accepts NULL, vector, data.frame and file", {
  expect_equal(BIGpopA:::.read_founders(NULL), character(0))
  expect_equal(BIGpopA:::.read_founders(c("A", "B", "A")), c("A", "B"))
  expect_equal(BIGpopA:::.read_founders(data.frame(id = c("A", "B"))), c("A", "B"))
  f <- tempfile(fileext = ".txt")
  writeLines(c("A", "B"), f)
  expect_equal(BIGpopA:::.read_founders(f), c("A", "B"))
  expect_error(BIGpopA:::.read_founders(list(1)), "founders_file")
})
