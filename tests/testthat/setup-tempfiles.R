#### Remove temporary files written by the tests ####
# Test helpers write small VCF, .ped, .map and .txt files with tempfile().
# Everything added to tempdir() while the tests run is deleted afterwards, so
# the test suite leaves no files behind.

.tempdir_before <- list.files(tempdir(), full.names = TRUE, all.files = TRUE, no.. = TRUE)

withr::defer({
  after <- list.files(tempdir(), full.names = TRUE, all.files = TRUE, no.. = TRUE)
  unlink(setdiff(after, .tempdir_before), recursive = TRUE)
}, testthat::teardown_env())
