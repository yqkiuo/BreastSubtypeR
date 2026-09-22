## Tests for the cohort preflight of the bundled Shiny app. get_methods()
## requires both ER and HER2 columns, so the preflight must not report "ready"
## when only one of them (or only TN) is present.

test_that("the cohort preflight requires both ER and HER2 columns", {
  skip_if_not_installed("shiny")
  env <- new.env(parent = asNamespace("shiny"))
  sys.source(system.file("ShinyBreastSubtypeR/server.R", package = "BreastSubtypeR"), env)

  both <- data.frame(ER = c("ER+", "ER-"), HER2 = c("HER2-", "HER2+"), stringsAsFactors = FALSE)
  res <- env$.summarize_cohort(both)
  expect_identical(res$kind, "ERHER2")
  expect_true(res$ok)

  er_only <- both[, "ER", drop = FALSE]
  res <- env$.summarize_cohort(er_only)
  expect_identical(res$kind, "ER")
  expect_false(res$ok)
  expect_match(res$msg, "missing: HER2", fixed = TRUE)

  her2_only <- both[, "HER2", drop = FALSE]
  res <- env$.summarize_cohort(her2_only)
  expect_identical(res$kind, "HER2")
  expect_false(res$ok)
  expect_match(res$msg, "missing: ER", fixed = TRUE)

  tn_only <- data.frame(TN = c("TN", "TN"), stringsAsFactors = FALSE)
  res <- env$.summarize_cohort(tn_only)
  expect_identical(res$kind, "TN")
  expect_false(res$ok)
  expect_match(res$msg, "missing: ER, HER2", fixed = TRUE)

  tn_full <- cbind(both, TN = c("nonTN", "TN"), stringsAsFactors = FALSE)
  res <- env$.summarize_cohort(tn_full)
  expect_identical(res$kind, "TN")
  expect_true(res$ok)

  invalid <- data.frame(ER = c("pos", "ER-"), HER2 = c("HER2-", "HER2+"), stringsAsFactors = FALSE)
  res <- env$.summarize_cohort(invalid)
  expect_false(res$ok)
  expect_match(res$msg, "Invalid ER: pos", fixed = TRUE)
})
