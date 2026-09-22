## Tests for the iBreastSubtypeR() launcher helper. .load_app_dependencies()
## replaces the former .attach_if(), whose requireNamespace() call used
## unsupported arguments and therefore loaded nothing.

test_that("the app launcher loads dependency namespaces and reports missing packages", {
  expect_error(
    BreastSubtypeR:::.load_app_dependencies(c("stringr", "no_such_package_for_ibreastsubtyper")),
    "Please install required package(s) before launching the app: no_such_package_for_ibreastsubtyper",
    fixed = TRUE
  )
  expect_null(BreastSubtypeR:::.load_app_dependencies("stringr"))
  expect_true(isNamespaceLoaded("stringr"))
  expect_null(BreastSubtypeR:::.load_app_dependencies(character(0)))
})
