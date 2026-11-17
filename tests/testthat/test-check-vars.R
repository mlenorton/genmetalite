test_that("check_vars flags extra variables", {
  ref <- data.frame(a = 1:3, b = 1:3)
  coh <- data.frame(a = 1:3, b = 1:3, c = 1:3)
  expect_warning(check_vars(ref, coh), "not the first")
  expect_error(check_vars_stop(ref, coh), "not the first")
})

test_that("check_vars flags factor levels missing from the template", {
  ref <- data.frame(g = factor(c("a", "b")))
  coh <- data.frame(g = factor(c("a", "c")))
  expect_warning(check_vars(ref, coh), "level c present")
})

test_that("ordered factors are detected", {
  ## class(x) is c("ordered", "factor") for an ordered factor, so the original
  ## `which(col_classes == "factor")` test silently skipped these columns.
  ref <- data.frame(g = factor(c("a", "b"), ordered = TRUE))
  coh <- data.frame(g = factor(c("a", "c"), ordered = TRUE))
  expect_warning(check_vars(ref, coh), "level c present")
})

test_that("check_vars is silent when everything matches", {
  ref <- data.frame(g = factor(c("a", "b", "c")), x = 1:3)
  coh <- data.frame(g = factor(c("a", "b")), x = 1:2)
  expect_silent(out <- check_vars(ref, coh))
  expect_length(out, 0L)
})

test_that("multiple missing levels are reported together", {
  ref <- data.frame(g = factor("a"))
  coh <- data.frame(g = factor(c("b", "c")))
  expect_warning(check_vars(ref, coh), "levels b, c")
})

test_that("action = 'message' reports without warning", {
  ref <- data.frame(g = factor("a"))
  coh <- data.frame(g = factor(c("a", "b")))
  expect_message(check_vars(ref, coh, action = "message"), "level b")
})
