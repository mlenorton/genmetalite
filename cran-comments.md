## Test environments

* local Ubuntu 24.04, R 4.3.3

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Notes

* `ebal` and `cobalt` are in Suggests and are used only by `ebal_weights()`,
  which checks for them with `requireNamespace()` and errors with an
  installation hint if absent. Tests that need them use
  `skip_if_not_installed()`.
