# Contributing to fastDiversity

Thanks for your interest in fastDiversity. Bug reports, questions, and
contributions are all welcome.

## Reporting bugs

Please open an issue at
<https://github.com/eilishmcmaster/fastDiversity/issues> and include:

- a short description of the problem and what you expected to happen
- a minimal reproducible example, ideally using the bundled data
  (`data(example_data)`) or a small simulated genotype matrix
- the output of
  [`sessionInfo()`](https://rdrr.io/r/utils/sessionInfo.html)

## Asking for help

If you have a question about how to use the package, or about how a
statistic is calculated, please open an issue. Check the [documentation
site](https://eilishmcmaster.github.io/fastDiversity/) and existing
issues first, as your question may already be answered.

## Contributing code or documentation

1.  Fork the repository and create a branch for your change.

2.  Make your changes. If you change a function’s documentation, edit
    the roxygen comments in `R/` (not the files in `man/`) and run
    `devtools::document()`.

3.  Add or update tests in `tests/testthat/` for any change in
    behaviour.

4.  Check that the package still passes:

    ``` r

    devtools::test()   # run the test suite
    devtools::check()  # full R CMD check
    ```

5.  Add a bullet describing the change to `NEWS.md`.

6.  Open a pull request describing what you changed and why. The R CMD
    check GitHub Action runs automatically on every pull request.

For larger changes, please open an issue to discuss the idea before
starting work.
