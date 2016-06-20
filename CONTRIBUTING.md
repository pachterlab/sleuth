# contributing to `sleuth`

Firstly -- thank you for being interested in helping us!
What follows is a short document on how to contribute to `sleuth`.

## development cycle

The basic development cycle look something like this:

1. Make some modifications
2. Make sure all of the tests work: `devtools::test()`
3. Under most circumstances you should write your own tests

## style

Style is checked using [lintr](https://github.com/jimhester/lintr).
The linters that are set can be found in the configuration file `.lintr`.
Style can be checked for the entire project using `lintr::lint_package()`.

**All proposed code must pass the lint tests** otherwise it will not be accepted into the main branch.

## types of pull requests

Please note that this project is still in alpha stages and the statistics are constantly being developed while the paper is being written.
Because of this, we are currently not accepting changes to the statistics.
If you are interested in adding features, please check with us first by opening a Github issue.
Thanks!
