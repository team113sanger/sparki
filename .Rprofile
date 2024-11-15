# Typically we want to activate the project's renv environment when we start an
# R but if the renv directory is not present, we want to create it and
# initialize a new renv environment. When working with Docker & bind mounts, the
# renv directory may be masked by the bind mount, so we need to ensure that the
# renv directory is present and initialized when we start an R session.
#
# Due to the configurations of the RENV_PATHS_CACHE and RENV_PATHS_LIBRARY
# libaries do not need to be reinstalled when the renv directory is recreated.
if (!file.exists("renv/activate.R")) {
  renv::init(bare = TRUE)
} else {
  source("renv/activate.R")
}

# During development, you can use the following code to load the devtools and
# usethis packages without printing the package loading messages to the console:
if (interactive()) {
  suppressMessages(require(devtools))
  suppressMessages(require(usethis))
}
