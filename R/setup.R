# Sourced by every script. Run scripts from the repository root.

if (!file.exists("config.yml")) stop("Run scripts from the repository root")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(data.table)
  library(vroom)
})

Sys.setenv(VROOM_SHOW_PROGRESS = "false")

`%||%` <- function(x, y) if (is.null(x)) y else x

cfg <- yaml::read_yaml("config.yml")
if (file.exists("config.local.yml")) cfg <- modifyList(cfg, yaml::read_yaml("config.local.yml"))
traits <- names(cfg$traits)

for (f in list.files("R", pattern = "^fn_.*\\.R$", full.names = TRUE)) source(f)

res_dir <- function(...) {
  d <- file.path(cfg$results_dir, ...)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

set.seed(20250101)
