## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE
)

## -----------------------------------------------------------------------------
#  gert::git_branch_checkout("main")
#  gert::git_pull()

## -----------------------------------------------------------------------------
#  ver <- desc::desc_get_version()[1, 1:3]

## -----------------------------------------------------------------------------
#  gert::git_branch_create(paste0("pkgdown-v", ver), paste0("v", ver))

## -----------------------------------------------------------------------------
#  files <- c(
#    # overall site config
#    "_pkgdown.yml",
#    # the workflow that builds the site
#    ".github/workflows/pkgdown.yaml",
#    # readme and vignettes
#    "README.md", "vignettes",
#    # logo and favicon
#    "man/figures/", "pkgdown/",
#    # Author metadata and Config/Needs/Website
#    "DESCRIPTION"
#  )
#  
#  glue::glue("git checkout v{ver} -- {files}")

## -----------------------------------------------------------------------------
#  desc::desc_set_version(ver)

## -----------------------------------------------------------------------------
#  pkgdown::build_site()

## -----------------------------------------------------------------------------
#  usethis:::git_push_first()

