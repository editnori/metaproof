# Installs the R packages the reference engine was checked with into reference_r/lib.
#
#   Rscript reference_r/install_packages.R            # into reference_r/lib
#   Rscript reference_r/install_packages.R <folder>   # into another folder
#
# The packages come from a dated Posit Package Manager snapshot of CRAN, so the same
# versions install every time. They build from source, which needs a C compiler
# (xcode-select --install on a Mac). After installing, the script checks every version
# against reference_r/R_PACKAGES_LOCK.txt and stops with an error if one differs.

SNAPSHOT <- "https://packagemanager.posit.co/cran/2026-11-01"
WANT <- c("metafor", "metadat", "jsonlite", "yaml")

file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
here <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
args <- commandArgs(TRUE)
lib <- if (length(args)) args[1] else file.path(here, "lib")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
lib <- normalizePath(lib)

lock_lines <- readLines(file.path(here, "R_PACKAGES_LOCK.txt"))
lock <- regmatches(lock_lines, regexec("^([A-Za-z][A-Za-z0-9.]*) +([0-9][0-9.-]*) +", lock_lines))
lock <- do.call(rbind, lapply(Filter(length, lock), function(m) m[2:3]))
lock <- setNames(lock[, 2], lock[, 1])
lock <- lock[names(lock) != "R"]
stopifnot(all(WANT %in% names(lock)))

r_lock <- sub("^R ([0-9.]+).*", "\\1", grep("^R [0-9]", lock_lines, value = TRUE)[1])
if (!identical(as.character(getRversion()), r_lock))
  message("Note: this is R ", getRversion(), "; the kit was checked with R ", r_lock, ".")

# Only lib and R's own base library are visible, so a dependency already sitting in a
# personal library still gets installed here (the engine runs with R_LIBS_USER=lib).
.libPaths(lib, include.site = FALSE)
install.packages(WANT, lib = lib, repos = SNAPSHOT, type = "source",
                 dependencies = c("Depends", "Imports", "LinkingTo"))

got <- installed.packages(lib.loc = lib)[, "Version"]
bad <- character(0)
for (p in names(lock)) {
  have <- if (p %in% names(got)) got[[p]] else "missing"
  ok <- identical(have, lock[[p]])
  cat(sprintf("%-4s %-9s want %-11s have %s\n", if (ok) "OK" else "BAD", p, lock[[p]], have))
  if (!ok) bad <- c(bad, p)
}
for (p in WANT) suppressPackageStartupMessages(library(p, lib.loc = lib, character.only = TRUE))
if (length(bad)) stop("versions differ from R_PACKAGES_LOCK.txt: ", paste(bad, collapse = ", "))
cat("All", length(lock), "packages match R_PACKAGES_LOCK.txt and load from", lib, "\n")
