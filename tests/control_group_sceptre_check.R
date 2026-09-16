#!/usr/bin/env Rscript
# R-level check of the SCEPTRE driver's control-group resolution.
#
# bin/inference_sceptre.R needs sceptre, MuData and the Bioconductor stack to
# source, so this extracts the two self-contained functions and exercises them
# directly. tests/test_control_group_resolution.py runs this and skips when
# Rscript is unavailable.

args <- commandArgs(trailingOnly = TRUE)
script_path <- if (length(args) >= 1) args[1] else "bin/inference_sceptre.R"

src <- readLines(script_path, warn = FALSE)

# Take a top-level `name <- function(...) { ... }` definition by brace balance.
extract_function <- function(name) {
  start <- grep(sprintf("^%s <- function", name), src)
  if (length(start) != 1) {
    stop(sprintf("expected exactly one top-level definition of %s, found %d", name, length(start)))
  }
  depth <- 0L
  end <- start
  repeat {
    line <- src[end]
    depth <- depth + lengths(regmatches(line, gregexpr("\\{", line))) -
      lengths(regmatches(line, gregexpr("\\}", line)))
    if (depth == 0L && end > start) break
    end <- end + 1L
    if (end > length(src)) stop(sprintf("unbalanced braces in %s", name))
  }
  paste(src[start:end], collapse = "\n")
}

eval(parse(text = extract_function("resolve_sceptre_control_group")))
eval(parse(text = extract_function("write_control_group_metadata")))

failures <- character()
check <- function(label, expr) {
  ok <- tryCatch(isTRUE(expr), error = function(e) {
    failures <<- c(failures, sprintf("%s: error %s", label, conditionMessage(e)))
    FALSE
  })
  if (!ok) {
    failures <<- c(failures, label)
  }
}
check_error <- function(label, expr, pattern) {
  message_text <- tryCatch({
    force(expr)
    NA_character_
  }, error = function(e) conditionMessage(e))
  if (is.na(message_text)) {
    failures <<- c(failures, sprintf("%s: expected an error, got none", label))
  } else if (!grepl(pattern, message_text, fixed = TRUE)) {
    failures <<- c(failures, sprintf("%s: error did not mention '%s': %s", label, pattern, message_text))
  }
}

# --- nothing requested: decide from the analysis ---------------------------
r <- resolve_sceptre_control_group(NULL, TRUE, 500L)
check("auto low-MOI with NT cells -> nt_cells", r$control_group == "nt_cells")
check("auto low-MOI provenance", r$provenance == "auto")

r <- resolve_sceptre_control_group(NULL, FALSE, 500L)
check("auto high-MOI -> complement", r$control_group == "complement")

r <- resolve_sceptre_control_group(NULL, TRUE, 0L)
check("auto low-MOI without NT cells -> complement", r$control_group == "complement")

r <- resolve_sceptre_control_group("auto", TRUE, 500L)
check("the literal 'auto' means unrequested", r$provenance == "auto")
r <- resolve_sceptre_control_group("", TRUE, 500L)
check("the empty string means unrequested", r$provenance == "auto")

# --- requested, and honoured ----------------------------------------------
r <- resolve_sceptre_control_group("nt_cells", TRUE, 500L)
check("nt_cells honoured at low MOI", r$control_group == "nt_cells")
check("nt_cells provenance", r$provenance == "requested")

r <- resolve_sceptre_control_group("complement", TRUE, 500L)
check("complement honoured at low MOI, not overridden to nt_cells", r$control_group == "complement")
check("complement provenance", r$provenance == "requested")

r <- resolve_sceptre_control_group("complement", FALSE, 0L)
check("complement honoured at high MOI", r$control_group == "complement")

r <- resolve_sceptre_control_group(" NT_Cells ", TRUE, 500L)
check("case and whitespace tolerated", r$control_group == "nt_cells")

# --- requested, and refused ------------------------------------------------
check_error(
  "nt_cells at high MOI stops",
  resolve_sceptre_control_group("nt_cells", FALSE, 500L),
  "this is a high-MOI analysis"
)
check_error(
  "nt_cells without a control population stops",
  resolve_sceptre_control_group("nt_cells", TRUE, 0L),
  "no control population to compare against"
)
check_error(
  "an unknown control group stops",
  resolve_sceptre_control_group("control-anchored", TRUE, 500L),
  "is not a SCEPTRE control group"
)

# --- the recorded metadata -------------------------------------------------
tmp <- tempfile(fileext = ".json")
write_control_group_metadata("nt_cells", "nt_cells", "requested", TRUE, 500L, path = tmp)
recorded <- paste(readLines(tmp), collapse = " ")
check('metadata names the control group', grepl('"control_group": "nt_cells"', recorded, fixed = TRUE))
check('metadata names the provenance', grepl('"provenance": "requested"', recorded, fixed = TRUE))
check('metadata names the analysis MOI', grepl('"analysis_moi": "low"', recorded, fixed = TRUE))
check('metadata names the NT cell count', grepl('"non_targeting_cells": 500', recorded, fixed = TRUE))

write_control_group_metadata("complement", NA_character_, "auto", FALSE, NA_integer_, path = tmp)
recorded <- paste(readLines(tmp), collapse = " ")
check('an unrequested group records null', grepl('"control_group_requested": null', recorded, fixed = TRUE))
check('an unknown count records null', grepl('"non_targeting_cells": null', recorded, fixed = TRUE))

if (length(failures)) {
  cat("FAILED:\n")
  cat(paste0(" - ", failures, collapse = "\n"), "\n")
  quit(status = 1)
}
cat("control_group_sceptre_check: all checks passed\n")
