## Build the nhanes_bmi dataset shipped with svyrcs.
##
## Source files live outside the package, in the analysis repository this package grew out of:
##   ../data/derived/analytic_dataset.rds        merged NHANES demographics / exam / lab files
##   ../data/mortality/nhanes_mortality.sas7bdat NCHS 2019 Linked Mortality File
##
## Both are US Government public domain, so a derived extract may be redistributed.
##
## Two survey cycles (2005-2006 and 2007-2008) are kept in full. Taking whole cycles keeps every
## sampling stratum and PSU intact, which is what makes the shipped design valid; sampling rows or
## PSUs would not. Run from the package root:
##
##   Rscript data-raw/make_nhanes_bmi.R

stopifnot(requireNamespace("haven", quietly = TRUE))

cycles <- c("2005-2006", "2007-2008")
n_cycles <- length(cycles)

demo <- as.data.frame(readRDS(file.path("..", "data", "derived", "analytic_dataset.rds")))
mort <- haven::read_sas(file.path("..", "data", "mortality", "nhanes_mortality.sas7bdat"))
mort <- data.frame(
  SEQN = mort$seqn,
  eligstat = mort$eligstat,
  mortstat = mort$mortstat,
  permth = mort$permth_exm
)

d <- merge(demo, mort, by = "SEQN")
d <- d[d$cycle_label %in% cycles, ]

keep <- with(d,
  examined_adult_nonpregnant == 1 &
    eligstat == 1 &
    !is.na(bmi) & !is.na(age) & !is.na(sex) & !is.na(race) &
    !is.na(permth) & !is.na(mortstat) & permth > 0
)
d <- d[keep, ]

nhanes_bmi <- data.frame(
  seqn       = as.integer(d$SEQN),
  cycle      = factor(d$cycle_label, levels = cycles),
  psu        = as.integer(d$SDMVPSU),
  strata     = as.integer(d$SDMVSTRA),
  ## Two 2-year cycles combined: the MEC weight is divided by the number of cycles pooled,
  ## per the NHANES analytic guidelines.
  weight     = d$raw_mec_weight / n_cycles,
  age        = as.integer(d$age),
  sex        = d$sex,
  race       = d$race,
  bmi        = round(d$bmi, 2),
  tchol      = as.integer(d$tchol),
  hdl        = as.integer(d$hdl),
  glucose    = round(d$glucose, 1),
  statin_use = as.integer(d$statin_use),
  ## Follow-up from the MEC exam to death or 31 December 2019, in years.
  time       = round(d$permth / 12, 3),
  event      = as.integer(d$mortstat),
  stringsAsFactors = FALSE
)
nhanes_bmi$high_chol <- as.integer(nhanes_bmi$tchol >= 240)
rownames(nhanes_bmi) <- NULL
nhanes_bmi <- nhanes_bmi[order(nhanes_bmi$strata, nhanes_bmi$psu, nhanes_bmi$seqn), ]
rownames(nhanes_bmi) <- NULL

dir.create("data", showWarnings = FALSE)
save(nhanes_bmi, file = file.path("data", "nhanes_bmi.rda"), compress = "xz", version = 3)

cat("rows:", nrow(nhanes_bmi),
    " strata:", length(unique(nhanes_bmi$strata)),
    " psu-in-strata:", nrow(unique(nhanes_bmi[c("strata", "psu")])),
    " events:", sum(nhanes_bmi$event), "\n")
cat("file size:", format(file.size(file.path("data", "nhanes_bmi.rda")), big.mark = ","), "bytes\n")
