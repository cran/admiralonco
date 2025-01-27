## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(admiraldev)

## ----echo=FALSE---------------------------------------------------------------
list_cases <- tibble::tribble(
  ~"Use case",
  ~"Use Recommended by GCIG",
  ~"Not Standard and Needs Further Validation",
  ~"Not Recommended by GCIG",
  "First-line trials", "CA-125 progression", "", "CA-125 response",
  "Maintenance or consolidation trials", "", "CA-125 response and progression ", "",
  "Relapse trials", "CA-125 response and progression", "", ""
)

library(magrittr)

list_cases %>%
  gt::gt() %>%
  gt::cols_label_with(fn = ~ gt::md(paste0("**", .x, "**"))) %>%
  gt::tab_header(
    title = "GCIG recommendations for CA-125 criteria for response and progression in various clinical situations"
  )

## ----echo=FALSE---------------------------------------------------------------
# styler: off
list_resp <- tibble::tribble(
  ~"CA-125 response per GCIG",            ~"CA-125 response mapped",
  "Response within Normal Range",         "CR",
  "Response but not within Normal Range", "PR",
  "Non-Response/Non-PD",                  "SD",
  "PD",                                   "PD",
  "NE",                                   "NE"
)

knitr::kable(list_resp)
# styler: on

## ----echo=FALSE---------------------------------------------------------------
# styler: off
list_supp <- tibble::tribble(
  ~"QNAM", ~"QLABEL", ~"QVAL", ~"Purpose", ~"Use case",
  "`CA125EFL`", "CA-125 response evaluable", "Y/N", "Indicates population evaluable for CA-125 response (baseline CA-125 >= 2 * ULRR and no mouse antibodies)", "`CA125EFL` variable",
  "`CAELEPRE`", "Elevated pre-treatment CA-125", "Y/N", "Indicates CA-125 level at baseline (Y - elevated, N - not elevated)", "Derivation of PD category (`MCRIT1`/`MCRIT1ML`/`MCRIT1MN`)",
  "`MOUSEANT`", "Received mouse antibodies", "Y", "Indicates if a prohibited therapy was received", "Derivation of `ANL02FL`",
  "`CA50RED`", ">=50% reduction from baseline", "Y", "Indicates response, but does not distinguish between CR and PR", "Not used in further derivations",
  "`CANORM2X`", "CA125 normal, lab increased >=2x ULRR", "Y", "Indicates PD category A or C", "Derivation of PD category (`MCRIT1`/`MCRIT1ML`/`MCRIT1MN`)",
  "`CNOTNORM`", "CA125 not norm, lab increased >=2x nadir", "Y", "Indicates PD category B", "Derivation of PD category (`MCRIT1`/`MCRIT1ML`/`MCRIT1MN`)"
)

knitr::kable(list_supp)
# styler: on

## ----warning=FALSE, message=FALSE---------------------------------------------
library(admiral)
library(admiralonco)
library(pharmaversesdtm)
library(pharmaverseadam)
library(metatools)
library(dplyr)
library(tibble)

## ----message=FALSE------------------------------------------------------------
adsl <- pharmaverseadam::adsl
# GCIG sdtm data
rs <- pharmaversesdtm::rs_onco_ca125
supprs <- pharmaversesdtm::supprs_onco_ca125

rs <- combine_supp(rs, supprs)
rs <- convert_blanks_to_na(rs)

## ----eval=TRUE, echo=FALSE----------------------------------------------------
dataset_vignette(
  rs,
  display_vars = exprs(USUBJID, RSTESTCD, RSCAT, RSSTRESC, VISIT, CA125EFL, CAELEPRE, CA50RED, CANORM2X, CNOTNORM, MOUSEANT)
)

## ----echo=FALSE---------------------------------------------------------------
# select subjects from adsl such that there is one subject without RS data
rs_subjects <- unique(rs$USUBJID)
adsl_subjects <- unique(adsl$USUBJID)
adsl <- filter(
  adsl,
  USUBJID %in% union(rs_subjects, setdiff(adsl_subjects, rs_subjects)[1])
)

## ----eval=TRUE----------------------------------------------------------------
adsl_vars <- exprs(RANDDT, TRTSDT)
adrs <- derive_vars_merged(
  rs,
  dataset_add = adsl,
  new_vars = adsl_vars,
  by_vars = get_admiral_option("subject_keys")
)

## ----echo=TRUE, message=FALSE-------------------------------------------------
param_lookup <- tribble(
  ~RSCAT, ~RSTESTCD, ~RSEVAL,
  ~PARAMCD, ~PARAM, ~PARAMN,
  ~PARCAT1, ~PARCAT1N, ~PARCAT2, ~PARCAT2N,

  # CA-125
  "CA125", "OVRLRESP", "INVESTIGATOR",
  "OVRCA125", "CA-125 Overall Response by Investigator", 1,
  "CA-125", 1, "Investigator", 1,

  # RECIST 1.1
  "RECIST 1.1", "OVRLRESP", "INVESTIGATOR",
  "OVRR11", "RECIST 1.1 Overall Response by Investigator", 2,
  "RECIST 1.1", 2, "Investigator", 1,

  # Combined
  "RECIST 1.1 - CA125", "OVRLRESP", "INVESTIGATOR",
  "OVRR11CA", "Combined Overall Response by Investigator", 3,
  "Combined", 3, "Investigator", 1
)

## ----eval=TRUE, include=TRUE, message=FALSE-----------------------------------
adrs <- derive_vars_merged_lookup(
  adrs,
  dataset_add = param_lookup,
  by_vars = exprs(RSCAT, RSTESTCD, RSEVAL)
)

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs %>%
    arrange(!!!get_admiral_option("subject_keys"), PARAMN, RSSEQ),
  display_vars = exprs(USUBJID, VISIT, RSCAT, RSTESTCD, RSEVAL, PARAMCD, PARAM, PARCAT1, PARCAT2)
)

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  derive_vars_dt(
    dtc = RSDTC,
    new_vars_prefix = "A",
    highest_imputation = "D",
    date_imputation = "last"
  ) %>%
  derive_vars_dy(
    reference_date = TRTSDT,
    source_vars = exprs(ADT)
  ) %>%
  derive_vars_dtm(
    dtc = RSDTC,
    new_vars_prefix = "A",
    highest_imputation = "D",
    date_imputation = "last",
    flag_imputation = "time"
  ) %>%
  mutate(AVISIT = VISIT)

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, PARAMCD, VISIT, AVISIT, RSDTC, ADT, ADTF, ADY)
)

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  mutate(
    AVALC = RSSTRESC,
    AVAL = aval_resp(AVALC)
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, PARAMCD, AVISIT, ADT, AVAL, AVALC)
)

## -----------------------------------------------------------------------------
worst_resp <- function(arg) {
  case_when(
    arg == "NE" ~ 1,
    arg == "CR" ~ 2,
    arg == "PR" ~ 3,
    arg == "SD" ~ 4,
    arg == "NON-CR/NON-PD" ~ 5,
    arg == "PD" ~ 6,
    TRUE ~ 0
  )
}

adrs <- adrs %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = c(get_admiral_option("subject_keys"), exprs(PARAMCD, ADT)),
      order = exprs(worst_resp(AVALC), RSSEQ),
      new_var = ANL01FL,
      mode = "last"
    ),
    filter = !is.na(AVAL) & ADT >= RANDDT
  )

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  derive_var_relative_flag(
    by_vars = c(get_admiral_option("subject_keys"), exprs(PARAMCD)),
    order = exprs(ADT, RSSEQ),
    new_var = ANL02FL,
    condition = (AVALC == "PD" | MOUSEANT == "Y"),
    mode = "first",
    selection = "before",
    inclusive = TRUE
  )

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  select(-CA125EFL) %>%
  derive_var_merged_exist_flag(
    dataset_add = adrs,
    by_vars = get_admiral_option("subject_keys"),
    new_var = CA125EFL,
    condition = (CA125EFL == "Y")
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, AVALC, ADT, ANL01FL, ANL02FL, CA125EFL)
)

## -----------------------------------------------------------------------------
# used for derivation of CA-125 PD
ovr_pd <- filter(adrs, PARAMCD == "OVRCA125" & ANL01FL == "Y" & ANL02FL == "Y")

# used for derivation of CA-125 response parameters
ovr_ca125 <- filter(adrs, PARAMCD == "OVRCA125" & CA125EFL == "Y" & ANL01FL == "Y" & ANL02FL == "Y")

# used for derivation of unconfirmed best overall response from RECIST 1.1 and confirmed CA-125 together
ovr_ubor <- filter(adrs, PARAMCD == "OVRR11CA" & CA125EFL == "Y" & ANL01FL == "Y" & ANL02FL == "Y")

# used for derivation of confirmed best overall response from RECIST 1.1 and confirmed CA-125 together
ovr_r11 <- filter(adrs, PARAMCD == "OVRR11" & CA125EFL == "Y" & ANL01FL == "Y" & ANL02FL == "Y")

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  derive_extreme_records(
    dataset_ref = adsl,
    dataset_add = ovr_pd,
    by_vars = get_admiral_option("subject_keys"),
    filter_add = AVALC == "PD",
    order = exprs(ADT),
    mode = "first",
    keep_source_vars = exprs(everything()),
    set_values_to = exprs(
      PARAMCD = "PDCA125",
      PARAM = "CA-125 Disease Progression by Investigator",
      PARAMN = 4,
      PARCAT1 = "CA-125",
      PARCAT1N = 1,
      PARCAT2 = "Investigator",
      PARCAT2N = 1,
      ANL01FL = "Y",
      ANL02FL = "Y"
    )
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs %>%
    arrange(!!!get_admiral_option("subject_keys"), PARAMN, ADT),
  display_vars = exprs(USUBJID, PARAMCD, AVISIT, ADT, AVALC),
  filter = PARAMCD == "PDCA125"
)

## -----------------------------------------------------------------------------
definition_mcrit <- exprs(
  ~PARAMCD, ~condition,
  ~MCRIT1ML, ~MCRIT1MN,
  "PDCA125", CAELEPRE == "Y" & CANORM2X == "Y",
  "Patients with elevated CA-125 before treatment and normalization of CA-125 (A)", 1,
  "PDCA125", CAELEPRE == "Y" & CNOTNORM == "Y",
  "Patients with elevated CA-125 before treatment, which never normalizes (B)", 2,
  "PDCA125", CAELEPRE == "N" & CANORM2X == "Y",
  "Patients with CA-125 in the reference range before treatment (C)", 3
)

adrs <- adrs %>%
  mutate(MCRIT1 = if_else(PARAMCD == "PDCA125", "PD Category Group", NA_character_)) %>%
  derive_vars_cat(
    definition = definition_mcrit,
    by_vars = exprs(PARAMCD)
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, PARAMCD, AVISIT, ADT, AVALC, MCRIT1, MCRIT1ML, MCRIT1MN),
  filter = PARAMCD == "PDCA125"
)

## -----------------------------------------------------------------------------
bor_sd_gcig <- event(
  description = "Define stable disease (SD) for best overall response (BOR)",
  dataset_name = "ovr",
  condition = AVALC == "SD",
  set_values_to = exprs(AVALC = "SD")
)

bor_ne_gcig <- event(
  description = "Define not evaluable (NE) for best overall response (BOR)",
  dataset_name = "ovr",
  condition = AVALC == "NE",
  set_values_to = exprs(AVALC = "NE")
)

adrs <- adrs %>%
  derive_extreme_event(
    by_vars = get_admiral_option("subject_keys"),
    tmp_event_nr_var = event_nr,
    order = exprs(event_nr, ADT),
    mode = "first",
    source_datasets = list(
      ovr = ovr_ca125,
      adsl = adsl
    ),
    events = list(
      bor_cr, bor_pr, bor_sd_gcig, bor_pd, bor_ne_gcig, no_data_missing
    ),
    set_values_to = exprs(
      PARAMCD = "CBORCA",
      PARAM = "CA-125 Best Confirmed Overall Response by Investigator",
      PARAMN = 5,
      PARCAT1 = "CA-125",
      PARCAT1N = 1,
      PARCAT2 = "Investigator",
      PARCAT2N = 1,
      AVAL = aval_resp(AVALC),
      ANL01FL = "Y",
      ANL02FL = "Y"
    )
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, AVALC, ADT),
  filter = PARAMCD == "CBORCA"
)

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  derive_extreme_event(
    by_vars = get_admiral_option("subject_keys"),
    tmp_event_nr_var = event_nr,
    order = exprs(event_nr, ADT),
    mode = "first",
    source_datasets = list(
      ovr = ovr_ubor,
      adsl = adsl
    ),
    events = list(
      bor_cr, bor_pr, bor_sd_gcig, bor_pd, bor_ne_gcig, no_data_missing
    ),
    set_values_to = exprs(
      PARAMCD = "BORCA11",
      PARAM = "Combined Best Unconfirmed Overall Response by Investigator",
      PARAMN = 6,
      PARCAT1 = "Combined",
      PARCAT1N = 3,
      PARCAT2 = "Investigator",
      PARCAT2N = 1,
      AVAL = aval_resp(AVALC),
      ANL01FL = "Y",
      ANL02FL = "Y"
    )
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, AVALC, ADT),
  filter = PARAMCD == "BORCA11"
)

