## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(admiraldev)
library(gt)

## ----warning=FALSE, message=FALSE---------------------------------------------
library(admiral)
library(admiralonco)
library(pharmaversesdtm)
library(pharmaverseadam)
library(dplyr)
library(tibble)

## ----message=FALSE------------------------------------------------------------
# PCWG3 SDTM data
rs <- pharmaversesdtm::rs_onco_pcwg3
rs <- convert_blanks_to_na(rs)

# Exclude PSA records
rs <- rs %>% filter(RSTEST != "Tumor Marker Response")

# ADaM data
adsl <- pharmaverseadam::adsl

## ----echo=FALSE---------------------------------------------------------------
# select subjects from adsl such that there is one subject without RS data
rs_subjects <- unique(rs$USUBJID)
adsl_subjects <- unique(adsl$USUBJID)
adsl <- filter(
  adsl,
  USUBJID %in% union(rs_subjects, setdiff(adsl_subjects, rs_subjects)[1])
)

## ----eval=TRUE, echo=FALSE----------------------------------------------------
################### To be Removed in next release ######
rs <- rs %>%
  mutate(
    RSSTRESC = case_when(
      USUBJID == "01-701-1275" & VISIT == "WEEK 8" & RSTESTCD %in% c("BONERESP", "OVRLRESP") ~ "PD",
      USUBJID == "01-701-1097" & VISIT == "WEEK 16" & RSTESTCD %in% c("BONERESP", "OVRLRESP") ~ "PD",
      TRUE ~ RSSTRESC
    )
  )
########################

dataset_vignette(
  rs,
  display_vars = exprs(USUBJID, RSCAT, RSTESTCD, RSSTRESC, VISIT, VISITNUM, RSDTC)
)

## ----eval=TRUE----------------------------------------------------------------
adsl_vars <- exprs(TRTSDT)
adrs <- derive_vars_merged(
  rs,
  dataset_add = adsl,
  new_vars = adsl_vars,
  by_vars = get_admiral_option("subject_keys")
)

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  derive_vars_dtm(
    dtc = RSDTC,
    new_vars_prefix = "A",
    highest_imputation = "D",
    date_imputation = "last"
  ) %>%
  derive_vars_dtm_to_dt(exprs(ADTM)) %>%
  derive_vars_dy(
    reference_date = TRTSDT,
    source_vars = exprs(ADT)
  ) %>%
  mutate(
    AVISIT = VISIT,
    AVISITN = VISITNUM
  )

## ----eval=TRUE, include=TRUE, message=FALSE-----------------------------------
# Prepare param_lookup for SDTM RSTESTCD to add metadata
param_lookup <- tibble::tribble(
  ~RSTESTCD,  ~PARAMCD,   ~PARAM,                                   ~PARAMN,
  "SFTSRESP", "SFTSRESP", "Soft Tissue Response by Investigator",         1,
  "BONERESP", "BONERESP", "Bone Response by Investigator",                2,
  "OVRLRESP", "OVRLRESP", "Overall Tumor Response by Investigator",       3
)

adrs <- adrs %>%
  derive_vars_merged_lookup(
    dataset_add = param_lookup,
    by_vars = exprs(RSTESTCD)
  ) %>%
  mutate(
    PARCAT1 = RSCAT,
    AVALC = RSSTRESC
  )

## ----eval=TRUE, include=TRUE, message=FALSE,echo=FALSE------------------------
overall_tpr_table <- tribble(
  ~`Soft Tissue (RECIST 1.1) TPR`, ~`Bone Lesion (PCWG3) TPR`, ~`Overall PCWG TPR`,
  "PD",                            "Any",                      "PD",
  "Any",                           "PD",                       "PD",
  "NE",                            "Non-PD, PDu, NED or NE",   "NE",
  "NED",                           "Non-PD",                   "Non-CR/Non-PD",
  "NED",                           "PDu",                      "PDu",
  "NED",                           "NED",                      "NE",
  "NED",                           "NE",                       "NE",
  "SD",                            "Non-PD, PDu, NED or NE",   "SD",
  "Non-CR/Non-PD",                 "Non-PD, PDu, NED or NE",   "Non-CR/Non-PD",
  "PR",                            "Non-PD, PDu, NED or NE",   "PR",
  "CR",                            "Non-PD, PDu, or NE",       "PR (1)",
  "CR",                            "Non-PD, PDu, or NE",       "Non-CR/Non-PD (2)",
  "CR",                            "NED",                      "CR"
)

overall_tpr_table %>%
  gt() %>%
  tab_header(
    title = "Table 1: Overall Time Point Response",
    subtitle = "Soft Tissue (RECIST 1.1) TPR, Bone Lesion (PCWG3) TPR, and PCWG Combined TPR"
  ) %>%
  cols_label(
    `Soft Tissue (RECIST 1.1) TPR` = "Soft Tissue (RECIST 1.1)",
    `Bone Lesion (PCWG3) TPR` = "Bone Lesion (PCWG3)",
    `Overall PCWG TPR` = "Overall PCWG"
  ) %>%
  tab_footnote(
    footnote = "* When no target and non-target lesions are identified at baseline, and no new lesions are identified on-study, the response will be No Evidence of Disease (NED)."
  ) %>%
  tab_footnote(
    footnote = "** Progressive Disease Unconfirmed (PDu): Temporary marker of possible PD where at least 2 new bone lesions are present, but an additional scan is required for confirmation. To be updated to PD or Non-PD once a subsequent scan is available. If this is the final visit, the response remains as PDu."
  ) %>%
  tab_footnote(
    footnote = "(1) The overall TPR will be PR if target lesions were present at screening."
  ) %>%
  tab_footnote(
    footnote = "(2) The overall TPR will be Non-CR/Non-PD if no target lesions were present at screening."
  )

## ----eval=TRUE, message=FALSE, include=TRUE-----------------------------------
adrs <- derive_param_computed(
  dataset = adrs,
  by_vars = exprs(
    !!!get_admiral_option("subject_keys"), !!!adsl_vars, DOMAIN, RSEVAL, ADT,
    ADY, ADTM, ADTF, VISIT, VISITNUM, AVISIT, AVISITN
  ),
  parameters = c("SFTSRESP", "BONERESP"),
  set_values_to = exprs(
    AVALC = case_when(
      # Scenario 1 & 2: Soft Tissue PD or Bone Lesion PD -> Overall response = PD
      AVALC.SFTSRESP == "PD" | AVALC.BONERESP == "PD" ~ "PD",

      # Scenario 3: Soft Tissue = NE + Bone Lesion = NON-PD, PDu, NED, or NE -> Overall response = NE
      AVALC.SFTSRESP == "NE" & AVALC.BONERESP %in% c("NON-PD", "PDu", "NED", "NE") ~ "NE",

      # Scenario 4: Soft Tissue = NED + Bone Lesion = NON-PD -> Overall response = Non-CR/NON-PD
      AVALC.SFTSRESP == "NED" & AVALC.BONERESP == "NON-PD" ~ "Non-CR/NON-PD",

      # Scenario 5: Soft Tissue = NED + Bone Lesion = PDu -> Overall response = PDu
      AVALC.SFTSRESP == "NED" & AVALC.BONERESP == "PDu" ~ "PDu",

      # Scenario 6: Soft Tissue = NED + Bone Lesion = NED -> Overall response = NE
      AVALC.SFTSRESP == "NED" & AVALC.BONERESP == "NED" ~ "NE",

      # Scenario 7: Soft Tissue = NED + Bone Lesion = NE -> Overall response = NE
      AVALC.SFTSRESP == "NED" & AVALC.BONERESP == "NE" ~ "NE",

      # Scenario 8: Soft Tissue = SD + Bone Lesion = NON-PD, PDu, NED, or NE -> Overall response = SD
      AVALC.SFTSRESP == "SD" & AVALC.BONERESP %in% c("NON-PD", "PDu", "NED", "NE") ~ "SD",

      # Scenario 9: Soft Tissue = Non-CR/NON-PD + Bone Lesion = NON-PD, PDu, NED, or NE -> Overall response = Non-CR/NON-PD
      AVALC.SFTSRESP == "Non-CR/NON-PD" & AVALC.BONERESP %in% c("NON-PD", "PDu", "NED", "NE") ~ "Non-CR/NON-PD",

      # Scenario 10: Soft Tissue = PR + Bone Lesion = NON-PD, PDu, NED, or NE -> Overall response = PR
      AVALC.SFTSRESP == "PR" & AVALC.BONERESP %in% c("NON-PD", "PDu", "NED", "NE") ~ "PR",

      # Scenario 11: Soft Tissue = CR + Bone Lesion = NON-PD, PDu, NE -> Overall response = PR
      # ((1) The overall TPR will be PR if target lesions were present at screening.)
      AVALC.SFTSRESP == "CR" & AVALC.BONERESP %in% c("NON-PD", "PDu", "NE") ~ "PR",

      # Soft Tissue = CR + Bone Lesion = NON-PD, PDu, NE -> Overall response =Non-CR/NON-PD
      # (2) The overall TPR will be Non-CR/NON-PD if no target lesions were present at screening.)
      # AVALC.SFTSRESP == "CR" & AVALC.BONERESP %in% c("NON-PD", "PDu", "NE") ~ "Non-CR/NON-PD",

      # Scenario 12: Soft Tissue = CR + Bone Lesion = NED -> Overall response = CR
      AVALC.SFTSRESP == "CR" & AVALC.BONERESP == "NED" ~ "CR",

      # Default: If conditions are not met, assign NA
      TRUE ~ NA_character_
    ),
    PARAMCD = "OVRLRESC",
    PARAM = "Overall Tumor Response by Investigator - Derived",
    PARAMN = 4,
    PARCAT1 = "PCWG3 and RECIST 1.1"
  )
)

## ----eval=TRUE, echo=FALSE----------------------------------------------------
dataset_vignette(
  adrs %>%
    arrange(!!!get_admiral_option("subject_keys"), AVISITN, PARAMN),
  display_vars = exprs(USUBJID, PARAM, PARAMCD, PARCAT1, AVALC, AVISIT, ADT)
)

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  mutate(
    AVAL = case_when(
      AVALC == "CR" ~ 1, # Complete Response
      AVALC == "PR" ~ 2, # Partial Response
      AVALC == "SD" ~ 3, # Stable Disease
      AVALC == "PD" ~ 4, # Progressive Disease
      AVALC == "Non-CR/NON-PD" ~ 5, # Neither Complete Response nor Progressive Disease
      AVALC == "NON-PD" ~ 6, # Non-Progressive Disease
      AVALC == "PDu" ~ 7, # Progressive Disease Unconfirmed
      AVALC == "NE" ~ 8, # Not Evaluable
      AVALC == "NED" ~ 9, # No Evidence of Disease
      TRUE ~ NA_real_ # Default for unexpected/missing AVALC values
    )
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs %>%
    arrange(!!!get_admiral_option("subject_keys"), AVISITN, PARAMN),
  display_vars = exprs(USUBJID, PARAMCD, PARAM, AVISIT, ADT, AVALC, AVAL)
)

## ----eval=TRUE, include=TRUE, message=FALSE-----------------------------------
bor_cr <- event(
  description = "Complete Response (CR)",
  dataset_name = "adrs",
  condition = AVALC == "CR",
  set_values_to = exprs(AVALC = "CR")
)

bor_pr <- event(
  description = "Partial Response (PR)",
  dataset_name = "adrs",
  condition = AVALC == "PR",
  set_values_to = exprs(AVALC = "PR")
)

bor_non_crpd <- event(
  description = "Non-CR/Non-PD",
  dataset_name = "adrs",
  condition = AVALC == "Non-CR/NON-PD",
  set_values_to = exprs(AVALC = "Non-CR/Non-PD")
)

bor_sd <- event(
  description = "Stable Disease (SD)",
  dataset_name = "adrs",
  # CR and PR are included for CBOR when CR or PR couldn't be confirmed
  # PDu can occur only as last assessment and is considered as SD
  condition = AVALC %in% c("CR", "PR", "SD", "PDu"),
  set_values_to = exprs(AVALC = "SD")
)

bor_pd <- event(
  description = "Progressive Disease (PD)",
  dataset_name = "adrs",
  condition = AVALC == "PD",
  set_values_to = exprs(AVALC = "PD")
)

bor_ne <- event(
  description = "Not Evaluable (NE)",
  dataset_name = "adrs",
  condition = AVALC == "NE",
  set_values_to = exprs(AVALC = "NE")
)

bor_ned <- event(
  description = "No Evidence of Disease (NED)",
  dataset_name = "adrs",
  condition = AVALC == "NED",
  set_values_to = exprs(AVALC = "NED")
)

no_data_missing <- event(
  description = paste(
    "Define missing response (MISSING) for all patients in adsl (should be used",
    "as last event)"
  ),
  dataset_name = "adsl",
  condition = TRUE,
  set_values_to = exprs(AVALC = "MISSING"),
  keep_source_vars = adsl_vars
)

## ----eval=TRUE, include=TRUE, message=FALSE-----------------------------------
adrs <- adrs %>%
  derive_extreme_event(
    by_vars = get_admiral_option("subject_keys"),
    events = list(
      bor_cr, bor_pr, bor_sd, bor_non_crpd, bor_pd, bor_ne, bor_ned, no_data_missing
    ),
    source_datasets = list(
      adsl = adsl,
      adrs = adrs %>% filter(PARAMCD == "OVRLRESC") # Use derived responses (OVRLRESC)
    ),
    order = exprs(event_nr, ADT), # Prioritize earliest valid event
    tmp_event_nr_var = event_nr,
    mode = "first", # Retain the best response observed at the first occurrence
    set_values_to = exprs(
      PARAMCD = "BOR",
      PARAM = "Best Overall Response",
      PARAMN = 5,
      PARCAT1 = "PCWG3 and RECIST 1.1"
    )
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs %>%
    filter(PARAMCD == "BOR"),
  display_vars = exprs(USUBJID, PARAM, PARAMCD, AVISIT, AVISITN, ADT, AVALC, AVAL)
)

## ----eval=TRUE, include=TRUE, message=FALSE-----------------------------------
# Confirmed CR Event with 28-day persistence
cbor_cr <- event_joined(
  description = "Confirmed Complete Response (CR)",
  dataset_name = "adrs",
  join_vars = exprs(AVALC, ADT),
  join_type = "after",
  first_cond_upper = AVALC.join == "CR" & ADT.join >= ADT + 28, # Follow-up within 28-day window
  condition = AVALC == "CR" & all(AVALC.join == "CR"), # All linked records must also be CR
  set_values_to = exprs(AVALC = "CR") # Set response as Confirmed CR
)

# Confirmed PR Event with 28-day persistence
cbor_pr <- event_joined(
  description = "Confirmed Partial Response (PR)",
  dataset_name = "adrs",
  join_vars = exprs(AVALC, ADT),
  join_type = "after",
  first_cond_upper = AVALC.join %in% c("CR", "PR") & ADT.join >= ADT + 28, # Include CR as confirmation
  condition = AVALC == "PR" & all(AVALC.join %in% c("CR", "PR")), # Ensure no events other than CR or PR in between
  set_values_to = exprs(AVALC = "PR")
)

adrs <- adrs %>%
  derive_extreme_event(
    by_vars = get_admiral_option("subject_keys"),
    events = list(
      cbor_cr, cbor_pr, bor_sd, bor_non_crpd, bor_pd, bor_ne, bor_ned, no_data_missing
    ),
    source_datasets = list(
      adsl = adsl,
      adrs = adrs %>% filter(PARAMCD == "OVRLRESC")
    ),
    tmp_event_nr_var = event_nr,
    order = exprs(event_nr, ADT),
    mode = "first",
    set_values_to = exprs(
      PARAMCD = "CBOR",
      PARAM = "Confirmed Best Overall Response",
      PARAMN = 6,
      PARCAT1 = "PCWG3 and RECIST 1.1"
    )
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs %>%
    filter(PARAMCD == "CBOR"),
  display_vars = exprs(USUBJID, PARAM, PARAMCD, AVISIT, AVISITN, ADT, AVALC, AVAL)
)

