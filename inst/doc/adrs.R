## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(admiraldev)

## ----message=FALSE------------------------------------------------------------
library(admiral)
library(admiralonco)
library(dplyr)
library(pharmaversesdtm)
library(pharmaverseadam)
library(lubridate)
library(stringr)
data("adsl")
data("rs_onco_recist")
data("tu_onco_recist")

rs <- rs_onco_recist
tu <- tu_onco_recist

rs <- convert_blanks_to_na(rs)
tu <- convert_blanks_to_na(tu)

## ----echo=FALSE---------------------------------------------------------------
# select subjects from adsl such that there is one subject without RS data
rs_subjects <- unique(rs$USUBJID)
adsl_subjects <- unique(adsl$USUBJID)
adsl <- filter(
  adsl,
  USUBJID %in% union(rs_subjects, setdiff(adsl_subjects, rs_subjects)[1])
)

## ----eval=TRUE----------------------------------------------------------------
adsl_vars <- exprs(RANDDT)
adrs <- derive_vars_merged(
  rs,
  dataset_add = adsl,
  new_vars = adsl_vars,
  by_vars = get_admiral_option("subject_keys")
)

## ----eval=TRUE, echo=FALSE----------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, RSTESTCD, RSDTC, VISIT, RANDDT),
  filter = RSTESTCD == "OVRLRESP"
)

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  filter(RSEVAL == "INVESTIGATOR" & RSTESTCD == "OVRLRESP") %>%
  mutate(
    PARAMCD = "OVR",
    PARAM = "Overall Response by Investigator",
    PARCAT1 = "Tumor Response",
    PARCAT2 = "Investigator",
    PARCAT3 = "RECIST 1.1"
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, VISIT, RSTESTCD, RSEVAL, PARAMCD, PARAM, PARCAT1, PARCAT2, PARCAT3)
)

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  derive_vars_dt(
    dtc = RSDTC,
    new_vars_prefix = "A",
    highest_imputation = "D",
    date_imputation = "last"
  ) %>%
  mutate(AVISIT = VISIT)

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, RSSTRESC, RSDTC, ADT, ADTF)
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
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, RSSTRESC, AVALC, AVAL)
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
      by_vars = c(get_admiral_option("subject_keys"), exprs(ADT)),
      order = exprs(worst_resp(AVALC), RSSEQ),
      new_var = ANL01FL,
      mode = "last"
    ),
    filter = !is.na(AVAL) & ADT >= RANDDT
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, RANDDT, ANL01FL)
)

## ----eval=FALSE---------------------------------------------------------------
#  adrs <- adrs %>%
#    mutate(
#      ANL01FL = case_when(
#        !is.na(AVAL) & ADT >= RANDDT & ADT < NACTDT ~ "Y",
#        TRUE ~ NA_character_
#      )
#    )

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  derive_var_relative_flag(
    by_vars = get_admiral_option("subject_keys"),
    order = exprs(ADT, RSSEQ),
    new_var = ANL02FL,
    condition = AVALC == "PD",
    mode = "first",
    selection = "before",
    inclusive = TRUE
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, AVALC, ADT, ANL01FL, ANL02FL)
)

## -----------------------------------------------------------------------------
ovr <- filter(adrs, PARAMCD == "OVR" & ANL01FL == "Y" & ANL02FL == "Y")

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  ovr,
  display_vars = exprs(USUBJID, AVISIT, AVALC, ADT, RANDDT)
)

## -----------------------------------------------------------------------------
confirmation_period <- 21

crsp_y_cr <- event_joined(
  description = paste(
    "Define confirmed response as CR followed by CR at least",
    confirmation_period,
    "days later and at most one NE in between"
  ),
  dataset_name = "ovr",
  join_vars = exprs(AVALC, ADT),
  join_type = "after",
  order = exprs(ADT),
  first_cond_upper = AVALC.join == "CR" &
    ADT.join >= ADT + days(confirmation_period),
  condition = AVALC == "CR" &
    all(AVALC.join %in% c("CR", "NE")) &
    count_vals(var = AVALC.join, val = "NE") <= 1,
  set_values_to = exprs(AVALC = "Y")
)

crsp_y_pr <- event_joined(
  description = paste(
    "Define confirmed response as PR followed by CR or PR at least",
    confirmation_period,
    "days later, at most one NE in between, and no PR after CR"
  ),
  dataset_name = "ovr",
  join_vars = exprs(AVALC, ADT),
  join_type = "after",
  order = exprs(ADT),
  first_cond_upper = AVALC.join %in% c("CR", "PR") &
    ADT.join >= ADT + days(confirmation_period),
  condition = AVALC == "PR" &
    all(AVALC.join %in% c("CR", "PR", "NE")) &
    count_vals(var = AVALC.join, val = "NE") <= 1 &
    (
      min_cond(
        var = ADT.join,
        cond = AVALC.join == "CR"
      ) > max_cond(var = ADT.join, cond = AVALC.join == "PR") |
        count_vals(var = AVALC.join, val = "CR") == 0 |
        count_vals(var = AVALC.join, val = "PR") == 0
    ),
  set_values_to = exprs(AVALC = "Y")
)

cbor_cr <- event_joined(
  description = paste(
    "Define complete response (CR) for confirmed best overall response (CBOR) as",
    "CR followed by CR at least",
    confirmation_period,
    "days later and at most one NE in between"
  ),
  dataset_name = "ovr",
  join_vars = exprs(AVALC, ADT),
  join_type = "after",
  first_cond_upper = AVALC.join == "CR" &
    ADT.join >= ADT + confirmation_period,
  condition = AVALC == "CR" &
    all(AVALC.join %in% c("CR", "NE")) &
    count_vals(var = AVALC.join, val = "NE") <= 1,
  set_values_to = exprs(AVALC = "CR")
)

cbor_pr <- event_joined(
  description = paste(
    "Define partial response (PR) for confirmed best overall response (CBOR) as",
    "PR followed by CR or PR at least",
    confirmation_period,
    "28 days later, at most one NE in between, and no PR after CR"
  ),
  dataset_name = "ovr",
  join_vars = exprs(AVALC, ADT),
  join_type = "after",
  first_cond_upper = AVALC.join %in% c("CR", "PR") &
    ADT.join >= ADT + confirmation_period,
  condition = AVALC == "PR" &
    all(AVALC.join %in% c("CR", "PR", "NE")) &
    count_vals(var = AVALC.join, val = "NE") <= 1 &
    (
      min_cond(
        var = ADT.join,
        cond = AVALC.join == "CR"
      ) > max_cond(var = ADT.join, cond = AVALC.join == "PR") |
        count_vals(var = AVALC.join, val = "CR") == 0 |
        count_vals(var = AVALC.join, val = "PR") == 0
    ),
  set_values_to = exprs(AVALC = "PR")
)

no_data_n <- event(
  description = "Define no response for all patients in adsl (should be used as last event)",
  dataset_name = "adsl",
  condition = TRUE,
  set_values_to = exprs(AVALC = "N"),
  keep_source_vars = adsl_vars
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

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  derive_extreme_records(
    dataset_ref = adsl,
    dataset_add = adrs,
    by_vars = get_admiral_option("subject_keys"),
    filter_add = PARAMCD == "OVR" & AVALC == "PD" & ANL01FL == "Y",
    order = exprs(ADT, RSSEQ),
    mode = "first",
    exist_flag = AVALC,
    false_value = "N",
    set_values_to = exprs(
      PARAMCD = "PD",
      PARAM = "Disease Progression by Investigator",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "RECIST 1.1",
      AVAL = yn_to_numeric(AVALC),
      ANL01FL = "Y"
    )
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, ANL01FL),
  filter = PARAMCD == "PD"
)

## -----------------------------------------------------------------------------
rsp_y

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  derive_extreme_event(
    by_vars = get_admiral_option("subject_keys"),
    order = exprs(event_nr, ADT),
    tmp_event_nr_var = event_nr,
    mode = "first",
    events = list(rsp_y, no_data_n),
    source_datasets = list(
      ovr = ovr,
      adsl = adsl
    ),
    set_values_to = exprs(
      PARAMCD = "RSP",
      PARAM = "Response by Investigator (confirmation not required)",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "RECIST 1.1",
      AVAL = yn_to_numeric(AVALC),
      ANL01FL = "Y"
    )
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, ANL01FL),
  filter = PARAMCD == "RSP"
)

## -----------------------------------------------------------------------------
cb_y

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  derive_extreme_event(
    by_vars = get_admiral_option("subject_keys"),
    order = exprs(desc(AVALC), ADT, event_nr),
    tmp_event_nr_var = event_nr,
    mode = "first",
    events = list(rsp_y, cb_y, no_data_n),
    source_datasets = list(
      ovr = ovr,
      adsl = adsl
    ),
    set_values_to = exprs(
      PARAMCD = "CB",
      PARAM = "Clinical Benefit by Investigator (confirmation for response not required)",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "RECIST 1.1",
      AVAL = yn_to_numeric(AVALC),
      ANL01FL = "Y"
    )
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, RANDDT, ANL01FL),
  filter = PARAMCD == "CB"
)

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  derive_extreme_event(
    by_vars = get_admiral_option("subject_keys"),
    order = exprs(event_nr, ADT),
    tmp_event_nr_var = event_nr,
    mode = "first",
    source_datasets = list(
      ovr = ovr,
      adsl = adsl
    ),
    events = list(bor_cr, bor_pr, bor_sd, bor_non_crpd, bor_pd, bor_ne, no_data_missing),
    set_values_to = exprs(
      PARAMCD = "BOR",
      PARAM = "Best Overall Response by Investigator (confirmation not required)",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "RECIST 1.1",
      AVAL = aval_resp(AVALC),
      ANL01FL = "Y"
    )
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, RANDDT, ANL01FL),
  filter = PARAMCD == "BOR"
)

## ----eval=FALSE---------------------------------------------------------------
#  aval_resp_new <- function(arg) {
#    case_when(
#      arg == "CR" ~ 7,
#      arg == "PR" ~ 6,
#      arg == "SD" ~ 5,
#      arg == "NON-CR/NON-PD" ~ 4,
#      arg == "PD" ~ 3,
#      arg == "NE" ~ 2,
#      arg == "MISSING" ~ 1,
#      TRUE ~ NA_real_
#    )
#  }

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  derive_extreme_records(
    dataset_ref = adsl,
    dataset_add = adrs,
    by_vars = get_admiral_option("subject_keys"),
    filter_add = PARAMCD == "BOR" & AVALC %in% c("CR", "PR"),
    exist_flag = AVALC,
    false_value = "N",
    set_values_to = exprs(
      PARAMCD = "BCP",
      PARAM = "Best Overall Response of CR/PR by Investigator (confirmation not required)",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "RECIST 1.1",
      AVAL = yn_to_numeric(AVALC),
      ANL01FL = "Y"
    )
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, ANL01FL),
  filter = PARAMCD == "BCP"
)

## -----------------------------------------------------------------------------
crsp_y_cr

## -----------------------------------------------------------------------------
crsp_y_pr

## -----------------------------------------------------------------------------
cbor_cr

## -----------------------------------------------------------------------------
cbor_pr

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  derive_extreme_event(
    by_vars = get_admiral_option("subject_keys"),
    order = exprs(desc(AVALC), ADT, event_nr),
    tmp_event_nr_var = event_nr,
    mode = "first",
    source_datasets = list(
      ovr = ovr,
      adsl = adsl
    ),
    events = list(crsp_y_cr, crsp_y_pr, no_data_n),
    set_values_to = exprs(
      PARAMCD = "CRSP",
      PARAM = "Confirmed Response by Investigator",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "RECIST 1.1",
      AVAL = yn_to_numeric(AVALC),
      ANL01FL = "Y"
    )
  )

adrs <- adrs %>%
  derive_extreme_event(
    by_vars = get_admiral_option("subject_keys"),
    order = exprs(desc(AVALC), ADT, event_nr),
    tmp_event_nr_var = event_nr,
    mode = "first",
    events = list(crsp_y_cr, crsp_y_pr, cb_y, no_data_n),
    source_datasets = list(
      ovr = ovr,
      adsl = adsl
    ),
    set_values_to = exprs(
      PARAMCD = "CCB",
      PARAM = "Confirmed Clinical Benefit by Investigator",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "RECIST 1.1",
      AVAL = yn_to_numeric(AVALC),
      ANL01FL = "Y"
    )
  )

adrs <- adrs %>%
  derive_extreme_event(
    by_vars = get_admiral_option("subject_keys"),
    order = exprs(event_nr, ADT),
    tmp_event_nr_var = event_nr,
    mode = "first",
    events = list(cbor_cr, cbor_pr, bor_sd, bor_non_crpd, bor_pd, bor_ne, no_data_missing),
    source_datasets = list(
      ovr = ovr,
      adsl = adsl
    ),
    set_values_to = exprs(
      PARAMCD = "CBOR",
      PARAM = "Best Confirmed Overall Response by Investigator",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "RECIST 1.1",
      AVAL = aval_resp(AVALC),
      ANL01FL = "Y"
    )
  ) %>%
  derive_extreme_records(
    dataset_ref = adsl,
    dataset_add = adrs,
    by_vars = get_admiral_option("subject_keys"),
    filter_add = PARAMCD == "CBOR" & AVALC %in% c("CR", "PR"),
    exist_flag = AVALC,
    false_value = "N",
    set_values_to = exprs(
      PARAMCD = "CBCP",
      PARAM = "Best Confirmed Overall Response of CR/PR by Investigator",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "RECIST 1.1",
      AVAL = yn_to_numeric(AVALC),
      ANL01FL = "Y"
    )
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, RANDDT, ANL01FL),
  filter = PARAMCD %in% c("CRSP", "CCB", "CBOR", "CBCP")
)

## -----------------------------------------------------------------------------
cb_y_pd <- event(
  description = paste(
    "Define PD occuring more than 42 days after",
    "randomization as clinical benefit"
  ),
  dataset_name = "ovr",
  condition = AVALC == "PD" & ADT > RANDDT + 42,
  set_values_to = exprs(AVALC = "Y")
)

adrs <- adrs %>%
  derive_extreme_event(
    by_vars = get_admiral_option("subject_keys"),
    order = exprs(desc(AVALC), ADT, event_nr),
    tmp_event_nr_var = event_nr,
    mode = "first",
    events = list(crsp_y_cr, crsp_y_pr, cb_y, cb_y_pd, no_data_n),
    source_datasets = list(
      ovr = ovr,
      adsl = adsl
    ),
    set_values_to = exprs(
      PARAMCD = "ACCB",
      PARAM = "Alternative Confirmed Clinical Benefit by Investigator",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "RECIST 1.1",
      AVAL = yn_to_numeric(AVALC),
      ANL01FL = "Y"
    )
  )

## -----------------------------------------------------------------------------
bor_ned <- event(
  description = paste(
    "Define no evidence of disease (NED) for best overall response (BOR) as NED",
    "occuring at least 42 days after randomization"
  ),
  dataset_name = "ovr",
  condition = AVALC == "NED" & ADT >= RANDDT + 42,
  set_values_to = exprs(AVALC = "NED")
)

adrs <- adrs %>%
  derive_extreme_event(
    by_vars = get_admiral_option("subject_keys"),
    order = exprs(event_nr, ADT),
    tmp_event_nr_var = event_nr,
    mode = "first",
    source_datasets = list(
      ovr = ovr,
      adsl = adsl
    ),
    events = list(bor_cr, bor_pr, bor_sd, bor_non_crpd, bor_ned, bor_pd, bor_ne, no_data_missing),
    set_values_to = exprs(
      PARAMCD = "A1BOR",
      PARAM = paste(
        "Best Overall Response by Investigator (confirmation not required)",
        "- RECIST 1.1 adjusted for NED at Baseline"
      ),
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "RECIST 1.1 adjusted for NED at Baseline",
      AVAL = aval_resp(AVALC),
      ANL01FL = "Y"
    )
  )

## -----------------------------------------------------------------------------
adrs_bicr <- rs %>%
  filter(
    RSEVAL == "INDEPENDENT ASSESSOR" & RSACPTFL == "Y" & RSTESTCD == "OVRLRESP"
  ) %>%
  mutate(
    PARAMCD = "OVRB",
    PARAM = "Overall Response by BICR",
    PARCAT1 = "Tumor Response",
    PARCAT2 = "Blinded Independent Central Review",
    PARCAT3 = "RECIST 1.1"
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs_bicr,
  display_vars = exprs(USUBJID, VISIT, RSTESTCD, RSEVAL, PARAMCD, PARAM, PARCAT1, PARCAT2, PARCAT3),
  filter = PARAMCD == "OVRB"
)

## -----------------------------------------------------------------------------
adsldth <- adsl %>%
  select(!!!get_admiral_option("subject_keys"), DTHDT, !!!adsl_vars)

adrs <- adrs %>%
  derive_extreme_records(
    dataset_ref = adsldth,
    dataset_add = adsldth,
    by_vars = get_admiral_option("subject_keys"),
    filter_add = !is.na(DTHDT),
    exist_flag = AVALC,
    false_value = "N",
    set_values_to = exprs(
      PARAMCD = "DEATH",
      PARAM = "Death",
      PARCAT1 = "Reference Event",
      AVAL = yn_to_numeric(AVALC),
      ANL01FL = "Y",
      ADT = DTHDT
    )
  ) %>%
  select(-DTHDT)

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, ANL01FL),
  filter = PARAMCD == "DEATH"
)

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  derive_extreme_records(
    dataset_ref = adsl,
    dataset_add = adrs,
    by_vars = get_admiral_option("subject_keys"),
    filter_add = PARAMCD == "OVR" & ANL01FL == "Y",
    order = exprs(ADT, RSSEQ),
    mode = "last",
    set_values_to = exprs(
      PARAMCD = "LSTA",
      PARAM = "Last Disease Assessment by Investigator",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "RECIST 1.1",
      ANL01FL = "Y"
    )
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, ANL01FL),
  filter = PARAMCD == "LSTA"
)

## -----------------------------------------------------------------------------
adslmdis <- adsl %>%
  select(!!!get_admiral_option("subject_keys"), !!!adsl_vars)

adrs <- adrs %>%
  derive_param_exist_flag(
    dataset_ref = adslmdis,
    dataset_add = tu,
    condition = TUEVAL == "INVESTIGATOR" & TUSTRESC == "TARGET" & VISIT == "SCREENING",
    false_value = "N",
    missing_value = "N",
    set_values_to = exprs(
      PARAMCD = "MDIS",
      PARAM = "Measurable Disease at Baseline by Investigator",
      PARCAT2 = "Investigator",
      PARCAT3 = "RECIST 1.1",
      AVAL = yn_to_numeric(AVALC),
      ANL01FL = "Y"
    )
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, ANL01FL),
  filter = PARAMCD == "MDIS"
)

## ----eval=TRUE----------------------------------------------------------------
adrs <- adrs %>%
  derive_var_obs_number(
    by_vars = get_admiral_option("subject_keys"),
    order = exprs(PARAMCD, ADT, VISITNUM, RSSEQ),
    check_type = "error"
  )

## ----eval=TRUE, echo=FALSE----------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, PARAMCD, ADT, VISITNUM, AVISIT, ASEQ),
  filter = USUBJID == "01-701-1015"
)

## ----eval=TRUE----------------------------------------------------------------
adrs <- adrs %>%
  derive_vars_merged(
    dataset_add = select(adsl, !!!negate_vars(adsl_vars)),
    by_vars = get_admiral_option("subject_keys")
  )

## ----eval=TRUE, echo=FALSE----------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, RFSTDTC, RFENDTC, DTHDTC, DTHFL, AGE, AGEU),
  filter = USUBJID == "01-701-1015"
)

