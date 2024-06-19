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
library(pharmaverseadam)
library(pharmaversesdtm)
library(lubridate)
library(stringr)
data("adsl")
# iRECIST oncology sdtm data
data("rs_onco_irecist")

rs <- rs_onco_irecist

rs <- convert_blanks_to_na(rs)

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
    PARCAT3 = "iRECIST"
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
aval_resp_new <- function(arg) {
  case_when(
    arg == "NE" ~ 8,
    arg == "MISSING" ~ 7,
    arg == "iCR" ~ 6,
    arg == "iPR" ~ 5,
    arg == "iSD" ~ 4,
    arg == "NON-iCR/NON-iUPD" ~ 3,
    arg == "iUPD" ~ 2,
    arg == "iCPD" ~ 1,
    TRUE ~ NA_real_
  )
}

adrs <- adrs %>%
  mutate(
    AVALC = RSSTRESC,
    AVAL = aval_resp_new(AVALC)
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, RSSTRESC, AVALC, AVAL)
)

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = c(get_admiral_option("subject_keys"), exprs(ADT)),
      order = exprs(AVAL, RSSEQ),
      new_var = ANL01FL,
      mode = "first"
    ),
    filter = !is.na(AVAL) & AVALC != "MISSING" & ADT >= RANDDT
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
    condition = AVALC == "iCPD",
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
icpd_y <- event_joined(
  description = paste(
    "Define confirmed progressive disease (iCPD) as",
    "iUPD followed by iCPD with only other iUPD and NE responses in between"
  ),
  dataset_name = "ovr",
  join_vars = exprs(AVALC, ADT),
  join_type = "after",
  first_cond_upper = AVALC.join == "iCPD",
  condition = AVALC == "iUPD" &
    all(AVALC.join %in% c("iCPD", "iUPD", "NE")),
  set_values_to = exprs(AVALC = "Y")
)

iupd_y <- event_joined(
  description = paste(
    "Define unconfirmed progressive disease (iUPD) as",
    "iUPD followed only by other iUPD or NE responses"
  ),
  dataset_name = "ovr",
  join_vars = exprs(AVALC, ADT),
  join_type = "all",
  condition = ADT <= ADT.join & AVALC == "iUPD" & all(AVALC.join %in% c("iUPD", "NE")),
  set_values_to = exprs(AVALC = "Y")
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
irsp_y <- event(
  description = "Define CR, iCR, PR, or iPR as (unconfirmed) response",
  dataset_name = "ovr",
  condition = AVALC %in% c("CR", "iCR", "PR", "iPR"),
  set_values_to = exprs(AVALC = "Y")
)

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  derive_extreme_event(
    by_vars = get_admiral_option("subject_keys"),
    order = exprs(ADT),
    mode = "first",
    source_datasets = list(
      ovr = ovr,
      adsl = adsl
    ),
    events = list(icpd_y, no_data_n),
    set_values_to = exprs(
      PARAMCD = "ICPD",
      PARAM = "iRECIST Confirmation of Disease Progression by Investigator",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "iRECIST",
      AVAL = yn_to_numeric(AVALC),
      ANL01FL = "Y"
    )
  )

ovr_orig <- ovr
ovr <- ovr %>%
  group_by(!!!get_admiral_option("subject_keys")) %>%
  filter(ADT >= max_cond(var = ADT, cond = AVALC == "iUPD")) %>%
  ungroup(!!!get_admiral_option("subject_keys"))

adrs <- adrs %>%
  derive_extreme_event(
    by_vars = get_admiral_option("subject_keys"),
    order = exprs(ADT),
    mode = "first",
    source_datasets = list(
      ovr = ovr,
      adsl = adsl
    ),
    events = list(iupd_y, no_data_n),
    set_values_to = exprs(
      PARAMCD = "IUPD",
      PARAM = "iRECIST Unconfirmed Disease Progression by Investigator",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "iRECIST",
      AVAL = yn_to_numeric(AVALC),
      ANL01FL = "Y"
    )
  )
ovr <- ovr_orig

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, ANL01FL),
  filter = PARAMCD %in% c("ICPD", "IUPD")
)

## -----------------------------------------------------------------------------
irsp_y <- event(
  description = "Define iCR or iPR as (unconfirmed) response",
  dataset_name = "ovr",
  condition = AVALC %in% c("iCR", "iPR"),
  set_values_to = exprs(AVALC = "Y")
)

adrs <- adrs %>%
  derive_extreme_event(
    by_vars = get_admiral_option("subject_keys"),
    order = exprs(ADT),
    mode = "first",
    events = list(irsp_y, no_data_n),
    source_datasets = list(
      ovr = ovr,
      adsl = adsl
    ),
    set_values_to = exprs(
      PARAMCD = "IRSP",
      PARAM = "iRECIST Response by Investigator (confirmation not required)",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "iRECIST",
      AVAL = yn_to_numeric(AVALC),
      ANL01FL = "Y"
    )
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, ANL01FL),
  filter = PARAMCD == "IRSP"
)

## -----------------------------------------------------------------------------
icb_y <- event(
  description = paste(
    "Define iCR, iPR, iSD, or NON-iCR/NON-iUPD occuring at least 42 days after",
    "randomization as clinical benefit"
  ),
  dataset_name = "ovr",
  condition = AVALC %in% c("iCR", "iPR", "iSD", "NON-iCR/NON-iUPD") &
    ADT >= RANDDT + 42,
  set_values_to = exprs(AVALC = "Y")
)

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  derive_extreme_event(
    by_vars = get_admiral_option("subject_keys"),
    order = exprs(desc(AVALC), ADT),
    mode = "first",
    events = list(irsp_y, icb_y, no_data_n),
    source_datasets = list(
      ovr = ovr,
      adsl = adsl
    ),
    set_values_to = exprs(
      PARAMCD = "ICB",
      PARAM = "iRECIST Clinical Benefit by Investigator (confirmation for response not required)",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "iRECIST",
      AVAL = yn_to_numeric(AVALC),
      ANL01FL = "Y"
    ),
    check_type = "none"
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, RANDDT, ANL01FL),
  filter = PARAMCD == "ICB"
)

## -----------------------------------------------------------------------------
ibor_icr <- event(
  description = "Define complete response (iCR) for best overall response (iBOR)",
  dataset_name = "ovr",
  condition = AVALC == "iCR",
  set_values_to = exprs(AVALC = "iCR")
)

ibor_ipr <- event(
  description = "Define partial response (iPR) for best overall response (iBOR)",
  dataset_name = "ovr",
  condition = AVALC == "iPR",
  set_values_to = exprs(AVALC = "iPR")
)

ibor_isd <- event(
  description = paste(
    "Define stable disease (iSD) for best overall response (iBOR) as iCR, iPR, or iSD",
    "occurring at least 42 days after randomization"
  ),
  dataset_name = "ovr",
  condition = AVALC %in% c("iCR", "iPR", "iSD") & ADT >= RANDDT + 42,
  set_values_to = exprs(AVALC = "iSD")
)

ibor_non_icriupd <- event(
  description = paste(
    "Define NON-iCR/NON-iUPD for best overall response (iBOR) as NON-iCR/NON-iUPD",
    "occuring at least 42 days after randomization"
  ),
  dataset_name = "ovr",
  condition = AVALC == "NON-iCR/NON-iUPD" & ADT >= RANDDT + 42,
  set_values_to = exprs(AVALC = "NON-iCR/NON-iUPD")
)


ibor_icpd <- event_joined(
  description = paste(
    "Define confirmed progressive disease (iCPD) for best overall response (iBOR) as",
    "iUPD followed by iCPD with only other iUPD and NE responses in between"
  ),
  dataset_name = "ovr",
  join_vars = exprs(AVALC, ADT),
  join_type = "after",
  first_cond_upper = AVALC.join == "iCPD",
  condition = AVALC == "iUPD" &
    all(AVALC.join %in% c("iCPD", "iUPD", "NE")),
  set_values_to = exprs(AVALC = "iCPD")
)


ibor_iupd <- event(
  description = "Define unconfirmed progressive disease (iUPD) for best overall response (iBOR)",
  dataset_name = "ovr",
  condition = AVALC == "iUPD",
  set_values_to = exprs(AVALC = "iUPD")
)

ibor_ne <- event(
  description = paste(
    "Define not evaluable (NE) for best overall response (iBOR) as iCR, iPR, iSD,",
    "NON-iCR/NON-iUPD, or NE (should be specified after ibor_isd and ibor_non_icriupd)"
  ),
  dataset_name = "ovr",
  condition = AVALC %in% c("iCR", "iPR", "iSD", "NON-iCR/NON-iUPD", "NE"),
  set_values_to = exprs(AVALC = "NE")
)

adrs <- adrs %>%
  derive_extreme_event(
    by_vars = get_admiral_option("subject_keys"),
    tmp_event_nr_var = event_nr,
    order = exprs(event_nr, ADT),
    mode = "first",
    source_datasets = list(
      ovr = ovr,
      adsl = adsl
    ),
    events = list(ibor_icr, ibor_ipr, ibor_isd, ibor_non_icriupd, ibor_icpd, ibor_iupd, ibor_ne, no_data_missing),
    set_values_to = exprs(
      PARAMCD = "IBOR",
      PARAM = "iRECIST Best Overall Response by Investigator (confirmation not required)",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "iRECIST",
      AVAL = aval_resp_new(AVALC),
      ANL01FL = "Y"
    )
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, RANDDT, ANL01FL),
  filter = PARAMCD == "IBOR"
)

## -----------------------------------------------------------------------------
confirmation_period <- 28

icrsp_y_icr <- event_joined(
  description = paste(
    "Define confirmed response as iCR followed by iCR at least",
    confirmation_period,
    "days later and at most one NE in between"
  ),
  dataset_name = "ovr",
  join_vars = exprs(AVALC, ADT),
  join_type = "after",
  order = exprs(ADT),
  first_cond_upper = AVALC.join == "iCR" &
    ADT.join >= ADT + days(confirmation_period),
  condition = AVALC == "iCR" &
    all(AVALC.join %in% c("iCR", "NE")) &
    count_vals(var = AVALC.join, val = "NE") <= 1,
  set_values_to = exprs(AVALC = "Y")
)

icrsp_y_ipr <- event_joined(
  description = paste(
    "Define confirmed response as iPR followed by iCR or iPR at least",
    confirmation_period,
    "days later at most one NE in between, and no iPR after iCR"
  ),
  dataset_name = "ovr",
  join_vars = exprs(AVALC, ADT),
  join_type = "after",
  order = exprs(ADT),
  first_cond_upper = AVALC.join %in% c("iCR", "iPR") &
    ADT.join >= ADT + days(confirmation_period),
  condition = AVALC == "iPR" &
    all(AVALC.join %in% c("iCR", "iPR", "NE")) &
    count_vals(var = AVALC.join, val = "NE") <= 1 &
    (
      min_cond(
        var = ADT.join,
        cond = AVALC.join == "iCR"
      ) > max_cond(var = ADT.join, cond = AVALC.join == "iPR") |
        count_vals(var = AVALC.join, val = "iCR") == 0 |
        count_vals(var = AVALC.join, val = "iPR") == 0
    ),
  set_values_to = exprs(AVALC = "Y")
)

icbor_icr <- event_joined(
  description = paste(
    "Define complete response (iCR) for confirmed best overall response (iCBOR) as",
    "iCR followed by iCR at least",
    confirmation_period,
    "days later and at most one NE in between"
  ),
  dataset_name = "ovr",
  join_vars = exprs(AVALC, ADT),
  join_type = "after",
  first_cond_upper = AVALC.join == "iCR" &
    ADT.join >= ADT + confirmation_period,
  condition = AVALC == "iCR" &
    all(AVALC.join %in% c("iCR", "NE")) &
    count_vals(var = AVALC.join, val = "NE") <= 1,
  set_values_to = exprs(AVALC = "iCR")
)

icbor_ipr <- event_joined(
  description = paste(
    "Define partial response (iPR) for confirmed best overall response (iCBOR) as",
    "iPR followed by iCR or iPR at least",
    confirmation_period,
    "days later, at most one NE in between and no iPR after iCR"
  ),
  dataset_name = "ovr",
  join_vars = exprs(AVALC, ADT),
  join_type = "after",
  first_cond_upper = AVALC.join %in% c("iCR", "iPR") &
    ADT.join >= ADT + confirmation_period,
  condition = AVALC == "iPR" &
    all(AVALC.join %in% c("iCR", "iPR", "NE")) &
    count_vals(var = AVALC.join, val = "NE") <= 1 &
    (
      min_cond(
        var = ADT.join,
        cond = AVALC.join == "iCR"
      ) > max_cond(var = ADT.join, cond = AVALC.join == "iPR") |
        count_vals(var = AVALC.join, val = "iCR") == 0 |
        count_vals(var = AVALC.join, val = "iPR") == 0
    ),
  set_values_to = exprs(AVALC = "iPR")
)

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  derive_extreme_event(
    by_vars = get_admiral_option("subject_keys"),
    order = exprs(desc(AVALC), ADT),
    mode = "first",
    source_datasets = list(
      ovr = ovr,
      adsl = adsl
    ),
    events = list(icrsp_y_icr, icrsp_y_ipr, no_data_n),
    set_values_to = exprs(
      PARAMCD = "ICRSP",
      PARAM = "iRECIST Confirmed Response by Investigator",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "iRECIST",
      AVAL = yn_to_numeric(AVALC),
      ANL01FL = "Y"
    )
  )

adrs <- adrs %>%
  derive_extreme_event(
    by_vars = get_admiral_option("subject_keys"),
    order = exprs(desc(AVALC), ADT),
    mode = "first",
    events = list(icrsp_y_icr, icrsp_y_ipr, icb_y, no_data_n),
    source_datasets = list(
      ovr = ovr,
      adsl = adsl
    ),
    set_values_to = exprs(
      PARAMCD = "ICCB",
      PARAM = "iRECIST Confirmed Clinical Benefit by Investigator",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "iRECIST",
      AVAL = yn_to_numeric(AVALC),
      ANL01FL = "Y"
    ),
    check_type = "none"
  )

adrs <- adrs %>%
  derive_extreme_event(
    by_vars = get_admiral_option("subject_keys"),
    tmp_event_nr_var = event_nr,
    order = exprs(event_nr, ADT),
    mode = "first",
    events = list(icbor_icr, icbor_ipr, ibor_isd, ibor_non_icriupd, ibor_icpd, ibor_iupd, ibor_ne, no_data_missing),
    source_datasets = list(
      ovr = ovr,
      adsl = adsl
    ),
    set_values_to = exprs(
      PARAMCD = "ICBOR",
      PARAM = "iRECIST Best Confirmed Overall Response by Investigator",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "iRECIST",
      AVAL = aval_resp(AVALC),
      ANL01FL = "Y"
    )
  )

## ----eval=TRUE, echo=FALSE----------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, RANDDT, ANL01FL),
  filter = PARAMCD %in% c("ICRSP", "ICCB", "ICBOR")
)

