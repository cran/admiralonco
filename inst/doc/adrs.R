## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(admiral)

link <- function(text, url) {
  return(
    paste0(
      "[", text, "]",
      "(", url, ")"
    )
  )
}
dyn_link <- function(text,
                     base_url,
                     relative_url = "",
                     # Change to TRUE when admiral adopts multiversion docs
                     is_multiversion = FALSE,
                     multiversion_default_ref = "main") {
  url <- paste(base_url, relative_url, sep = "/")
  if (is_multiversion) {
    url <- paste(
      base_url,
      Sys.getenv("BRANCH_NAME", multiversion_default_ref),
      relative_url,
      sep = "/"
    )
  }
  return(link(text, url))
}
# Other variables
admiral_homepage <- "https://pharmaverse.github.io/admiral/cran-release"

library(admiraldev)

## ----message=FALSE------------------------------------------------------------
library(admiral)
library(admiralonco)
library(dplyr)
library(admiral.test)
library(lubridate)
library(stringr)
data("admiral_adsl")
data("admiral_rs")
data("admiral_tu")

adsl <- admiral_adsl
rs <- admiral_rs
tu <- admiral_tu

rs <- convert_blanks_to_na(rs)
tu <- convert_blanks_to_na(tu)

## ----echo=FALSE---------------------------------------------------------------
rs <- filter(rs, USUBJID %in% c("01-701-1015", "01-701-1023", "01-703-1086", "01-703-1096", "01-707-1037", "01-716-1024")) %>%
  ungroup()
tu <- filter(tu, USUBJID %in% c("01-701-1015", "01-701-1023", "01-703-1086", "01-703-1096", "01-707-1037", "01-716-1024")) %>%
  ungroup()
adsl <- filter(adsl, USUBJID %in% c("01-701-1015", "01-701-1023", "01-703-1086", "01-703-1096", "01-707-1037", "01-716-1024"))

## ----eval=TRUE----------------------------------------------------------------
adsl_vars <- exprs(RANDDT)
adrs <- derive_vars_merged(
  rs,
  dataset_add = adsl,
  new_vars = adsl_vars,
  by_vars = exprs(STUDYID, USUBJID)
)

## ---- eval=TRUE, echo=FALSE---------------------------------------------------
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
    PARCAT3 = "Recist 1.1"
  )

## ---- echo=FALSE--------------------------------------------------------------
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

## ---- echo=FALSE--------------------------------------------------------------
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

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, RSSTRESC, AVALC, AVAL)
)

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = exprs(STUDYID, USUBJID, ADT),
      order = exprs(AVAL, RSSEQ),
      new_var = ANL01FL,
      mode = "last"
    ),
    filter = !is.na(AVAL) & ADT >= RANDDT
  )

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, RANDDT, ANL01FL)
)

## ---- eval=FALSE--------------------------------------------------------------
#  adrs <- adrs %>%
#    mutate(
#      ANL01FL = case_when(
#        !is.na(AVAL) & ADT >= RANDDT & ADT < NACTDT ~ "Y",
#        TRUE ~ NA_character_
#      )
#    )

## ---- eval=FALSE--------------------------------------------------------------
#  adrs <- adrs %>%
#    derive_var_relative_flag(
#      by_vars = exprs(USUBJID),
#      order = exprs(ADT, AVISITN),
#      new_var = ANL02FL,
#      condition = AVALC == "PD",
#      mode = "first",
#      selection = "before",
#      inclusive = TRUE
#    )

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  derive_param_extreme_event(
    dataset_adsl = adsl,
    dataset_source = adrs,
    filter_source = PARAMCD == "OVR" & AVALC == "PD" & ANL01FL == "Y",
    order = exprs(ADT, RSSEQ),
    set_values_to = exprs(
      PARAMCD = "PD",
      PARAM = "Disease Progression by Investigator",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "Recist 1.1",
      ANL01FL = "Y"
    )
  )

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, ANL01FL),
  filter = PARAMCD == "PD"
)

## -----------------------------------------------------------------------------
pd <- date_source(
  dataset_name = "adrs",
  date = ADT,
  filter = PARAMCD == "PD" & AVALC == "Y"
)

## ---- eval=FALSE--------------------------------------------------------------
#  pd <- date_source(
#    dataset_name = "adsl",
#    date = PDDT
#  )

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  derive_param_response(
    dataset_adsl = adsl,
    filter_source = PARAMCD == "OVR" & AVALC %in% c("CR", "PR") & ANL01FL == "Y",
    source_pd = pd,
    source_datasets = list(adrs = adrs),
    set_values_to = exprs(
      PARAMCD = "RSP",
      PARAM = "Response by Investigator (confirmation not required)",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "Recist 1.1",
      ANL01FL = "Y"
    )
  )

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, ANL01FL),
  filter = PARAMCD == "RSP"
)

## -----------------------------------------------------------------------------
resp <- date_source(
  dataset_name = "adrs",
  date = ADT,
  filter = PARAMCD == "RSP" & AVALC == "Y"
)

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  derive_param_clinbenefit(
    dataset_adsl = adsl,
    filter_source = PARAMCD == "OVR" & ANL01FL == "Y",
    source_resp = resp,
    source_pd = pd,
    source_datasets = list(adrs = adrs),
    reference_date = RANDDT,
    ref_start_window = 42,
    set_values_to = exprs(
      PARAMCD = "CB",
      PARAM = "Clinical Benefit by Investigator (confirmation for response not required)",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "Recist 1.1",
      ANL01FL = "Y"
    )
  )

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, RANDDT, ANL01FL),
  filter = PARAMCD == "CB"
)

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  derive_param_bor(
    dataset_adsl = adsl,
    filter_source = PARAMCD == "OVR" & ANL01FL == "Y",
    source_pd = pd,
    source_datasets = list(adrs = adrs),
    reference_date = RANDDT,
    ref_start_window = 42,
    set_values_to = exprs(
      PARAMCD = "BOR",
      PARAM = "Best Overall Response by Investigator (confirmation not required)",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "Recist 1.1",
      ANL01FL = "Y"
    )
  )

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, RANDDT, ANL01FL),
  filter = PARAMCD == "BOR"
)

## ---- eval=FALSE--------------------------------------------------------------
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
  derive_param_extreme_event(
    dataset_adsl = adsl,
    dataset_source = adrs,
    filter_source = PARAMCD == "BOR" & AVALC %in% c("CR", "PR") & ANL01FL == "Y",
    order = exprs(ADT, RSSEQ),
    set_values_to = exprs(
      PARAMCD = "BCP",
      PARAM = "Best Overall Response of CR/PR by Investigator (confirmation not required)",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "Recist 1.1",
      ANL01FL = "Y"
    )
  )

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, ANL01FL),
  filter = PARAMCD == "BCP"
)

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  derive_param_confirmed_resp(
    dataset_adsl = adsl,
    filter_source = PARAMCD == "OVR" & AVALC %in% c("CR", "PR") & ANL01FL == "Y",
    source_pd = pd,
    source_datasets = list(adrs = adrs),
    ref_confirm = 28,
    set_values_to = exprs(
      PARAMCD = "CRSP",
      PARAM = "Confirmed Response by Investigator",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "Recist 1.1",
      ANL01FL = "Y"
    )
  )

confirmed_resp <- date_source(
  dataset_name = "adrs",
  date = ADT,
  filter = PARAMCD == "CRSP" & AVALC == "Y"
)

adrs <- adrs %>%
  derive_param_clinbenefit(
    dataset_adsl = adsl,
    filter_source = PARAMCD == "OVR" & ANL01FL == "Y",
    source_resp = confirmed_resp,
    source_pd = pd,
    source_datasets = list(adrs = adrs),
    reference_date = RANDDT,
    ref_start_window = 42,
    set_values_to = exprs(
      PARAMCD = "CCB",
      PARAM = "Confirmed Clinical Benefit by Investigator",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "Recist 1.1",
      ANL01FL = "Y"
    )
  ) %>%
  derive_param_confirmed_bor(
    dataset_adsl = adsl,
    filter_source = PARAMCD == "OVR" & ANL01FL == "Y",
    source_pd = pd,
    source_datasets = list(adrs = adrs),
    reference_date = RANDDT,
    ref_start_window = 42,
    ref_confirm = 28,
    set_values_to = exprs(
      PARAMCD = "CBOR",
      PARAM = "Best Confirmed Overall Response by Investigator",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "Recist 1.1",
      ANL01FL = "Y"
    )
  ) %>%
  derive_param_extreme_event(
    dataset_adsl = adsl,
    dataset_source = adrs,
    filter_source = PARAMCD == "CBOR" & AVALC %in% c("CR", "PR") & ANL01FL == "Y",
    order = exprs(ADT, RSSEQ),
    set_values_to = exprs(
      PARAMCD = "CBCP",
      PARAM = "Best Confirmed Overall Response of CR/PR by Investigator",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "Recist 1.1",
      ANL01FL = "Y"
    )
  )

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, RANDDT, ANL01FL),
  filter = PARAMCD %in% c("CRSP", "CCB", "CBOR", "CBCP")
)

## -----------------------------------------------------------------------------
adrsirf <- rs %>%
  filter(RSEVAL == "INDEPENDENT ASSESSOR" & RSEVALID == "RADIOLOGIST 1" & RSTESTCD == "OVRLRESP") %>%
  mutate(
    PARAMCD = "OVRR1",
    PARAM = "Overall Response by Radiologist 1",
    PARCAT1 = "Tumor Response",
    PARCAT2 = "Radiologist",
    PARCAT3 = "Recist 1.1"
  )

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(
  adrsirf,
  display_vars = exprs(USUBJID, VISIT, RSTESTCD, RSEVAL, PARAMCD, PARAM, PARCAT1, PARCAT2, PARCAT3),
  filter = PARAMCD == "OVRR1"
)

## -----------------------------------------------------------------------------
adsldth <- adsl %>%
  select(STUDYID, USUBJID, DTHDT, !!!adsl_vars)

adrs <- adrs %>%
  derive_param_extreme_event(
    dataset_adsl = adsldth,
    dataset_source = adsldth,
    filter_source = !is.na(DTHDT),
    set_values_to = exprs(
      PARAMCD = "DEATH",
      PARAM = "Death",
      PARCAT1 = "Reference Event",
      ANL01FL = "Y",
      ADT = DTHDT
    )
  ) %>%
  select(-DTHDT)

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, ANL01FL),
  filter = PARAMCD == "DEATH"
)

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  derive_param_extreme_event(
    dataset_adsl = adsl,
    dataset_source = adrs,
    filter_source = PARAMCD == "OVR" & ANL01FL == "Y",
    order = exprs(ADT, RSSEQ),
    mode = "last",
    new_var = dummy,
    set_values_to = exprs(
      PARAMCD = "LSTA",
      PARAM = "Last Disease Assessment by Investigator",
      PARCAT1 = "Tumor Response",
      PARCAT2 = "Investigator",
      PARCAT3 = "Recist 1.1",
      ANL01FL = "Y"
    )
  ) %>%
  select(-dummy)

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, ANL01FL),
  filter = PARAMCD == "LSTA"
)

## -----------------------------------------------------------------------------
adslmdis <- adsl %>%
  select(STUDYID, USUBJID, !!!adsl_vars)

adrs <- adrs %>%
  derive_param_exist_flag(
    dataset_adsl = adslmdis,
    dataset_add = tu,
    condition = TUEVAL == "INVESTIGATOR" & TUSTRESC == "TARGET" & VISIT == "BASELINE",
    false_value = "N",
    missing_value = "N",
    set_values_to = exprs(
      PARAMCD = "MDIS",
      PARAM = "Measurable Disease at Baseline by Investigator",
      PARCAT2 = "Investigator",
      PARCAT3 = "Recist 1.1",
      ANL01FL = "Y"
    )
  )

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, PARAM, AVALC, ADT, ANL01FL),
  filter = PARAMCD == "MDIS"
)

## ----eval=TRUE----------------------------------------------------------------
adrs <- adrs %>%
  mutate(
    AVAL = case_when(
      AVALC == "Y" ~ 1,
      AVALC == "N" ~ 0,
      TRUE ~ AVAL
    )
  )

## ---- eval=TRUE, echo=FALSE---------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, PARAMCD, AVALC, AVAL),
  filter = USUBJID == "01-701-1015" & AVALC %in% c("Y", "N")
)

## ----eval=TRUE----------------------------------------------------------------
adrs <- adrs %>%
  derive_var_obs_number(
    by_vars = exprs(STUDYID, USUBJID),
    order = exprs(PARAMCD, ADT, VISITNUM, RSSEQ),
    check_type = "error"
  )

## ---- eval=TRUE, echo=FALSE---------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, PARAMCD, ADT, VISITNUM, AVISIT, ASEQ),
  filter = USUBJID == "01-701-1015"
)

## ----eval=TRUE----------------------------------------------------------------
adrs <- adrs %>%
  derive_vars_merged(
    dataset_add = select(adsl, !!!negate_vars(adsl_vars)),
    by_vars = exprs(STUDYID, USUBJID)
  )

## ---- eval=TRUE, echo=FALSE---------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, RFSTDTC, RFENDTC, DTHDTC, DTHFL, AGE, AGEU),
  filter = USUBJID == "01-701-1015"
)

