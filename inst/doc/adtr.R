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
admiral_homepage <- "https://pharmaverse.github.io/admiral"
library(admiraldev)

## ---- warning=FALSE, message=FALSE--------------------------------------------
library(admiral)
library(dplyr)
library(admiral.test)
library(lubridate)
library(stringr)
library(admiralonco)

## -----------------------------------------------------------------------------
data("admiral_adsl")
data("admiral_adrs")
data("admiral_rs")
data("admiral_tu")
data("admiral_tr")
adsl <- admiral_adsl
adrs <- admiral_adrs
tu <- admiral_tu
tr <- admiral_tr
rs <- admiral_rs

tu <- convert_blanks_to_na(tu) %>%
  filter(TUEVAL == "INVESTIGATOR")
tr <- convert_blanks_to_na(tr) %>%
  filter(
    TREVAL == "INVESTIGATOR" & TRGRPID == "TARGET" & TRTESTCD %in% c("LDIAM", "LPERP")
  )
rs <- convert_blanks_to_na(rs)

## ----echo=FALSE---------------------------------------------------------------
subjects <- c("01-701-1015", "01-701-1023", "01-701-1440")
# exclude assessments after PD
tr <- filter(tr, USUBJID %in% subjects & (USUBJID != "01-701-1440" | VISITNUM <= 7))
tu <- filter(tu, USUBJID %in% subjects)
rs <- filter(rs, USUBJID %in% subjects)
adsl <- filter(adsl, USUBJID %in% subjects)
adrs <- filter(adrs, USUBJID %in% subjects)

## ----eval=TRUE----------------------------------------------------------------
adsl_vars <- exprs(RANDDT)
tr <- derive_vars_merged(
  tr,
  dataset_add = adsl,
  new_vars = adsl_vars,
  by_vars = exprs(STUDYID, USUBJID)
)

## -----------------------------------------------------------------------------
tr <- derive_vars_merged(
  tr,
  dataset_add = tu,
  new_vars = exprs(TULOC),
  by_vars = exprs(STUDYID, USUBJID, TRLNKID = TULNKID)
) %>% mutate(
  TULOCGR1 = if_else(
    TULOC == "LYMPH NODE",
    "NODAL",
    "NON-NODAL"
  )
)

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(tr, display_vars = exprs(USUBJID, VISIT, TRLNKID, TULOC, TULOCGR1))

## -----------------------------------------------------------------------------
tr <- mutate(
  tr,
  LSEXP = TRLNKID,
  LSASS = if_else(!is.na(TRSTRESN), TRLNKID, NA_character_)
)

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(tr, display_vars = exprs(USUBJID, TRLNKID, VISIT, LSEXP, LSASS))

## -----------------------------------------------------------------------------
tr <- derive_vars_dt(
  tr,
  dtc = TRDTC,
  new_vars_prefix = "A",
  highest_imputation = "D",
  date_imputation = "first"
) %>%
  derive_vars_dy(
    reference_date = RANDDT,
    source_vars = exprs(ADT)
  )

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(tr, display_vars = exprs(USUBJID, RANDDT, TRLNKID, TRDTC, ADT, ADTF, ADY))

## -----------------------------------------------------------------------------
tr <- mutate(
  tr,
  AVISIT = VISIT,
  AVISITN = VISITNUM
)

## -----------------------------------------------------------------------------
tr <- mutate(tr, tmp_lesion_nr = str_sub(TRLNKID, 3))
adtr <- bind_rows(
  tr %>%
    filter(TRTESTCD == "LDIAM") %>%
    mutate(
      PARAMCD = paste0("LDIAM", tmp_lesion_nr),
      PARAM = paste("Target Lesion", tmp_lesion_nr, "Analysis Diameter")
    ),
  tr %>%
    filter(TRTESTCD == "LPERP") %>%
    mutate(
      PARAMCD = paste0("NLDIAM", tmp_lesion_nr),
      PARAM = paste("Target Lesion", tmp_lesion_nr, "Analysis Perpendicular")
    )
) %>%
  mutate(
    PARCAT1 = "Target Lesion(s)",
    PARCAT2 = "Investigator",
    PARCAT3 = "Recist 1.1",
    AVAL = TRSTRESN,
    ANL01FL = if_else(!is.na(AVAL), "Y", NA_character_)
  ) %>%
  select(-tmp_lesion_nr)

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(
  arrange(adtr, USUBJID, TRLNKID, AVISITN),
  display_vars = exprs(USUBJID, TRLNKID, AVISIT, PARAMCD, PARAM, AVAL, ANL01FL)
)

## -----------------------------------------------------------------------------
adtr_sum <- get_summary_records(
  adtr,
  by_vars = exprs(STUDYID, USUBJID, !!!adsl_vars, AVISIT, AVISITN),
  filter = (str_starts(PARAMCD, "LDIAM") & TULOCGR1 == "NON-NODAL") |
    (str_starts(PARAMCD, "NLDIAM") & TULOCGR1 == "NODAL"),
  analysis_var = AVAL,
  summary_fun = function(x) sum(x, na.rm = TRUE),
  set_values_to = exprs(
    PARAMCD = "SDIAM",
    PARAM = "Target Lesions Sum of Diameters by Investigator",
    PARCAT1 = "Target Lesion(s)",
    PARCAT2 = "Investigator",
    PARCAT3 = "Recist 1.1"
  )
)

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(
  adtr_sum %>%
    arrange(USUBJID, AVISITN) %>%
    select(USUBJID, PARAMCD, PARAM, AVISIT, AVAL, everything()),
  display_vars = exprs(USUBJID, PARAMCD, PARAM, AVISIT, AVAL)
)

## -----------------------------------------------------------------------------
adtr_sum <- adtr_sum %>%
  derive_var_merged_summary(
    dataset_add = adtr,
    by_vars = exprs(USUBJID),
    filter_add = AVISIT == "BASELINE" &
      ((str_starts(PARAMCD, "LDIAM") & TULOCGR1 == "NON-NODAL") |
        (str_starts(PARAMCD, "NLDIAM") & TULOCGR1 == "NODAL")),
    new_var = LSEXP,
    analysis_var = TRLNKID,
    summary_fun = function(x) paste(sort(x), collapse = ", ")
  ) %>%
  derive_var_merged_summary(
    dataset_add = adtr,
    by_vars = exprs(USUBJID, AVISIT),
    filter_add = ((str_starts(PARAMCD, "LDIAM") & TULOCGR1 == "NON-NODAL") |
      (str_starts(PARAMCD, "NLDIAM") & TULOCGR1 == "NODAL")) & ANL01FL == "Y",
    new_var = LSASS,
    analysis_var = TRLNKID,
    summary_fun = function(x) paste(sort(x), collapse = ", ")
  ) %>%
  mutate(
    ANL01FL = if_else(LSEXP == LSASS, "Y", NA_character_)
  )

## -----------------------------------------------------------------------------
adtr_sum <- adtr_sum %>%
  derive_var_merged_summary(
    dataset_add = adtr,
    by_vars = exprs(USUBJID, AVISIT),
    filter_add = (str_starts(PARAMCD, "LDIAM") & TULOCGR1 == "NON-NODAL") |
      (str_starts(PARAMCD, "NLDIAM") & TULOCGR1 == "NODAL"),
    new_var = ADY,
    analysis_var = ADY,
    summary_fun = function(x) min(x, na.rm = TRUE)
  ) %>%
  derive_var_merged_summary(
    dataset_add = adtr,
    by_vars = exprs(USUBJID, AVISIT),
    filter_add = (str_starts(PARAMCD, "LDIAM") & TULOCGR1 == "NON-NODAL") |
      (str_starts(PARAMCD, "NLDIAM") & TULOCGR1 == "NODAL"),
    new_var = ADT,
    analysis_var = ADT,
    summary_fun = function(x) min(x, na.rm = TRUE)
  )

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(
  arrange(adtr_sum, USUBJID, AVISITN),
  display_vars = exprs(USUBJID, PARAMCD, AVISIT, LSEXP, LSASS, ANL01FL, ADT, ADY)
)

## -----------------------------------------------------------------------------
adtr_sum <- adtr_sum %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = exprs(USUBJID),
      order = exprs(ADY),
      new_var = ABLFL,
      mode = "last"
    ),
    filter = ADY <= 1
  ) %>%
  derive_var_base(
    by_vars = exprs(USUBJID)
  )

## -----------------------------------------------------------------------------
adtr_sum <- adtr_sum %>%
  derive_vars_joined(
    dataset_add = adtr_sum,
    by_vars = exprs(USUBJID),
    order = exprs(AVAL),
    new_vars = exprs(NADIR = AVAL),
    join_vars = exprs(ADY),
    filter_add = ANL01FL == "Y",
    filter_join = ADY.join < ADY,
    mode = "first",
    check_type = "none"
  )

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(
  adtr_sum %>%
    arrange(USUBJID, AVISITN) %>%
    select(USUBJID, PARAMCD, AVISIT, ADY, AVAL, NADIR, everything()),
  display_vars = exprs(USUBJID, PARAMCD, AVISIT, ADY, AVAL, NADIR)
)

## -----------------------------------------------------------------------------
adtr_sum <- adtr_sum %>%
  derive_var_chg() %>%
  derive_var_pchg() %>%
  mutate(
    CHGNAD = AVAL - NADIR,
    PCHGNAD = if_else(NADIR == 0, NA_real_, 100 * CHGNAD / NADIR)
  )

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(
  arrange(adtr_sum, USUBJID, AVISITN),
  display_vars = exprs(USUBJID, AVISIT, AVAL, NADIR, BASE, CHG, PCHG, CHGNAD, PCHGNAD)
)

## ----eval=TRUE----------------------------------------------------------------
adtr_sum <- adtr_sum %>%
  derive_var_merged_exist_flag(
    dataset_add = adrs,
    filter_add = PARAMCD == "PD",
    by_vars = exprs(USUBJID, ADT),
    new_var = PDFL,
    condition = AVALC == "Y"
  )

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(
  arrange(adtr_sum, USUBJID, AVISITN),
  display_vars = exprs(USUBJID, AVISIT, AVAL, NADIR, CHGNAD, PCHGNAD, PDFL)
)

adtr_sum <- select(adtr_sum, -PDFL)

## ----eval=TRUE----------------------------------------------------------------
adtr_sum <- adtr_sum %>%
  derive_var_merged_exist_flag(
    dataset_add = derive_vars_dt(
      rs,
      dtc = RSDTC,
      new_vars_prefix = "A",
      highest_imputation = "D",
      flag_imputation = "none"
    ),
    filter_add = RSTESTCD == "OVRLRESP" & RSEVAL == "INVESTIGATOR",
    by_vars = exprs(USUBJID, ADT),
    new_var = PDFL,
    condition = RSSTRESC == "PD"
  )

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(
  arrange(adtr_sum, USUBJID, AVISITN),
  display_vars = exprs(USUBJID, AVISIT, AVAL, NADIR, CHGNAD, PCHGNAD, PDFL)
)

adtr_sum <- select(adtr_sum, -PDFL)

## -----------------------------------------------------------------------------
adtr_sum <- adtr_sum %>% mutate(
  CRFL = if_else(AVAL == 0 & ANL01FL == "Y", "Y", NA_character_),
  CRNFL = if_else(NADIR == 0, "Y", NA_character_),
  PDFL = if_else(is.na(CRFL) & CRNFL == "Y", "Y", NA_character_)
)

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(
  arrange(adtr_sum, USUBJID, AVISITN),
  display_vars = exprs(USUBJID, AVISIT, AVAL, NADIR, CHGNAD, PCHGNAD, PDFL)
)

adtr_sum <- select(adtr_sum, -PDFL, -CRFL, -CRNFL)

## -----------------------------------------------------------------------------
adtr_sum <- adtr_sum %>% mutate(
  CRFL = if_else(AVAL == 0 & ANL01FL == "Y", "Y", NA_character_),
  CRNFL = if_else(NADIR == 0, "Y", NA_character_)
)

## -----------------------------------------------------------------------------
adtr_sum <- adtr_sum %>%
  mutate(
    PDFL = if_else(
      PCHGNAD >= 20 & CHGNAD >= 5 | is.na(CRFL) & CRNFL == "Y",
      "Y",
      NA_character_
    )
  )

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(
  arrange(adtr_sum, USUBJID, AVISITN),
  display_vars = exprs(USUBJID, AVISIT, AVAL, NADIR, CHGNAD, PCHGNAD, PDFL)
)

## -----------------------------------------------------------------------------
adtr_sum <- adtr_sum %>%
  derive_var_relative_flag(
    by_vars = exprs(USUBJID),
    order = exprs(ADT),
    new_var = POSTRNDFL,
    condition = ADT > RANDDT,
    mode = "first",
    selection = "after",
    inclusive = TRUE,
    flag_no_ref_groups = FALSE
  ) %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = exprs(USUBJID),
      new_var = ANL02FL,
      order = exprs(PCHG),
      mode = "first",
      check_type = "none"
    ),
    filter = ANL01FL == "Y" & POSTRNDFL == "Y"
  ) %>%
  select(-POSTRNDFL)

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(
  adtr_sum %>%
    arrange(USUBJID, AVISITN) %>%
    select(USUBJID, AVISIT, AVAL, NADIR, PCHG, ANL01FL, ANL02FL, everything()),
  display_vars = exprs(USUBJID, AVISIT, AVAL, NADIR, PCHG, ANL01FL, ANL02FL)
)

## -----------------------------------------------------------------------------
adtr_sum <- adtr_sum %>%
  restrict_derivation(
    derivation = derive_var_relative_flag,
    args = params(
      by_vars = exprs(USUBJID),
      new_var = ANL03FL,
      condition = PDFL == "Y",
      order = exprs(ADY),
      mode = "first",
      selection = "before",
      inclusive = FALSE
    ),
    filter = ANL01FL == "Y" | PDFL == "Y"
  )

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(
  adtr_sum %>%
    arrange(USUBJID, AVISITN) %>%
    select(USUBJID, AVISIT, AVAL, NADIR, PDFL, ANL01FL, ANL02FL, ANL03FL, everything()),
  display_vars = exprs(USUBJID, AVISIT, AVAL, NADIR, PDFL, ANL01FL, ANL03FL)
)

## -----------------------------------------------------------------------------
adtr_sum <- adtr_sum %>%
  mutate(
    ANL04FL = if_else(ANL01FL == "Y" | PDFL == "Y", "Y", NA_character_)
  )

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(
  adtr_sum %>%
    arrange(USUBJID, AVISITN) %>%
    select(USUBJID, AVISIT, AVAL, NADIR, PDFL, ANL01FL, ANL02FL, ANL03FL, ANL04FL, everything()),
  display_vars = exprs(USUBJID, AVISIT, AVAL, NADIR, PDFL, ANL01FL, ANL04FL)
)

## -----------------------------------------------------------------------------
adtr <- bind_rows(adtr, adtr_sum)

## ----eval=TRUE----------------------------------------------------------------
adtr <- adtr %>%
  derive_var_obs_number(
    by_vars = exprs(STUDYID, USUBJID),
    order = exprs(PARAMCD, AVISITN, TRSEQ),
    check_type = "error"
  )

## ---- echo=FALSE--------------------------------------------------------------
dataset_vignette(
  arrange(adtr, USUBJID, ASEQ),
  display_vars = exprs(USUBJID, PARAMCD, AVISIT, ASEQ, AVAL)
)

## ----eval=TRUE----------------------------------------------------------------
adtr <- adtr %>%
  derive_vars_merged(
    dataset_add = select(adsl, !!!negate_vars(adsl_vars)),
    by_vars = exprs(STUDYID, USUBJID)
  )

## ---- eval=TRUE, echo=FALSE---------------------------------------------------
dataset_vignette(
  adtr,
  display_vars = exprs(USUBJID, RFSTDTC, RFENDTC, DTHDTC, DTHFL, AGE, AGEU)
)

