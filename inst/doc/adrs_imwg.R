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
library(metatools)
library(cli)
data("adsl")
# IMWG sdtm data
data("rs_onco_imwg")
data("supprs_onco_imwg")

rs <- rs_onco_imwg
supprs <- supprs_onco_imwg

rs <- combine_supp(rs, supprs)

rs <- convert_blanks_to_na(rs)

## ----eval=TRUE, echo=FALSE----------------------------------------------------
dataset_vignette(
  rs,
  display_vars = exprs(USUBJID, RSTESTCD, RSSTRESC, VISIT, RSDTC, RSDY, PDOFL, DTHPDFL, NACTDT, PDIFL)
)

## ----eval=TRUE----------------------------------------------------------------
adsl_vars <- exprs(RANDDT, TRTSDT)
adrs <- derive_vars_merged(
  rs,
  dataset_add = adsl,
  new_vars = adsl_vars,
  by_vars = get_admiral_option("subject_keys")
)

## ----eval=TRUE, echo=FALSE----------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, RSTESTCD, VISIT, RSDTC, RANDDT, TRTSDT)
)

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  filter(RSEVAL == "INVESTIGATOR" & RSTESTCD == "OVRLRESP") %>%
  mutate(
    PARAMCD = "OVR",
    PARAM = "Overall Response by Investigator",
    PARCAT1 = "Investigator",
    PARCAT2 = "IMWG"
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, VISIT, RSTESTCD, RSEVAL, PARAMCD, PARAM, PARCAT1, PARCAT2)
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
aval_resp_imwg <- function(arg) {
  case_match(
    arg,
    "NE" ~ 8,
    "sCR" ~ 7,
    "CR" ~ 6,
    "VGPR" ~ 5,
    "PR" ~ 4,
    "MR" ~ 3,
    "SD" ~ 2,
    "PD" ~ 1,
    NA ~ NA_real_
  )
}

adrs <- adrs %>%
  mutate(
    AVALC = RSSTRESC,
    AVAL = aval_resp_imwg(AVALC)
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, PARAMCD, AVISIT, ADT, AVAL, AVALC)
)

## -----------------------------------------------------------------------------
aval_resp_conf <- function(arg) {
  case_match(
    arg,
    "PD" ~ 8,
    "sCR" ~ 7,
    "CR" ~ 6,
    "VGPR" ~ 5,
    "PR" ~ 4,
    "MR" ~ 3,
    "SD" ~ 2,
    "NE" ~ 1,
    NA ~ 0
  )
}

## ----echo=FALSE---------------------------------------------------------------
list_resp <- tribble(
  ~"Response at 1st time point", ~"Response at 2nd time point",
  ~"Confirmed Response at 1st time point",
  "sCR", "sCR", "sCR",
  "CR", "sCR/CR", "max(CR, last Confirmed Response)",
  "VGPR", "sCR/CR/VGPR", "max(VGPR, last Confirmed Response)",
  "PR", "sCR/CR/VGPR/PR", "max(PR, last Confirmed Response)",
  "MR", "sCR/CR/VGPR/PR/MR", "max(MR, last Confirmed Response)",
  "sCR", "CR", "max(CR, last Confirmed Response)",
  "sCR/CR", "VGPR", "max(VGPR, last Confirmed Response)",
  "sCR/CR/VGPR", "PR", "max(PR, last Confirmed Response)",
  "sCR/CR/VGPR/PR", "MR", "max(MR, last Confirmed Response)",
  "sCR/CR/VGPR/PR/MR", "SD/PD/NE/NA", "max(SD, last Confirmed Response)",
  "SD", "any", "max(SD, last Confirmed Response)",
  "NE", "any", "max(NE, last Confirmed Response)",
  "PD reason imaging", "any", "PD",
  "PD reason serum/urine", "PD/death", "PD",
  "PD reason serum/urine", "sCR/CR/VGPR/PR/MR/SD/NE/NA", "max(NE, last Confirmed Response)",
)

knitr::kable(list_resp)

## ----echo=FALSE---------------------------------------------------------------
list_suppfl <- tribble(
  ~"Variable Name", ~"Variable Label",
  "PDOFL", "Progressive Disease: Other",
  "PDIFL", "Progressive Disease: Imaging",
  "DTHPDFL", "Death Due to Progressive Disease",
  "NACTDT", "New Anti-Cancer Therapy Date"
)

knitr::kable(list_suppfl)

## -----------------------------------------------------------------------------
confirmation_period <- 84

derive_confirmed_response <- function(datain) {
  data_adrs <- datain %>%
    arrange(USUBJID, ADTM) %>%
    filter(AVALC != "NE") %>%
    group_by(USUBJID) %>%
    mutate(
      AVALC.next = lead(AVALC),
      AVAL.next = lead(AVAL),
      ADT.next = lead(ADT)
    ) %>%
    ungroup() %>%
    mutate(AVALC.confirmed = case_when(
      # better response
      AVALC %in% c("sCR", "CR", "VGPR", "PR", "MR") &
        AVALC.next %in% c("sCR", "CR", "VGPR", "PR", "MR") &
        (is.na(NACTDT) | ADT.next <= NACTDT) &
        AVAL.next >= AVAL ~ AVALC,
      # worse response
      AVALC %in% c("sCR", "CR", "VGPR", "PR", "MR") &
        AVALC.next %in% c("sCR", "CR", "VGPR", "PR", "MR", "SD") &
        (is.na(NACTDT) | ADT.next <= NACTDT) &
        AVAL.next < AVAL ~ AVALC.next,
      # next assessment PD, NA or after subsequent therapy
      AVALC %in% c("sCR", "CR", "VGPR", "PR", "MR") &
        (AVALC.next == "PD" | is.na(AVALC.next) |
          !is.na(NACTDT) & ADT.next > NACTDT) ~ "SD",
      # no need to confirm SD
      AVALC %in% c("SD") ~ AVALC,
      # confirmed progression
      AVALC == "PD" &
        (PDIFL == "Y" | DTHPDFL == "Y" | PDOFL == "Y" & AVALC.next == "PD") ~
        AVALC,
      # unconfirmed progression
      AVALC == "PD" & is.na(DTHPDFL) & PDOFL == "Y" &
        (AVALC.next %in% c("sCR", "CR", "VGPR", "PR", "MR", "SD") |
          is.na(AVALC.next)) ~ "NE"
    ))

  data_adrs_check <- data_adrs %>%
    mutate(diff_days = as.numeric(difftime(ADT.next, ADT, units = "days"))) %>%
    filter(diff_days > confirmation_period) %>%
    mutate(warn = paste(
      "For USUBJID", USUBJID, "to confirm", AVISIT,
      "visit, a visit that took place", diff_days, "days later was used."
    )) %>%
    pull(warn)

  if (length(data_adrs_check) > 0) {
    cli_warn("{data_adrs_check}")
  }

  data_adrs_ne <- datain %>%
    filter(AVALC == "NE") %>%
    mutate(AVALC.confirmed = AVALC)

  data_adrs_all <- bind_rows(data_adrs, data_adrs_ne) %>%
    arrange(USUBJID, ADTM) %>%
    mutate(AVAL.confirmed = aval_resp_conf(AVALC.confirmed)) %>%
    group_by(USUBJID) %>%
    # best Confirmed Response so far
    mutate(AVAL.confirmed = cummax(AVAL.confirmed)) %>%
    ungroup(USUBJID)

  # char mapping to go back to AVALC values
  avalc_resp_conf <- function(arg) {
    case_match(
      arg,
      8 ~ "PD",
      7 ~ "sCR",
      6 ~ "CR",
      5 ~ "VGPR",
      4 ~ "PR",
      3 ~ "MR",
      2 ~ "SD",
      1 ~ "NE",
      NA_real_ ~ NA
    )
  }

  data_adrs_all <- data_adrs_all %>%
    mutate(AVALC.confirmed = avalc_resp_conf(AVAL.confirmed)) %>%
    select(
      -AVAL, -AVALC, -AVAL.next, -AVALC.next, -ADT.next,
      -AVAL.confirmed
    ) %>%
    rename(AVALC = AVALC.confirmed) %>%
    mutate(AVAL = aval_resp_imwg(AVALC)) %>%
    mutate(
      PARAMCD = "COVR",
      PARAM = "Confirmed Response at Time Point by Investigator",
      PARCAT1 = "Investigator",
      PARCAT2 = "IMWG"
    )
}
adrs_imwg <- derive_confirmed_response(adrs)

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs_imwg,
  display_vars = exprs(USUBJID, PARAMCD, AVISIT, ADT, AVALC)
)

## -----------------------------------------------------------------------------
filter_consecutive_vals <- function(dataset, by_vars, order, var, val, n) {
  var <- enexpr(var)
  var_join <- sym(paste0(as.character(var), ".join"))
  data <- derive_var_obs_number(dataset, by_vars = by_vars, order = order, new_var = temp_seq)

  left_join(
    data, select(data, !!!by_vars, temp_seq, !!var),
    by = admiraldev::vars2chr(by_vars),
    suffix = c("", ".join"),
    relationship = "many-to-many"
  ) %>%
    filter(temp_seq <= temp_seq.join & !!var == !!val) %>%
    filter_relative(
      by_vars = c(by_vars, expr(temp_seq)),
      order = exprs(temp_seq.join),
      condition = !!var_join != !!val,
      mode = "first",
      selection = "before",
      keep_no_ref_groups = TRUE,
      inclusive = FALSE
    ) %>%
    derive_var_merged_summary(
      dataset = .,
      dataset_add = .,
      by_vars = c(by_vars, expr(temp_seq)),
      new_vars = exprs(nr_vals = n())
    ) %>%
    filter(nr_vals >= n) %>%
    filter_extreme(
      by_vars = c(by_vars, expr(temp_seq)),
      order = exprs(temp_seq.join),
      mode = "first",
      check_type = "none"
    ) %>%
    select(-temp_seq, -temp_seq.join, -!!var_join, -nr_vals)
}

# filter on three or more NEs in a row
many_nes <- filter_consecutive_vals(
  adrs,
  by_vars = get_admiral_option("subject_keys"),
  order = exprs(ADTM),
  var = AVALC,
  val = "NE",
  n = 3
)

if (nrow(many_nes) > 0) {
  cli_warn("There are subjects with more than three NEs in a row.")
}

## -----------------------------------------------------------------------------
adrs_imwg <- adrs_imwg %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = exprs(!!!get_admiral_option("subject_keys"), ADT),
      order = exprs(AVAL, RSSEQ),
      new_var = ANL01FL,
      mode = "first"
    ),
    filter = !is.na(AVAL) & AVALC != "MISSING" & ADT >= RANDDT
  )

## ----eval=TRUE----------------------------------------------------------------
adrs_imwg <- adrs_imwg %>%
  mutate(
    ANL02FL = case_when(
      !is.na(AVAL) & ADT >= RANDDT & ADT < NACTDT ~ "Y",
      is.na(NACTDT) ~ "Y",
      TRUE ~ NA_character_
    )
  )

## -----------------------------------------------------------------------------
adrs_imwg <- adrs_imwg %>%
  derive_var_relative_flag(
    by_vars = get_admiral_option("subject_keys"),
    order = exprs(ADT, RSSEQ),
    new_var = ANL03FL,
    condition = AVALC == "PD",
    mode = "first",
    selection = "before",
    inclusive = TRUE
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs_imwg,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, AVALC, ADT, ANL01FL, ANL02FL, ANL03FL)
)

## -----------------------------------------------------------------------------
ovr <- filter(adrs_imwg, PARAMCD == "COVR" & ANL01FL == "Y" & ANL02FL == "Y" & ANL03FL == "Y")

adrs <- bind_rows(adrs, adrs_imwg)

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  ovr,
  display_vars = exprs(USUBJID, AVISIT, AVALC, ADT, RANDDT)
)

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  derive_extreme_records(
    dataset_ref = adsl,
    dataset_add = ovr,
    by_vars = get_admiral_option("subject_keys"),
    filter_add = PARAMCD == "COVR" & AVALC == "PD",
    order = exprs(ADT),
    mode = "first",
    exist_flag = AVALC,
    false_value = "N",
    set_values_to = exprs(
      PARAMCD = "PD",
      PARAM = "Disease Progression by Investigator",
      PARCAT1 = "Investigator",
      PARCAT2 = "IMWG",
      AVAL = yn_to_numeric(AVALC),
      ANL01FL = "Y",
      ANL02FL = "Y",
      ANL03FL = "Y"
    )
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, AVALC, ADT),
  filter = PARAMCD == "PD"
)

## -----------------------------------------------------------------------------
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
rsp_y_imwg <- event(
  description = "Define sCR, CR, VGPR or PR as response",
  dataset_name = "ovr",
  condition = AVALC %in% c("sCR", "CR", "VGPR", "PR"),
  set_values_to = exprs(AVALC = "Y")
)

adrs <- adrs %>%
  derive_extreme_event(
    by_vars = get_admiral_option("subject_keys"),
    order = exprs(ADT),
    mode = "first",
    events = list(rsp_y_imwg, no_data_n),
    source_datasets = list(
      ovr = ovr,
      adsl = adsl
    ),
    set_values_to = exprs(
      PARAMCD = "RSP",
      PARAM = "IMWG Response by Investigator",
      PARCAT1 = "Investigator",
      PARCAT2 = "IMWG",
      AVAL = yn_to_numeric(AVALC),
      ANL01FL = "Y",
      ANL02FL = "Y",
      ANL03FL = "Y"
    )
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, RSORRES, AVISIT, PARAMCD, AVALC, ADT),
  filter = PARAMCD == "RSP" & AVALC == "Y"
)

## -----------------------------------------------------------------------------
sustained_period <- 42

cb_y_imwg <- event(
  description = paste(
    "Define sCR, CR, VGPR, PR, MR or SD occuring at least",
    sustained_period,
    "days after randomization as clinical benefit"
  ),
  dataset_name = "ovr",
  condition = AVALC %in% c("sCR", "CR", "VGPR", "PR", "MR", "SD") &
    ADT >= RANDDT + days(sustained_period),
  set_values_to = exprs(AVALC = "Y")
)

## -----------------------------------------------------------------------------
adrs <- adrs %>%
  derive_extreme_event(
    by_vars = get_admiral_option("subject_keys"),
    order = exprs(desc(AVALC), ADT),
    mode = "first",
    events = list(rsp_y_imwg, cb_y_imwg, no_data_n),
    source_datasets = list(
      ovr = ovr,
      adsl = adsl
    ),
    set_values_to = exprs(
      PARAMCD = "CB",
      PARAM = "IMWG Clinical Benefit by Investigator",
      PARCAT1 = "Investigator",
      PARCAT2 = "IMWG",
      AVAL = yn_to_numeric(AVALC),
      ANL01FL = "Y",
      ANL02FL = "Y",
      ANL03FL = "Y"
    ),
    check_type = "none"
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, RSORRES, AVISIT, PARAMCD, AVALC, ADT),
  filter = PARAMCD == "CB" & AVALC == "Y"
)

## -----------------------------------------------------------------------------
cr_y_imwg <- event(
  description = "Define sCR or CR as response",
  dataset_name = "ovr",
  condition = AVALC %in% c("sCR", "CR"),
  set_values_to = exprs(AVALC = "Y")
)

adrs <- adrs %>%
  derive_extreme_event(
    by_vars = get_admiral_option("subject_keys"),
    order = exprs(desc(AVALC), ADT),
    mode = "first",
    events = list(cr_y_imwg, no_data_n),
    source_datasets = list(
      ovr = ovr,
      adsl = adsl
    ),
    set_values_to = exprs(
      PARAMCD = "CRRSP",
      PARAM = "IMWG Complete Response by Investigator",
      PARCAT1 = "Investigator",
      PARCAT2 = "IMWG",
      AVAL = yn_to_numeric(AVALC),
      ANL01FL = "Y",
      ANL02FL = "Y",
      ANL03FL = "Y"
    ),
    check_type = "none"
  )

vgpr_y_imwg <- event(
  description = "Define sCR, CR or VGPR as response",
  dataset_name = "ovr",
  condition = AVALC %in% c("sCR", "CR", "VGPR"),
  set_values_to = exprs(AVALC = "Y")
)

adrs <- adrs %>%
  derive_extreme_event(
    by_vars = get_admiral_option("subject_keys"),
    order = exprs(desc(AVALC), ADT),
    mode = "first",
    events = list(vgpr_y_imwg, no_data_n),
    source_datasets = list(
      ovr = ovr,
      adsl = adsl
    ),
    set_values_to = exprs(
      PARAMCD = "VGPRRSP",
      PARAM = "IMWG VGPR Response by Investigator",
      PARCAT1 = "Investigator",
      PARCAT2 = "IMWG",
      AVAL = yn_to_numeric(AVALC),
      ANL01FL = "Y",
      ANL02FL = "Y",
      ANL03FL = "Y"
    ),
    check_type = "none"
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, RSORRES, AVISIT, PARAMCD, AVALC, ADT),
  filter = PARAMCD %in% c("CRRSP", "VGPRRSP") & AVALC == "Y"
)

## -----------------------------------------------------------------------------
bor_scr <- event(
  description = "Define stringent complete response (sCR) for best overall response
  (BOR)",
  dataset_name = "ovr",
  condition = AVALC == "sCR",
  set_values_to = exprs(AVALC = "sCR")
)

bor_vgpr <- event(
  description = "Define very good partial response (VGPR) for best overall response
  (BOR)",
  dataset_name = "ovr",
  condition = AVALC == "VGPR",
  set_values_to = exprs(AVALC = "VGPR")
)

bor_mr <- event(
  description = "Define minimal response (MR) for best overall response (BOR)",
  dataset_name = "ovr",
  condition = AVALC == "MR",
  set_values_to = exprs(AVALC = "MR")
)

bor_sd_imwg <- event(
  description = "Define stable disease (SD) for best overall response (BOR)",
  dataset_name = "ovr",
  condition = AVALC == "SD",
  set_values_to = exprs(AVALC = "SD")
)

bor_ne_imwg <- event(
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
      ovr = ovr,
      adsl = adsl
    ),
    events = list(
      bor_scr, bor_cr, bor_vgpr, bor_pr, bor_mr, bor_sd_imwg, bor_pd, bor_ne_imwg,
      no_data_missing
    ),
    set_values_to = exprs(
      PARAMCD = "CBOR",
      PARAM = "IMWG Best Confirmed Overall Response by Investigator",
      PARCAT1 = "Investigator",
      PARCAT2 = "IMWG",
      AVAL = aval_resp_imwg(AVALC),
      ANL01FL = "Y",
      ANL02FL = "Y",
      ANL03FL = "Y"
    )
  )

## ----echo=FALSE---------------------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, AVISIT, PARAMCD, AVALC, ADT),
  filter = PARAMCD == "CBOR" & AVALC != "MISSING"
)

