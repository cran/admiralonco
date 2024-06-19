## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(admiraldev)

## ----warning=FALSE, message=FALSE---------------------------------------------
library(admiral)
library(pharmaverseadam)
library(dplyr)

data("adsl")
adsl_onco <- adsl
data("adrs_onco")

cm <- tribble(
  ~STUDYID, ~USUBJID, ~CMCAT, ~CMSCAT, ~CMTRT, ~CMSTDTC,
  "CDISCPILOT01", "01-701-1015", "PRIOR TREATMENT", "CHEMOTHERAPY", "DEXRAZOXANE", NA,
  "CDISCPILOT01", "01-701-1015", "ON TREATMENT", "CHEMOTHERAPY", "DEXROZOXANE", "2014-07-02",
  "CDISCPILOT01", "01-701-1015", "ON TREATMENT", "CHEMOTHERAPY", "DEXROZOXANE", "2014-06-19",
  "CDISCPILOT01", "01-701-1028", "PRIOR TREATMENT", "CHEMOTHERAPY", "METHOTREXATE", NA,
  "CDISCPILOT01", "01-701-1028", "ON TREATMENT", "CHEMOTHERAPY", "METHOTREXATE", "2014-01-14",
  "CDISCPILOT01", "01-701-1034", "PRIOR TREATMENT", "CHEMOTHERAPY", "OLAPARIB", NA,
  "CDISCPILOT01", "01-701-1034", "ON TREATMENT", "CHEMOTHERAPY", "OLAPARIB", "2014-12-30",
  "CDISCPILOT01", "01-701-1097", "PRIOR TREATMENT", "CHEMOTHERAPY", "TEMODAL", NA,
  "CDISCPILOT01", "01-701-1097", "ON TREATMENT", "CHEMOTHERAPY", "TEMODAL", "2013-12-31",
)

pr <- tribble(
  ~STUDYID, ~USUBJID, ~PRCAT, ~PRSCAT, ~PRTRT, ~PRSTDTC,
  "CDISCPILOT01", "01-701-1015", "CANCER RELATED", "ON TREATMENT", "SURGERY", "2014-06-18",
  "CDISCPILOT01", "01-701-1034", "CANCER RELATED", "ON TREATMENT", "SURGERY", "2014-12-16",
  "CDISCPILOT01", "01-701-1028", "CANCER RELATED", "PRIOR TREATMENT", "SURGERY", NA,
)

## ----message=FALSE------------------------------------------------------------
adsl <- derive_vars_merged(
  adsl_onco,
  dataset_add = cm,
  by_vars = get_admiral_option("subject_keys"),
  order = exprs(NACTDT),
  mode = "first",
  new_vars = exprs(NACTDT = convert_dtc_to_dt(CMSTDTC)),
  filter_add = CMSCAT == "CHEMOTHERAPY" & CMCAT == "ON TREATMENT"
)

## ----eval=TRUE, echo=FALSE----------------------------------------------------
dataset_vignette(
  adsl,
  display_vars = exprs(USUBJID, NACTDT),
  filter = !is.na(NACTDT)
)

## ----message=FALSE------------------------------------------------------------
cm_date <- event(
  dataset_name = "cm",
  condition = CMSCAT == "CHEMOTHERAPY" & CMCAT == "ON TREATMENT" & !is.na(CMSTDTC),
  set_values_to = exprs(NACTDT = convert_dtc_to_dt(CMSTDTC))
)

pr_date <- event(
  dataset_name = "pr",
  condition = PRCAT == "CANCER RELATED" & PRSCAT == "ON TREATMENT" & !is.na(PRSTDTC),
  set_values_to = exprs(NACTDT = convert_dtc_to_dt(PRSTDTC))
)

## ----message=FALSE------------------------------------------------------------
adsl <- adsl_onco %>%
  derive_vars_extreme_event(
    by_vars = get_admiral_option("subject_keys"),
    order = exprs(NACTDT),
    new_vars = exprs(NACTDT),
    events = list(cm_date, pr_date),
    source_datasets = list(
      cm = cm,
      pr = pr
    ),
    mode = "first"
  )

## ----eval=TRUE, echo=FALSE----------------------------------------------------
dataset_vignette(
  adsl,
  display_vars = exprs(USUBJID, NACTDT),
  filter = !is.na(NACTDT)
)

## ----eval=TRUE, echo=TRUE-----------------------------------------------------
adrs <- derive_extreme_event(
  dataset = adrs_onco,
  events = list(
    event(
      dataset_name = "cm",
      condition = CMSCAT == "CHEMOTHERAPY" & CMCAT == "ON TREATMENT" & !is.na(CMSTDTC),
      set_values_to = exprs(
        ADT = convert_dtc_to_dt(CMSTDTC),
        AVALC = CMTRT
      )
    ),
    event(
      dataset_name = "pr",
      condition = PRCAT == "CANCER RELATED" & PRSCAT == "ON TREATMENT" & !is.na(PRSTDTC),
      set_values_to = exprs(
        ADT = convert_dtc_to_dt(PRSTDTC),
        AVALC = PRTRT
      )
    )
  ),
  source_datasets = list(cm = cm, pr = pr),
  by_vars = get_admiral_option("subject_keys"),
  order = exprs(ADT),
  mode = "first",
  set_values_to = exprs(
    PARAMCD = "NACTDT",
    PARAM = "New Anti-Cancer Therapy Start Date"
  )
)

## ----eval=TRUE, echo=FALSE----------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, PARAMCD, PARAM, ADT, AVALC),
  filter = !is.na(ADT) & PARAMCD == "NACTDT"
)

