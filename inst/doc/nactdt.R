## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(admiralonco)
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

## ---- warning=FALSE, message=FALSE--------------------------------------------
library(admiral)
library(dplyr)

adsl <- admiral_adsl

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
  admiral_adsl,
  dataset_add = cm,
  by_vars = exprs(STUDYID, USUBJID),
  order = exprs(NACTDT),
  mode = "first",
  new_vars = exprs(NACTDT = convert_dtc_to_dt(CMSTDTC)),
  filter_add = CMSCAT == "CHEMOTHERAPY" & CMCAT == "ON TREATMENT"
)

## ---- eval=TRUE, echo=FALSE---------------------------------------------------
dataset_vignette(
  adsl,
  display_vars = exprs(USUBJID, NACTDT),
  filter = !is.na(NACTDT)
)

## ----message=FALSE------------------------------------------------------------
cm_date <- date_source(
  dataset_name = "cm",
  filter = CMSCAT == "CHEMOTHERAPY" & CMCAT == "ON TREATMENT" & !is.na(CMSTDTC),
  date = convert_dtc_to_dt(CMSTDTC)
)

pr_date <- date_source(
  dataset_name = "pr",
  filter = PRCAT == "CANCER RELATED" & PRSCAT == "ON TREATMENT" & !is.na(PRSTDTC),
  date = convert_dtc_to_dt(PRSTDTC)
)

## ----message=FALSE------------------------------------------------------------
adsl <- admiral_adsl %>%
  derive_var_extreme_dt(
    new_var = NACTDT,
    cm_date, pr_date,
    source_datasets = list(
      cm = cm,
      pr = pr
    ),
    mode = "first"
  )

## ---- eval=TRUE, echo=FALSE---------------------------------------------------
dataset_vignette(
  adsl,
  display_vars = exprs(USUBJID, NACTDT),
  filter = !is.na(NACTDT)
)

## ---- eval=TRUE, echo=TRUE----------------------------------------------------
library(admiralonco)

adrs <- derive_param_extreme_record(
  dataset = admiral_adrs,
  sources = list(
    records_source(
      dataset_name = "cm",
      filter = CMSCAT == "CHEMOTHERAPY" & CMCAT == "ON TREATMENT" & !is.na(CMSTDTC),
      new_vars = exprs(
        ADT = convert_dtc_to_dt(CMSTDTC),
        AVALC = CMTRT
      )
    ),
    records_source(
      dataset_name = "pr",
      filter = PRCAT == "CANCER RELATED" & PRSCAT == "ON TREATMENT" & !is.na(PRSTDTC),
      new_vars = exprs(
        ADT = convert_dtc_to_dt(PRSTDTC),
        AVALC = PRTRT
      )
    )
  ),
  source_datasets = list(cm = cm, pr = pr),
  by_vars = exprs(STUDYID, USUBJID),
  order = exprs(ADT),
  mode = "first",
  set_values_to = exprs(
    PARAMCD = "NACTDT",
    PARAM = "New Anti-Cancer Therapy Start Date"
  )
)

## ---- eval=TRUE, echo=FALSE---------------------------------------------------
dataset_vignette(
  adrs,
  display_vars = exprs(USUBJID, PARAMCD, PARAM, ADT, AVALC),
  filter = !is.na(ADT) & PARAMCD == "NACTDT"
)

