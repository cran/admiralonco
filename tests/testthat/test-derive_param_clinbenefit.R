adsl <- tibble::tribble(
  ~USUBJID, ~TRTSDT,      ~EOSDT,
  "01",     "2020-12-06", "2022-03-06",
  "02",     "2021-01-16", "2022-02-03",
  "03",     "2021-01-09", "2021-02-24",
  "04",     "2021-04-21", "2021-09-15",
  "05",     "2021-06-10", "2021-10-31",
  "06",     "2021-07-04", "2021-09-01",
  "07",     "2021-07-01", "2022-09-02"
) %>%
  mutate(
    STUDYID = "AB42",
    TRTSDT = lubridate::as_date(TRTSDT),
    EOSDT = lubridate::as_date(EOSDT),
  )

adrs <- tibble::tribble(
  ~USUBJID, ~PARAMCD, ~AVALC,          ~ADT,
  "01",     "RSP",    "Y",             "2021-04-08",
  "02",     "RSP",    "N",             "2021-05-07",
  "03",     "RSP",    "N",             NA,
  "04",     "RSP",    "N",             NA,
  "06",     "RSP",    "N",             NA,
  "07",     "RSP",    "N",             NA,
  "01",     "PD",     "N",             NA,
  "02",     "PD",     "Y",             "2021-05-07",
  "03",     "PD",     "N",             NA,
  "04",     "PD",     "N",             NA,
  "06",     "PD",     "Y",             "2021-08-20",
  "07",     "PD",     "N",             NA,
  "01",     "OVR",    "SD",            "2021-03-07",
  "01",     "OVR",    "PR",            "2021-04-08",
  "02",     "OVR",    "SD",            "2021-03-07",
  "02",     "OVR",    NA,              "2021-04-07",
  "02",     "OVR",    "PD",            "2021-05-07",
  "03",     "OVR",    "SD",            "2021-01-30",
  "04",     "OVR",    "NE",            "2021-05-21",
  "04",     "OVR",    "NA",            "2021-06-30",
  "04",     "OVR",    "NE",            "2021-07-24",
  "04",     "OVR",    "ND",            "2021-09-30",
  "06",     "OVR",    "PD",            "2021-08-20",
  "06",     "OVR",    "SD",            "2021-09-22",
  "07",     "OVR",    "NON-CR/NON-PD", "2021-12-05"
) %>%
  mutate(
    STUDYID = "AB42",
    ADT = lubridate::as_date(ADT),
    ANL01FL = "Y"
  ) %>%
  derive_vars_merged(
    dataset_add = adsl,
    by_vars = exprs(STUDYID, USUBJID),
    new_vars = exprs(TRTSDT)
  )

pd <- admiral::date_source(
  dataset_name = "adrs",
  date = ADT,
  filter = PARAMCD == "PD" & AVALC == "Y" & ANL01FL == "Y"
)

resp <- admiral::date_source(
  dataset_name = "adrs",
  date = ADT,
  filter = PARAMCD == "RSP" & AVALC == "Y" & ANL01FL == "Y"
)

## Test 1: ignore NON-CR/NON-PD ----
test_that("derive_param_clinbenefit Test 1: ignore NON-CR/NON-PD", {
  input_cbr <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVALC, ~AVAL, ~ADT,
    "01",     "CBR",    "Y",    1,     "2021-03-07",
    "02",     "CBR",    "Y",    1,     "2021-03-07",
    "03",     "CBR",    "N",    0,     NA,
    "04",     "CBR",    "N",    0,     NA,
    "05",     "CBR",    "N",    0,     NA,
    "06",     "CBR",    "N",    0,     NA,
    "07",     "CBR",    "N",    0,     NA
  ) %>%
    mutate(
      STUDYID = "AB42",
      ADT = lubridate::as_date(ADT),
      ANL01FL = "Y"
    ) %>%
    left_join(adsl, by = c("STUDYID", "USUBJID")) %>%
    select(-EOSDT)

  expected_output <- bind_rows(adrs, input_cbr)

  actual_output <- derive_param_clinbenefit(
    dataset = adrs,
    dataset_adsl = adsl,
    filter_source = PARAMCD == "OVR",
    source_resp = resp,
    source_pd = pd,
    source_datasets = list(adrs = adrs),
    reference_date = TRTSDT,
    ref_start_window = 28,
    clinben_vals = c("CR", "PR", "SD"),
    set_values_to = exprs(
      AVAL = yn_to_numeric(AVALC),
      PARAMCD = "CBR",
      ANL01FL = "Y"
    )
  )

  expect_dfs_equal(actual_output, expected_output,
    keys = c("USUBJID", "PARAMCD", "ADT")
  )
})

## Test 2: No source_pd ----
test_that("derive_param_clinbenefit Test 2: No source_pd", {
  input_cbr <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVALC, ~AVAL, ~ADT,
    "01",     "CBR",    "Y",    1,     "2021-03-07",
    "02",     "CBR",    "Y",    1,     "2021-03-07",
    "03",     "CBR",    "N",    0,     NA,
    "04",     "CBR",    "N",    0,     NA,
    "05",     "CBR",    "N",    0,     NA,
    "06",     "CBR",    "Y",    1,     "2021-09-22",
    "07",     "CBR",    "Y",    1,     "2021-12-05"
  ) %>%
    mutate(
      STUDYID = "AB42",
      ADT = lubridate::as_date(ADT),
      ANL01FL = "Y"
    ) %>%
    left_join(adsl, by = c("STUDYID", "USUBJID")) %>%
    select(-EOSDT)

  expected_output_no_source_pd <- bind_rows(adrs, input_cbr)

  actual_output_no_source_pd <- derive_param_clinbenefit(
    dataset = adrs,
    dataset_adsl = adsl,
    filter_source = PARAMCD == "OVR",
    source_resp = resp,
    source_pd = NULL,
    source_datasets = list(adrs = adrs),
    reference_date = TRTSDT,
    ref_start_window = 28,
    set_values_to = exprs(
      AVAL = yn_to_numeric(AVALC),
      PARAMCD = "CBR",
      ANL01FL = "Y"
    )
  )

  expect_dfs_equal(actual_output_no_source_pd,
    expected_output_no_source_pd,
    keys = c("USUBJID", "PARAMCD", "ADT")
  )
})

## Test 3: Deprecation warning for aval_fun ----
test_that("derive_param_clinbenefit Test 3: Deprecation warning for aval_fun", {
  input_cbr <- tibble::tribble(
    ~USUBJID, ~PARAMCD, ~AVALC, ~AVAL, ~ADT,
    "01",     "CBR",    "Y",    1,     "2021-03-07",
    "02",     "CBR",    "Y",    1,     "2021-03-07",
    "03",     "CBR",    "N",    0,     NA,
    "04",     "CBR",    "N",    0,     NA,
    "05",     "CBR",    "N",    0,     NA,
    "06",     "CBR",    "Y",    1,     "2021-09-22",
    "07",     "CBR",    "Y",    1,     "2021-12-05"
  ) %>%
    mutate(
      STUDYID = "AB42",
      ADT = lubridate::as_date(ADT),
      ANL01FL = "Y"
    ) %>%
    left_join(adsl, by = c("STUDYID", "USUBJID")) %>%
    select(-EOSDT)

  expected_output_no_source_pd <- bind_rows(adrs, input_cbr)

  expect_warning(
    actual_output_no_source_pd <- derive_param_clinbenefit(
      dataset = adrs,
      dataset_adsl = adsl,
      filter_source = PARAMCD == "OVR",
      source_resp = resp,
      source_pd = NULL,
      source_datasets = list(adrs = adrs),
      reference_date = TRTSDT,
      ref_start_window = 28,
      aval_fun = yn_to_numeric,
      set_values_to = exprs(
        PARAMCD = "CBR",
        ANL01FL = "Y"
      )
    ),
    class = "lifecycle_warning_deprecated"
  )

  expect_dfs_equal(actual_output_no_source_pd,
    expected_output_no_source_pd,
    keys = c("USUBJID", "PARAMCD", "ADT")
  )
})
