# filter_pd ----
## filter_pd Test 1: first PD in separate BDS dataset ----
test_that("filter_pd Test 1: first PD in separate BDS dataset", {
  adrs <- tibble::tribble(
    ~STUDYID,       ~USUBJID,      ~PARAMCD, ~AVALC, ~ADT,
    "CDISCPILOT01", "01-701-1015", "OVR",    "CR",   "2016-01-25",
    "CDISCPILOT01", "01-701-1015", "OVR",    "PD",   "2016-02-22",
    "CDISCPILOT01", "01-701-1015", "BOR",    "CR",   "2016-01-25",
    "CDISCPILOT01", "01-701-1034", "OVR",    "SD",   "2015-12-07",
    "CDISCPILOT01", "01-701-1034", "OVR",    "PD",   "2016-04-25",
    "CDISCPILOT01", "01-701-1034", "OVR",    "PD",   "2016-06-25",
    "CDISCPILOT01", "01-701-1034", "BOR",    "SD",   "2015-12-07",
    "CDISCPILOT01", "01-701-1035", "OVR",    "SD",   "2016-04-25",
    "CDISCPILOT01", "01-701-1035", "OVR",    "PR",   "2016-06-25",
    "CDISCPILOT01", "01-701-1035", "BOR",    "PR",   "2016-06-25"
  ) %>% mutate(
    ADT = lubridate::as_date(ADT)
  )

  adevent <- tibble::tribble(
    ~STUDYID,       ~USUBJID,      ~PARAMCD, ~AVALC, ~ADT,
    "CDISCPILOT01", "01-701-1015", "PD",     "Y",    "2016-02-22",
    "CDISCPILOT01", "01-701-1034", "PD",     "Y",    "2016-04-25"
  ) %>% mutate(
    ADT = lubridate::as_date(ADT)
  )

  expected_output <- tibble::tribble(
    ~STUDYID,       ~USUBJID,      ~PARAMCD, ~AVALC, ~ADT,
    "CDISCPILOT01", "01-701-1015", "OVR",    "CR",   "2016-01-25",
    "CDISCPILOT01", "01-701-1015", "OVR",    "PD",   "2016-02-22",
    "CDISCPILOT01", "01-701-1034", "OVR",    "SD",   "2015-12-07",
    "CDISCPILOT01", "01-701-1034", "OVR",    "PD",   "2016-04-25",
    "CDISCPILOT01", "01-701-1035", "OVR",    "SD",   "2016-04-25",
    "CDISCPILOT01", "01-701-1035", "OVR",    "PR",   "2016-06-25",
  ) %>% mutate(
    ADT = lubridate::as_date(ADT)
  )

  actual_output <- filter_pd(
    dataset = adrs,
    filter = PARAMCD == "OVR",
    source_pd = admiral::date_source(
      dataset_name = "adevent",
      date = ADT,
      filter = PARAMCD == "PD",
    ),
    source_datasets = list(adevent = adevent)
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID", "ADT")
  )
})

## filter_pd Test 2: first PD in ADSL dataset ----
test_that("filter_pd Test 2: first PD in ADSL dataset", {
  adrs <- tibble::tribble(
    ~STUDYID,       ~USUBJID,      ~PARAMCD, ~AVALC, ~ADT,
    "CDISCPILOT01", "01-701-1015", "OVR",    "CR",   "2016-01-25",
    "CDISCPILOT01", "01-701-1015", "OVR",    "PD",   "2016-02-22",
    "CDISCPILOT01", "01-701-1015", "BOR",    "CR",   "2016-01-25",
    "CDISCPILOT01", "01-701-1034", "OVR",    "SD",   "2015-12-07",
    "CDISCPILOT01", "01-701-1034", "OVR",    "PD",   "2016-04-25",
    "CDISCPILOT01", "01-701-1034", "OVR",    "PD",   "2016-06-25",
    "CDISCPILOT01", "01-701-1034", "BOR",    "SD",   "2015-12-07",
    "CDISCPILOT01", "01-701-1035", "OVR",    "SD",   "2016-04-25",
    "CDISCPILOT01", "01-701-1035", "OVR",    "PR",   "2016-06-25",
    "CDISCPILOT01", "01-701-1035", "BOR",    "PR",   "2016-06-25"
  ) %>% mutate(
    ADT = lubridate::as_date(ADT)
  )

  adsl <- tibble::tribble(
    ~STUDYID,       ~USUBJID,      ~PDDT,
    "CDISCPILOT01", "01-701-1015", "2016-02-22",
    "CDISCPILOT01", "01-701-1034", "2016-04-25"
  ) %>% mutate(
    PDDT = lubridate::as_date(PDDT)
  )

  expected_output <- tibble::tribble(
    ~STUDYID,       ~USUBJID,      ~PARAMCD, ~AVALC, ~ADT,
    "CDISCPILOT01", "01-701-1015", "OVR",    "CR",   "2016-01-25",
    "CDISCPILOT01", "01-701-1015", "OVR",    "PD",   "2016-02-22",
    "CDISCPILOT01", "01-701-1034", "OVR",    "SD",   "2015-12-07",
    "CDISCPILOT01", "01-701-1034", "OVR",    "PD",   "2016-04-25",
    "CDISCPILOT01", "01-701-1035", "OVR",    "SD",   "2016-04-25",
    "CDISCPILOT01", "01-701-1035", "OVR",    "PR",   "2016-06-25",
  ) %>% mutate(
    ADT = lubridate::as_date(ADT)
  )

  actual_output <- filter_pd(
    dataset = adrs,
    filter = PARAMCD == "OVR",
    source_pd = admiral::date_source(
      dataset_name = "adsl",
      date = PDDT
    ),
    source_datasets = list(adsl = adsl)
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID", "ADT")
  )
})

## filter_pd Test 3: first PD in input dataset ----
test_that("filter_pd Test 3: first PD in input dataset", {
  adrs <- tibble::tribble(
    ~STUDYID,       ~USUBJID,      ~PARAMCD, ~AVALC, ~ADT,
    "CDISCPILOT01", "01-701-1015", "OVR",    "CR",   "2016-01-25",
    "CDISCPILOT01", "01-701-1015", "OVR",    "PD",   "2016-02-22",
    "CDISCPILOT01", "01-701-1015", "BOR",    "CR",   "2016-01-25",
    "CDISCPILOT01", "01-701-1034", "OVR",    "SD",   "2015-12-07",
    "CDISCPILOT01", "01-701-1034", "OVR",    "PD",   "2016-04-25",
    "CDISCPILOT01", "01-701-1034", "OVR",    "PD",   "2016-06-25",
    "CDISCPILOT01", "01-701-1034", "BOR",    "SD",   "2015-12-07",
    "CDISCPILOT01", "01-701-1035", "OVR",    "SD",   "2016-04-25",
    "CDISCPILOT01", "01-701-1035", "OVR",    "PR",   "2016-06-25",
    "CDISCPILOT01", "01-701-1035", "BOR",    "PR",   "2016-06-25",
    "CDISCPILOT01", "01-701-1015", "PD",     "Y",    "2016-02-22",
    "CDISCPILOT01", "01-701-1034", "PD",     "Y",    "2016-04-25"
  ) %>% mutate(
    ADT = lubridate::as_date(ADT)
  )

  expected_output <- tibble::tribble(
    ~STUDYID,       ~USUBJID,      ~PARAMCD, ~AVALC, ~ADT,
    "CDISCPILOT01", "01-701-1015", "OVR",    "CR",   "2016-01-25",
    "CDISCPILOT01", "01-701-1015", "OVR",    "PD",   "2016-02-22",
    "CDISCPILOT01", "01-701-1034", "OVR",    "SD",   "2015-12-07",
    "CDISCPILOT01", "01-701-1034", "OVR",    "PD",   "2016-04-25",
    "CDISCPILOT01", "01-701-1035", "OVR",    "SD",   "2016-04-25",
    "CDISCPILOT01", "01-701-1035", "OVR",    "PR",   "2016-06-25",
  ) %>% mutate(
    ADT = lubridate::as_date(ADT)
  )

  actual_output <- filter_pd(
    dataset = adrs,
    filter = PARAMCD == "OVR",
    source_pd = admiral::date_source(
      dataset_name = "adrs",
      date = ADT,
      filter = PARAMCD == "PD",
    ),
    source_datasets = list(adrs = adrs)
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID", "ADT")
  )
})

## filter_pd Test 4: first PD derived from input dataset ----
test_that("filter_pd Test 4: first PD derived from input dataset", {
  adrs <- tibble::tribble(
    ~STUDYID,       ~USUBJID,      ~PARAMCD, ~AVALC, ~ADT,
    "CDISCPILOT01", "01-701-1015", "OVR",    "CR",   "2016-01-25",
    "CDISCPILOT01", "01-701-1015", "OVR",    "PD",   "2016-02-22",
    "CDISCPILOT01", "01-701-1015", "BOR",    "CR",   "2016-01-25",
    "CDISCPILOT01", "01-701-1034", "OVR",    "SD",   "2015-12-07",
    "CDISCPILOT01", "01-701-1034", "OVR",    "PD",   "2016-04-25",
    "CDISCPILOT01", "01-701-1034", "OVR",    "PD",   "2016-06-25",
    "CDISCPILOT01", "01-701-1034", "BOR",    "SD",   "2015-12-07",
    "CDISCPILOT01", "01-701-1035", "OVR",    "SD",   "2016-04-25",
    "CDISCPILOT01", "01-701-1035", "OVR",    "PR",   "2016-06-25",
    "CDISCPILOT01", "01-701-1035", "BOR",    "PR",   "2016-06-25"
  ) %>% mutate(
    ADT = lubridate::as_date(ADT)
  )

  expected_output <- tibble::tribble(
    ~STUDYID,       ~USUBJID,      ~PARAMCD, ~AVALC, ~ADT,
    "CDISCPILOT01", "01-701-1015", "OVR",    "CR",   "2016-01-25",
    "CDISCPILOT01", "01-701-1015", "OVR",    "PD",   "2016-02-22",
    "CDISCPILOT01", "01-701-1034", "OVR",    "SD",   "2015-12-07",
    "CDISCPILOT01", "01-701-1034", "OVR",    "PD",   "2016-04-25",
    "CDISCPILOT01", "01-701-1035", "OVR",    "SD",   "2016-04-25",
    "CDISCPILOT01", "01-701-1035", "OVR",    "PR",   "2016-06-25",
  ) %>% mutate(
    ADT = lubridate::as_date(ADT)
  )

  actual_output <- filter_pd(
    dataset = adrs,
    filter = PARAMCD == "OVR",
    source_pd = admiral::date_source(
      dataset_name = "adrs",
      date = ADT,
      filter = PARAMCD == "OVR" & AVALC == "PD",
    ),
    source_datasets = list(adrs = adrs)
  )

  expect_dfs_equal(
    expected_output,
    actual_output,
    keys = c("STUDYID", "USUBJID", "ADT")
  )
})

## filter_pd Test 5: error if invalid source_datasets ----
test_that("filter_pd Test 5: error if invalid source_datasets", {
  adrs <- tibble::tribble(
    ~STUDYID,       ~USUBJID,      ~PARAMCD, ~AVALC, ~ADT,
    "CDISCPILOT01", "01-701-1015", "OVR",    "CR",   "2016-01-25",
    "CDISCPILOT01", "01-701-1015", "OVR",    "PD",   "2016-02-22",
    "CDISCPILOT01", "01-701-1015", "BOR",    "CR",   "2016-01-25",
    "CDISCPILOT01", "01-701-1034", "OVR",    "SD",   "2015-12-07",
    "CDISCPILOT01", "01-701-1034", "OVR",    "PD",   "2016-04-25",
    "CDISCPILOT01", "01-701-1034", "OVR",    "PD",   "2016-06-25",
    "CDISCPILOT01", "01-701-1034", "BOR",    "SD",   "2015-12-07",
    "CDISCPILOT01", "01-701-1035", "OVR",    "SD",   "2016-04-25",
    "CDISCPILOT01", "01-701-1035", "OVR",    "PR",   "2016-06-25",
    "CDISCPILOT01", "01-701-1035", "BOR",    "PR",   "2016-06-25"
  ) %>% mutate(
    ADT = lubridate::as_date(ADT)
  )

  expect_error(
    filter_pd(
      dataset = adrs,
      filter = PARAMCD == "OVR",
      source_pd = admiral::date_source(
        dataset_name = "adrs",
        date = ADT,
        filter = PARAMCD == "OVR" & AVALC == "PD",
      ),
      source_datasets = list(ars = adrs)
    ),
    regexp = paste(
      "The dataset name specified for `source_pd` must be included in the list",
      " specified for the `source_datasets` parameter.",
      "Following names were provided by `source_datasets`:",
      "'ars'",
      sep = "\n"
    )
  )
})
