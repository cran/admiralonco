#' Adds a Parameter for Best Overall Response (without confirmation)
#'
#' @description
#' `r lifecycle::badge("superseded")` The `derive_param_bor()` function has been
#' superseded in favor of `derive_extreme_event()`.
#'
#' Adds a parameter for best overall response, without confirmation, optionally up to
#' first progressive disease
#'
#' @details
#'    Calculates the best overall response (BOR) parameter, as detailed below.
#'
#'    Records after PD can be removed using the source_pd and source_datasets
#'    arguments.
#'
#'   Note:
#'   1. All `CR`, `PR` and `PD` response records are considered for Best Overall Response.
#'
#'   2. All `SD` or `NON-CR/NON-PD` records where `ADT` >= `reference_date` +
#'      `ref_start_window` are also considered for Best Overall Response.
#'
#'   3. Subjects with **ONLY** an `SD` or `NON-CR/NON-PD` records where `ADT` <
#'      `reference_date` + `ref_start_window` are assigned a Best Overall Response of `NE`.
#'
#'   4. The Best Response, from the records in steps 1 to 3, is then selected in the following
#'      order of preference: CR, PR, SD, NON-CR/NON-PD, PD, NE, MISSING
#'
#'   5. The `AVAL` column is added and set using the `aval_fun(AVALC)` function
#'
#'   6. The columns specified by the `set_values_to` parameter and records
#'      are added to the dataframe passed into the `dataset` argument
#'
#'  Note: Any responses of SD or NON-CR/NON-PD that occur before `reference_date` +
#'  `ref_start_window` are ignored in the calculation of BOR. All other responses are included
#'  in the calculation of BOR, irrespective of the number of days from the reference date.
#'
#'  Also Note: All columns from the input dataset are kept. For subjects with no records in
#'  the input dataset (after the filter is applied) all columns are kept from ADSL which are
#'  also in the input dataset.  Columns which are not to be populated for the new parameter
#'  or populated differently (e.g. `RSSTRESC`, `VISIT`, `PARCATy`, `ANLzzFL`, ...) should be
#'  overwritten using the `set_values_to` parameter.
#'
#' @param dataset The input dataframe from which the Best Overall Response will
#'                be derived from and added to.
#'
#'   The columns `PARAMCD`, `ADT`, and `AVALC`and the columns specified in
#'   `subject_keys` and `reference_date` are expected.
#'
#'   After applying `filter_source` and/or `source_pd` the column `ADT` and the
#'   columns specified by `subject_keys` must be a unique key of the dataframe.
#'
#'    *Permitted Values:* a `data.frame()` object
#'
#' @param dataset_adsl ADSL input dataset.
#'
#'    The columns specified in the subject_keys argument are expected. For each subject in
#'    the passed `dataset` a new row is added to the input `dataset`. Columns
#'    in `dataset_adsl` that also appear in `dataset` will be populated with the
#'    appropriate subject-specific value for these new rows.
#'
#'    *Permitted Values:* a `data.frame()` object
#'
#' @param filter_source Filter to be applied to `dataset` to derive the
#'                       Best Overall Response
#'
#' @param source_pd Date of first progressive disease (PD)
#'
#'   If the parameter is specified, the observations of the input `dataset` for
#'   deriving the new parameter are restricted to observations up to the
#'   specified date. Observations at the specified date are included. For
#'   subjects without first PD date all observations are take into account.
#'
#'   *Permitted Values:* a `date_source` object (see `date_source()`
#'   for details)
#'
#' @param source_datasets Source dataframe to be used to calculate the
#'                        first PD date
#'
#'   A named list of dataframes is expected (although for BOR) only one dataframe is
#'   needed. It links the `dataset_name` from `source_pd` with an existing dataframe.
#'
#'   For example if `source_pd = pd_date` with
#'   ```{r, eval=FALSE}
#'   pd_date <- date_source(
#'     dataset_name = "adrs",
#'     date = ADT,
#'     filter = PARAMCD == PD
#'   )
#'   ```
#'   and the actual response dataframe in the script is `myadrs`, `source_datasets
#'   = list(adrs = myadrs)` should be specified.
#'
#' @param reference_date Reference date
#'
#'   The reference date is used along with `ref_start_window` to determine those
#'   records that occur before and after `ADT` (see Details section for further
#'   information). Usually it is treatment start date (`TRTSDT`) or
#'   randomization date (`RANDDT`).
#'
#'   *Permitted Values:* a numeric date column
#'
#' @param ref_start_window Stable disease time window
#'
#'   The ref_start_window is used along with `reference_date` to determine those
#'   records that occur before and after `ADT` (i.e. for a record determine
#'   whether `ADT` >= `reference_date` + `ref_start_window`),
#'   see Details section for further information.
#'
#'   *Permitted Values:* a non-negative numeric scalar
#'
#' @param missing_as_ne Consider no assessments as `"NE"`?
#'
#'   If the argument is set to `TRUE`, the response is set to `"NE"` for
#'   subjects in `dataset_adsl` without an assessment in the `dataset` after the
#'   filter has been applied. Otherwise, the response is set to `"MISSING"`
#'   for these subjects.
#'
#'   *Permitted Values:* a logical scalar
#'
#' @param aval_fun *Deprecated*, please use `set_values_to` instead.
#'
#'   Function to map character analysis value (`AVALC`) to numeric analysis
#'   value (`AVAL`)
#'
#'   The (first) argument of the function must expect a character vector and the
#'   function must return a numeric vector.
#'
#' @param set_values_to New columns to set
#'
#'   A named list returned by `exprs()` defining the columns to be set for the
#'   new parameter, e.g. `exprs(PARAMCD = "BOR", PARAM = "Best Overall
#'   Response")` is expected. The values must be symbols, character strings,
#'   numeric values, or `NA`.
#'
#' @param subject_keys Columns to uniquely identify a subject
#'
#'   *Permitted Values:* A list of symbols created using `exprs()`.
#'
#' @examples
#'
#' library(magrittr)
#' library(dplyr)
#' library(tibble)
#' library(lubridate)
#' library(admiral)
#'
#' # Create ADSL dataset
#' adsl <- tribble(
#'   ~USUBJID, ~TRTSDTC,
#'   "1",      "2020-01-01",
#'   "2",      "2019-12-12",
#'   "3",      "2019-11-11",
#'   "4",      "2019-12-30",
#'   "5",      "2020-01-01",
#'   "6",      "2020-02-02",
#'   "7",      "2020-02-02",
#'   "8",      "2020-04-01"
#' ) %>%
#'   mutate(
#'     TRTSDT = ymd(TRTSDTC),
#'     STUDYID = "XX1234"
#'   )
#'
#' # Create ADRS dataset
#' ovr_obs <- tribble(
#'   ~USUBJID, ~ADTC, ~AVALC, ~ANL01FL,
#'   "1", "2020-01-01", "PR", "Y",
#'   "1", "2020-02-01", "CR", "Y",
#'   "1", "2020-02-16", "NE", "Y",
#'   "1", "2020-03-01", "CR", "Y",
#'   "1", "2020-04-01", "SD", "Y",
#'   "2", "2020-01-01", "SD", "Y",
#'   "2", "2020-02-01", "PR", "Y",
#'   "2", "2020-03-01", "SD", "Y",
#'   "2", "2020-03-13", "CR", "Y",
#'   "3", "2019-11-12", "CR", "Y",
#'   "3", "2019-12-02", "CR", "Y",
#'   "3", "2020-01-01", "SD", "Y",
#'   "4", "2020-01-01", "PR", "Y",
#'   "4", "2020-03-01", "SD", "N",
#'   "4", "2020-04-01", "SD", "Y",
#'   "4", "2020-05-01", "PR", "Y",
#'   "4", "2020-05-15", "NON-CR/NON-PD", "Y",
#'   "5", "2020-01-01", "PR", "Y",
#'   "5", "2020-01-10", "SD", "Y",
#'   "5", "2020-01-20", "PR", "Y",
#'   "5", "2020-05-15", "NON-CR/NON-PD", "Y",
#'   "6", "2020-02-06", "PR", "Y",
#'   "6", "2020-02-16", "CR", "Y",
#'   "6", "2020-03-30", "PR", "Y",
#'   "6", "2020-04-12", "PD", "Y",
#'   "6", "2020-05-01", "CR", "Y",
#'   "6", "2020-06-01", "CR", "Y",
#'   "7", "2020-02-06", "PR", "Y",
#'   "7", "2020-02-16", "CR", "Y",
#'   "7", "2020-04-01", "NE", "N"
#' ) %>%
#'   mutate(PARAMCD = "OVR")
#'
#' pd_obs <-
#'   bind_rows(tribble(
#'     ~USUBJID, ~ADTC,        ~AVALC,
#'     "2",      "2020-03-01", "Y",
#'     "4",      "2020-02-01", "Y"
#'   ) %>%
#'     mutate(PARAMCD = "PD"))
#'
#' adrs <- bind_rows(ovr_obs, pd_obs) %>%
#'   mutate(
#'     ADT = ymd(ADTC),
#'     STUDYID = "XX1234"
#'   ) %>%
#'   select(-ADTC) %>%
#'   derive_vars_merged(
#'     dataset_add = adsl,
#'     by_vars     = exprs(STUDYID, USUBJID),
#'     new_vars    = exprs(TRTSDT)
#'   )
#'
#' pd_date <- date_source(
#'   dataset_name = "adrs",
#'   date         = ADT,
#'   filter       = PARAMCD == "PD"
#' )
#'
#' aval_fun_pass <- function(arg) {
#'   case_when(
#'     arg == "CR" ~ 11,
#'     arg == "PR" ~ 22,
#'     arg == "SD" ~ 33,
#'     arg == "NON-CR/NON-PD" ~ 44,
#'     arg == "PD" ~ 55,
#'     arg == "NE" ~ 66,
#'     arg == "MISSING" ~ 77,
#'     TRUE ~ NA_real_
#'   )
#' }
#'
#' # Derive best overall response parameter
#' derive_param_bor(
#'   adrs,
#'   dataset_adsl = adsl,
#'   filter_source = PARAMCD == "OVR" & ANL01FL == "Y",
#'   source_pd = pd_date,
#'   source_datasets = list(adrs = adrs),
#'   aval_fun = aval_fun_pass,
#'   reference_date = TRTSDT,
#'   ref_start_window = 28,
#'   set_values_to = exprs(
#'     PARAMCD = "BOR",
#'     PARAM = "Best Overall Response"
#'   )
#' ) %>%
#'   filter(PARAMCD == "BOR")
#' @export
#'
#' @author Stephen Gormley
#'
#' @family superseded
#' @keywords superseded
#'
#' @return The dataframe passed in the `dataset` argument with additional columns and/or
#'         rows as set in the `set_values_to` argument.

derive_param_bor <- function(dataset,
                             dataset_adsl,
                             filter_source,
                             source_pd = NULL,
                             source_datasets = NULL,
                             reference_date,
                             ref_start_window,
                             missing_as_ne = FALSE,
                             aval_fun,
                             set_values_to,
                             subject_keys = get_admiral_option("subject_keys")) {
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Assert statements (checked in order of signature) ----
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  reference_date <- assert_symbol(arg = enexpr(reference_date))

  assert_vars(arg = subject_keys)

  assert_data_frame(
    arg = dataset,
    required_vars = expr_c(
      subject_keys,
      reference_date,
      exprs(PARAMCD, ADT, AVALC)
    )
  )

  assert_data_frame(
    arg           = dataset_adsl,
    required_vars = expr_c(subject_keys)
  )

  filter_source <- assert_filter_cond(arg = enexpr(filter_source))

  assert_integer_scalar(
    arg      = ref_start_window,
    subset   = "non-negative"
  )

  assert_logical_scalar(arg = missing_as_ne)

  assert_varval_list(
    arg               = set_values_to,
    required_elements = c("PARAMCD")
  )

  assert_param_does_not_exist(
    dataset = dataset,
    param   = set_values_to$PARAMCD
  )

  if (!missing(aval_fun)) {
    deprecate_warn("0.4.0", "derive_param_bor(aval_fun = )", "derive_param_bor(set_values_to = )")
    set_values_to <- exprs(!!!set_values_to, AVAL = {{ aval_fun }}(AVALC))
  }

  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # filter_pd and filter_source: Filter source dataset using filter_source----
  # argument and also filter data after progressive disease with filter_pd
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  if (!is.null(source_pd)) {
    dataset_filter <- dataset %>%
      filter_pd(
        filter          = !!filter_source,
        source_pd       = source_pd,
        source_datasets = source_datasets,
        subject_keys    = subject_keys
      )
  } else {
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    # filter_source: Filter using filter_source argument ----
    # This would also be used to filter out records from dataset that are greater
    # than e.g. ADSL.TRTSDT
    # Not filtering data after progressive disease with filter_pd
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    dataset_filter <- dataset %>%
      filter(!!filter_source)
  }

  # Error if filter results in 0 records
  if (nrow(dataset_filter) == 0) {
    err_msg <- sprintf(
      "dataframe passed into %s argument with the filter %s has 0 records",
      "dataset",
      deparse(filter_source)
    )

    abort(err_msg)
  }

  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Create Sort Order for Selection of Minimum Later ----
  #
  #   Note:
  #   1. All `CR`, `PR` and `PD` response records are considered for Best Overall
  #      Response.
  #
  #   2. All `SD` or `NON-CR/NON-PD` records where `ADT` >= `reference_date` +
  #     `ref_start_window` are also considered for Best Overall Response.
  #
  #  3. Subjects with **ONLY** `SD` or `NON-CR/NON-PD` records where `ADT` <
  #     `reference_date` + `ref_start_window` would give Best Overall Response
  #     of `NE`.
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  dataset_ordered <- dataset_filter %>%
    mutate(
      AVALC = if_else(
        AVALC %in% c("SD", "NON-CR/NON-PD") & ADT < !!reference_date + days(ref_start_window),
        "NE",
        AVALC
      ),
      tmp_order = case_when(
        AVALC %in% c("CR") ~ 1,
        AVALC %in% c("PR") ~ 2,
        AVALC %in% c("SD") ~ 3,
        AVALC %in% c("NON-CR/NON-PD") ~ 4,
        AVALC %in% c("PD") ~ 5,
        AVALC %in% c("NE") ~ 6,
        is.null(AVALC) ~ 7
      )
    )

  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # adsl only subjects ----
  # Note Requirement: For subjects without observations in the input dataset
  # after the filter is applied, we keep all columns from ADSL which
  # are also in the input dataset.
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  adsl_data <- dataset_adsl %>%
    select(intersect(
      colnames(dataset_adsl),
      colnames(dataset)
    )) %>%
    mutate(
      AVALC = case_when(
        isTRUE(missing_as_ne) ~ "NE",
        TRUE ~ "MISSING"
      ),
      tmp_order = 999
    )

  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Bind two types of dataframes and select lowest value as BOR
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  param_bor <- bind_rows(
    dataset_ordered,
    adsl_data
  ) %>%
    filter_extreme(
      by_vars = subject_keys,
      order = exprs(tmp_order, ADT),
      mode = "first"
    ) %>%
    select(-tmp_order) %>%
    process_set_values_to(set_values_to)

  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # Bind back to passed dataframe in dataset argument and return ----
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  bind_rows(
    dataset,
    param_bor
  )
}
