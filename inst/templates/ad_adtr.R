# Name: ADTR
#
# Label: Tumor Response Analysis Dataset
#
# Input: adsl, rs, tr, tu
library(admiral)
library(admiralonco)
library(pharmaversesdtm) # Contains example datasets from the CDISC pilot project
# pharmaverseadam contains example datasets generated from the CDISC pilot
# project SDTM ran through admiral templates
library(pharmaverseadam)
library(dplyr)
library(lubridate)
library(stringr)

# Load source datasets ----

# Use e.g. haven::read_sas to read in .sas7bdat, or other suitable functions
# as needed and assign to the variables below.
# For illustration purposes read in pharmaverse test data
data("adsl")
data("rs_onco_recist")
data("tu_onco_recist")
data("tr_onco_recist")
tu <- tu_onco_recist
tr <- tr_onco_recist
rs <- rs_onco_recist

tu <- convert_blanks_to_na(tu) %>%
  filter(TUEVAL == "INVESTIGATOR")
tr <- convert_blanks_to_na(tr) %>%
  filter(
    TREVAL == "INVESTIGATOR" & TRGRPID == "TARGET" & TRTESTCD %in% c("LDIAM", "LPERP")
  )
rs <- convert_blanks_to_na(rs)

# Get list of ADSL vars required for derivations
adsl_vars <- exprs(RANDDT)

# Join ADSL vars to TR
tr <- derive_vars_merged(
  tr,
  dataset_add = adsl,
  new_vars = adsl_vars,
  by_vars = get_admiral_option("subject_keys")
)

# Add variables from TU (location of tumor) ----
tr <- derive_vars_merged(
  tr,
  dataset_add = tu,
  new_vars = exprs(TULOC),
  by_vars = c(get_admiral_option("subject_keys"), exprs(TRLNKID = TULNKID))
) %>% mutate(
  TULOCGR1 = if_else(
    TULOC == "LYMPH NODE",
    "NODAL",
    "NON-NODAL"
  )
)

tr <- mutate(
  tr,
  LSEXP = TRLNKID,
  LSASS = if_else(!is.na(TRSTRESN), TRLNKID, NA_character_)
)

# Derive timing variables (ADT, ADY, AVISIT, AVISITN) ----
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
  ) %>%
  mutate(
    AVISIT = if_else(
      VISIT == "SCREENING",
      "BASELINE",
      VISIT
    ),
    AVISITN = if_else(
      AVISIT == "BASELINE",
      0,
      VISITNUM
    )
  )

# Derive parameters for lesion diameters (LDIAMn & NLDIAMn) ----
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
    PARCAT3 = "RECIST 1.1",
    AVAL = TRSTRESN,
    ANL01FL = if_else(!is.na(AVAL), "Y", NA_character_)
  ) %>%
  select(-tmp_lesion_nr)

# Derive parameter for sum of diameter (SDIAM) ----
adtr_sum <- derive_summary_records(
  dataset_add = adtr,
  by_vars = exprs(!!!get_admiral_option("subject_keys"), !!!adsl_vars, AVISIT, AVISITN),
  filter_add = (str_starts(PARAMCD, "LDIAM") & TULOCGR1 == "NON-NODAL") |
    (str_starts(PARAMCD, "NLDIAM") & TULOCGR1 == "NODAL"),
  set_values_to = exprs(
    AVAL = sum(AVAL, na.rm = TRUE),
    ADY = min(ADY, na.rm = TRUE),
    ADT = min(ADT, na.rm = TRUE),
    PARAMCD = "SDIAM",
    PARAM = "Target Lesions Sum of Diameters by Investigator",
    PARCAT1 = "Target Lesion(s)",
    PARCAT2 = "Investigator",
    PARCAT3 = "RECIST 1.1"
  )
)

# Derive analysis flag (ANL01FL) ----
adtr_sum <- adtr_sum %>%
  derive_var_merged_summary(
    dataset_add = adtr,
    by_vars = get_admiral_option("subject_keys"),
    filter_add = AVISIT == "BASELINE" &
      ((str_starts(PARAMCD, "LDIAM") & TULOCGR1 == "NON-NODAL") |
        (str_starts(PARAMCD, "NLDIAM") & TULOCGR1 == "NODAL")),
    new_vars = exprs(LSEXP = paste(sort(TRLNKID), collapse = ", "))
  ) %>%
  derive_var_merged_summary(
    dataset_add = adtr,
    by_vars = c(get_admiral_option("subject_keys"), exprs(AVISIT)),
    filter_add = ((str_starts(PARAMCD, "LDIAM") & TULOCGR1 == "NON-NODAL") |
      (str_starts(PARAMCD, "NLDIAM") & TULOCGR1 == "NODAL")) & ANL01FL == "Y",
    new_vars = exprs(LSASS = paste(sort(TRLNKID), collapse = ", "))
  ) %>%
  mutate(
    ANL01FL = if_else(LSEXP == LSASS, "Y", NA_character_)
  )

# Derive baseline (ABLFL, BASE) ----
adtr_sum <- adtr_sum %>%
  restrict_derivation(
    derivation = derive_var_extreme_flag,
    args = params(
      by_vars = get_admiral_option("subject_keys"),
      order = exprs(ADY),
      new_var = ABLFL,
      mode = "last"
    ),
    filter = ADY <= 1
  ) %>%
  derive_var_base(
    by_vars = get_admiral_option("subject_keys")
  )

# Derive nadir (NADIR) ----
adtr_sum <- adtr_sum %>%
  derive_vars_joined(
    dataset_add = adtr_sum,
    by_vars = get_admiral_option("subject_keys"),
    order = exprs(AVAL),
    new_vars = exprs(NADIR = AVAL),
    join_vars = exprs(ADY),
    join_type = "all",
    filter_add = ANL01FL == "Y",
    filter_join = ADY.join < ADY,
    mode = "first",
    check_type = "none"
  )

# Derive change from baseline/nadir (CHG, PCHG, CHGNAD, PCHGNAD) ----
adtr_sum <- adtr_sum %>%
  derive_var_chg() %>%
  derive_var_pchg() %>%
  mutate(
    CHGNAD = AVAL - NADIR,
    PCHGNAD = if_else(NADIR == 0, NA_real_, 100 * CHGNAD / NADIR)
  )

# Derive additional flag variables (PDFL) ----
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
    by_vars = c(get_admiral_option("subject_keys"), exprs(ADT)),
    new_var = PDFL,
    condition = RSSTRESC == "PD"
  )

# Derive analysis flags (ANLzzFL) ----
adtr_sum <- adtr_sum %>%
  derive_var_relative_flag(
    by_vars = get_admiral_option("subject_keys"),
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
      by_vars = get_admiral_option("subject_keys"),
      new_var = ANL02FL,
      order = exprs(PCHG),
      mode = "first",
      check_type = "none"
    ),
    filter = ANL01FL == "Y" & POSTRNDFL == "Y"
  ) %>%
  select(-POSTRNDFL) %>%
  restrict_derivation(
    derivation = derive_var_relative_flag,
    args = params(
      by_vars = get_admiral_option("subject_keys"),
      new_var = ANL03FL,
      condition = PDFL == "Y",
      order = exprs(ADY),
      mode = "first",
      selection = "before",
      inclusive = FALSE
    ),
    filter = ANL01FL == "Y" | PDFL == "Y"
  ) %>%
  mutate(
    ANL04FL = if_else(ANL01FL == "Y" | PDFL == "Y", "Y", NA_character_)
  )

# Derive analysis sequence number (ASEQ) ----
adtr <- adtr %>%
  bind_rows(adtr_sum) %>%
  derive_var_obs_number(
    by_vars = get_admiral_option("subject_keys"),
    order = exprs(PARAMCD, AVISITN, TRSEQ),
    check_type = "error"
  )

# Add ADSL variables ----
adtr <- adtr %>%
  derive_vars_merged(
    dataset_add = select(adsl, !!!negate_vars(adsl_vars)),
    by_vars = get_admiral_option("subject_keys")
  )

# Change to whichever directory you want to save the dataset in
dir <- tools::R_user_dir("admiralonco_templates_data", which = "cache")
if (!file.exists(dir)) {
  # Create the folder
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
}
save(adtr, file = file.path(dir, "adtr.rda"), compress = "bzip2")
