# Calibrate MSA-based epigenetic-clock values to the
# EPICv2-referenced scale.
#
# Usage:
# calibrate_msa_to_epicv2(
#   msa_value,
#   clock,
#   model = 1,
#   age = NULL,
#   sex = NULL
# )
#
# Arguments:
# msa_value: Numeric vector of MSA-based clock values.
# clock: One of "Horvath", "Hannum", "PhenoAge", "GrimAge",
#        "GrimAge v2", or "DunedinPACE".
# model: Calibration model, specified as 1, 2, or 3.
#        "model1", "model2", and "model3" are also accepted.
# age: Chronological age. Required only for Model 3.
# sex: Sex coded as 0 for female and 1 for male.
#      Required only for Model 3.
#
# Value:
# A numeric vector of calibrated MSA values on the
# EPICv2-referenced scale.

calibrate_msa_to_epicv2 <- function(
    msa_value,
    clock,
    model = "model1",
    age = NULL,
    sex = NULL
) {
  clocks <- c(
    "Horvath",
    "Hannum",
    "PhenoAge",
    "GrimAge",
    "GrimAge v2",
    "DunedinPACE"
  )

  clock <- match.arg(clock, clocks)

  model <- tolower(gsub("\\s+", "", as.character(model)[1]))
  if (model %in% c("1", "2", "3")) {
    model <- paste0("model", model)
  }
  model <- match.arg(model, c("model1", "model2", "model3"))

  coefs <- list(
    Horvath = list(
      model1 = c(intercept = 7.67, msa = 1, age = 0, sex = 0),
      model2 = c(intercept = 11.6, msa = 0.916, age = 0, sex = 0),
      model3 = c(intercept = 10.3, msa = 0.587, age = 0.339, sex = 1.30)
    ),
    Hannum = list(
      model1 = c(intercept = 4.36, msa = 1, age = 0, sex = 0),
      model2 = c(intercept = 5.33, msa = 0.978, age = 0, sex = 0),
      model3 = c(intercept = 4.90, msa = 0.775, age = 0.189, sex = 0.968)
    ),
    PhenoAge = list(
      model1 = c(intercept = 5.09, msa = 1, age = 0, sex = 0),
      model2 = c(intercept = 7.03, msa = 0.946, age = 0, sex = 0),
      model3 = c(intercept = 0.713, msa = 0.638, age = 0.366, sex = 0.120)
    ),
    GrimAge = list(
      model1 = c(intercept = 1.30, msa = 1, age = 0, sex = 0),
      model2 = c(intercept = 1.20, msa = 1.00, age = 0, sex = 0),
      model3 = c(intercept = 2.17, msa = 0.862, age = 0.122, sex = -0.298)
    ),
    `GrimAge v2` = list(
      model1 = c(intercept = 0.634, msa = 1, age = 0, sex = 0),
      model2 = c(intercept = 1.41, msa = 0.985, age = 0, sex = 0),
      model3 = c(intercept = 2.83, msa = 0.863, age = 0.110, sex = -0.276)
    ),
    DunedinPACE = list(
      model1 = c(intercept = 0.117, msa = 1, age = 0, sex = 0),
      model2 = c(intercept = 0.304, msa = 0.793, age = 0, sex = 0),
      model3 = c(intercept = 0.289, msa = 0.744, age = 0.0011, sex = 0.0127)
    )
  )

  if (model == "model3") {
    if (is.null(age) || is.null(sex)) {
      stop("Model 3 requires chronological age and sex.")
    }

    n <- length(msa_value)

    if (!(length(age) %in% c(1L, n)) ||
        !(length(sex) %in% c(1L, n))) {
      stop("age and sex must have length 1 or the same length as msa_value.")
    }

    if (any(!is.na(sex) & !(sex %in% c(0, 1)))) {
      stop("sex must be coded as 0 for female and 1 for male.")
    }
  } else {
    age <- 0
    sex <- 0
  }

  b <- coefs[[clock]][[model]]

  unname(
    b["intercept"] +
      b["msa"] * msa_value +
      b["age"] * age +
      b["sex"] * sex
  )
}
