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
      model1 = c(intercept = 7.6655, msa = 1, age = 0, sex = 0),
      model2 = c(intercept = 11.6023, msa = 0.915845, age = 0, sex = 0),
      model3 = c(intercept = 10.3292, msa = 0.587056, age = 0.338519, sex = 1.29732)
    ),
    Hannum = list(
      model1 = c(intercept = 4.36099, msa = 1, age = 0, sex = 0),
      model2 = c(intercept = 5.32496, msa = 0.978194, age = 0, sex = 0),
      model3 = c(intercept = 4.89555, msa = 0.77497, age = 0.188571, sex = 0.968285)
    ),
    PhenoAge = list(
      model1 = c(intercept = 5.08917, msa = 1, age = 0, sex = 0),
      model2 = c(intercept = 7.03163, msa = 0.945517, age = 0, sex = 0),
      model3 = c(intercept = 0.712576, msa = 0.637665, age = 0.365991, sex = 0.120119)
    ),
    GrimAge = list(
      model1 = c(intercept = 1.30007, msa = 1, age = 0, sex = 0),
      model2 = c(intercept = 1.19803, msa = 1.00219, age = 0, sex = 0),
      model3 = c(intercept = 2.1728, msa = 0.861818, age = 0.121838, sex = -0.297574)
    ),
    `GrimAge v2` = list(
      model1 = c(intercept = 0.634309, msa = 1, age = 0, sex = 0),
      model2 = c(intercept = 1.40727, msa = 0.985322, age = 0, sex = 0),
      model3 = c(intercept = 2.83338, msa = 0.8632, age = 0.109576, sex = -0.275966)
    ),
    DunedinPACE = list(
      model1 = c(intercept = 0.117078, msa = 1, age = 0, sex = 0),
      model2 = c(intercept = 0.30409, msa = 0.792716, age = 0, sex = 0),
      model3 = c(intercept = 0.288731, msa = 0.744042, age = 0.0011102, sex = 0.0126529)
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
