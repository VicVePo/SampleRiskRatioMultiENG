#' Sample Size and Power Calculation for Cohort Studies
#' 
#' @param alpha Significance level
#' @param power Desired power (optional if n_exposed provided)
#' @param n_exposed Number of exposed subjects (optional if power provided)
#' @param p0 Incidence in unexposed group (required if method = "rho")
#' @param RR Relative Risk to detect (optional if p1 provided)
#' @param p1 Incidence in exposed group (optional if RR provided)
#' @param CI_upper Upper confidence interval limit of RR
#' @param CI_lower Lower confidence interval limit of RR (optional if SE provided)
#' @param SE Standard error of coefficient (optional if CI provided)
#' @param n_previous Sample size of previous study (required for method "se")
#' @param r Ratio of unexposed to exposed (n_unexposed / n_exposed)
#' @param method Calculation method: "rho" or "se" (standard error)
#' @param rho_values Vector of rho values for sample size calculation (if method = "rho")
#' @return A data frame with results according to chosen method
#' @export
SampleRiskRatioMultiENG <- function(alpha, power = NULL, n_exposed = NULL, 
                                    p0 = NULL, RR = NULL, p1 = NULL,
                                    CI_upper = NULL, CI_lower = NULL,
                                    SE = NULL, n_previous = NULL, r = 1,
                                    method = "rho",
                                    rho_values = seq(0, 0.9, by = 0.1)) {
  
  # Verifications according to method
  if (method == "rho") {
    if (is.null(p0)) {
      stop("For 'rho' method, p0 is required")
    }
    if (is.null(RR) && is.null(p1)) {
      stop("For 'rho' method, either RR or p1 must be provided")
    }
    if (!is.null(RR) && !is.null(p1)) {
      stop("Provide either RR or p1, not both")
    }
  } else if (method == "se") {
    if (is.null(RR)) {
      stop("For 'se' method, RR is required")
    }
    if (is.null(SE) && (is.null(CI_upper) || is.null(CI_lower))) {
      stop("For 'se' method, SE or both CI limits are required")
    }
    if (is.null(n_previous)) {
      stop("For 'se' method, n_previous is required")
    }
  } else {
    stop("Method must be 'rho' or 'se'")
  }
  
  # If CI provided, calculate SE
  if (!is.null(CI_upper) && !is.null(CI_lower)) {
    SE <- (log(CI_upper) - log(RR)) / qnorm(0.975)
    cat("Calculated Standard Error:", SE, "
")
  }
  
  if (method == "rho") {
    # Calculate p1 if RR provided
    if (!is.null(RR)) {
      p1 <- p0 * RR
    } else {
      RR <- p1/p0
    }
    
    z_alpha <- qnorm(1 - alpha/2)
    
    if (!is.null(power)) {
      # Calculate sample size
      z_beta <- qnorm(power)
      
      calculate_n <- function(rho) {
        numerator <- (z_alpha * sqrt((1 + 1/r) * p0 * (1 - p0)) + 
                       z_beta * sqrt(p1 * (1 - p1) + p0 * (1 - p0) / r))^2
        denominator <- (p1 - p0)^2 * (1 - rho^2)
        
        n_exposed <- ceiling(numerator / denominator)
        n_unexposed <- ceiling(n_exposed * r)
        
        return(c(n_exposed, n_unexposed))
      }
      
      results <- data.frame(
        rho = rho_values,
        n_exposed = sapply(rho_values, function(rho) calculate_n(rho)[1]),
        n_unexposed = sapply(rho_values, function(rho) calculate_n(rho)[2])
      )
      
      results$total_sample_size <- results$n_exposed + results$n_unexposed
    } else {
      # Calculate power
      n_unexposed <- ceiling(n_exposed * r)
      
      calculate_power <- function(rho) {
        numerator <- sqrt(n_exposed * (p1 - p0)^2 * (1 - rho^2)) - 
          z_alpha * sqrt((1 + 1/r) * p0 * (1 - p0))
        denominator <- sqrt(p1 * (1 - p1) + p0 * (1 - p0) / r)
        
        power <- pnorm(numerator / denominator)
        return(power)
      }
      
      results <- data.frame(
        rho = rho_values,
        power = sapply(rho_values, calculate_power)
      )
      
      results$n_exposed <- n_exposed
      results$n_unexposed <- n_unexposed
      results$total_sample_size <- n_exposed + n_unexposed
    }
  } else {
    # Standard error based method
    z_alpha <- qnorm(1 - alpha/2)
    if (!is.null(power)) {
      z_gamma <- qnorm(power)
      n <- ceiling((z_alpha + z_gamma)^2 * n_previous * SE^2 / (log(RR))^2)
      results <- data.frame(
        total_sample_size = n
      )
      results$n_exposed <- ceiling(n / (1 + r))
      results$n_unexposed <- ceiling(n * r / (1 + r))
    } else {
      z_gamma <- sqrt(n_exposed * (log(RR))^2 / (n_previous * SE^2)) - z_alpha
      power <- pnorm(z_gamma)
      results <- data.frame(
        power = power
      )
      results$n_exposed <- n_exposed
      results$n_unexposed <- ceiling(n_exposed * r)
      results$total_sample_size <- n_exposed + results$n_unexposed
    }
  }
  return(results)
}

#' Logistics Calculation for Cohort Studies
#'
#' @param n_exposed Number of exposed subjects required
#' @param n_unexposed Number of unexposed subjects required
#' @param recruitment_rate Recruitment rate per day
#' @param rejection_rate Expected rejection rate
#' @param loss_to_followup_rate Expected loss to follow-up rate
#' @param eligibility_rate Expected eligibility rate
#' @param working_days_month Number of working days per month
#' @return A list with logistics calculations for the study
#' @export
cohort_study_logistics <- function(n_exposed, 
                                 n_unexposed,
                                 recruitment_rate,
                                 rejection_rate,
                                 loss_to_followup_rate,
                                 eligibility_rate,
                                 working_days_month) {
  
  total_sample_size <- n_exposed + n_unexposed
  
  # Adjustment for losses to follow-up
  sample_with_losses <- total_sample_size / (1 - loss_to_followup_rate)
  
  # Adjustment for rejections
  sample_with_rejections <- sample_with_losses / (1 - rejection_rate)
  
  # Adjustment for eligibility
  sample_to_evaluate <- sample_with_rejections / eligibility_rate
  
  # Time calculations
  recruitment_days <- ceiling(sample_to_evaluate / recruitment_rate)
  recruitment_months <- recruitment_days / working_days_month
  
  # Results
  results <- list(
    final_sample = total_sample_size,
    sample_with_losses = ceiling(sample_with_losses),
    sample_to_enroll = ceiling(sample_with_rejections),
    sample_to_evaluate = ceiling(sample_to_evaluate),
    recruitment_days = recruitment_days,
    recruitment_months = round(recruitment_months, 2)
  )
  
  # Print summary
  cat("
Study logistics summary:
")
  cat("----------------------------------------
")
  cat("Required final sample:", total_sample_size, "
")
  cat("Sample considering losses:", ceiling(sample_with_losses),
      "(", loss_to_followup_rate*100, "% losses)
")
  cat("Sample to enroll:", ceiling(sample_with_rejections),
      "(", rejection_rate*100, "% rejection)
")
  cat("Sample to evaluate:", ceiling(sample_to_evaluate),
      "(", eligibility_rate*100, "% eligibility)
")
  cat("Days needed:", recruitment_days,
      "(", recruitment_rate, "persons per day)
")
  cat("Months needed:", round(recruitment_months, 2),
      "(", working_days_month, "working days per month)
")
  cat("----------------------------------------
")
  
  return(invisible(results))
}
