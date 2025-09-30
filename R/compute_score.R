compute_score <- function(Ug, Vu) {
  score <- as.numeric(Ug^2 / Vu)
  negative_log10_pval <- -pchisq(score, df = 1, lower.tail = FALSE, log.p = TRUE) / log(10)
  list(score = score, negative_log10_pval = negative_log10_pval)
}
