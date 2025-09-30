merge_results <- function(res_list) {
  nm <- names(res_list[[1]])
  setNames(
    lapply(nm, function(x) do.call(c, lapply(res_list, `[[`, x))),
    nm
  )
}
