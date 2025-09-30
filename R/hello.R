#' Say hello
#'
#' @param name Character, someone to greet.
#' @return A greeting string.
#' @export
hello <- function(name = "world") {
  paste0("Hello, ", name, "!")
}
