#' @importFrom dplyr mutate mutate_at select full_join left_join case_when
#' filter arrange desc slice everything group_by ungroup pull
#' @importFrom tibble tibble as_tibble
#' @importFrom purrr map map2 map_dbl map_dfr map2_dbl imap
#' discard keep compact some every
#' @importFrom tibble tibble as_tibble
#' @importFrom tidyr nest unnest pivot_longer pivot_wider
#' @importFrom rlang .data .env %||% set_names sym syms parse_expr :=
#' @keywords internal
"_PACKAGE"

# allowing for the use of the dot when piping
utils::globalVariables(".")

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
## usethis namespace: end
NULL
