#' Simulate null distribution of representational similarity
#'
#' Randomly split and recombine the data from adult and child conditions to
#' create two mixed datasets. Then compute response profiles and their
#' pairwise correlations. Finally, assess representational similarity between
#' the two correlation matrices.
#'
#' @param wide_assoc A list of two data frames, where each row contains a cue
#'     and three responses from one participant.
#' @return a correlation coefficient.
#'
#' @details
#' This function should be run many times to simulate a distribution of
#' correlation values under the null hypothesis that adult and child responses
#' are sampled from the same population.
simulate_null_cor <- function(wide_assoc) {
    require(dplyr)
    rsm <- wide_assoc %>%
            dplyr::bind_rows() %>%
            dplyr::group_by(CUE, COND) %>% # nolint
            dplyr::mutate(
                gmixed = sample(gl(2, ceiling(dplyr::n() / 2)))[seq_len(n())]
            ) %>%
            dplyr::ungroup() %>%
            tidyr::pivot_responses_long() %>% # nolint
            dplyr::group_by(gmixed) %>% # nolint
            dplyr::group_split() %>%
            purrr::map(generate_response_profiles) %>% # nolint
            purrr::map(correlate_response_profiles) # nolint
    return(cor(rsm[[1]], rsm[[2]], method = "spearman"))
}
