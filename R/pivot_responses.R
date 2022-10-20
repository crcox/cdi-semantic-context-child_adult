pivot_responses_wide <- function(d) {
    tidyr::pivot_wider(
        data = d,
        id = tidyselect::everything(),
        names_prefix = "R",
        names_from = "RESP_ID",
        values_from = "RESPONSE"
    )
}

pivot_responses_long <- function(d) {
    tidyr::pivot_longer(
        data = d,
        cols = c("R1", "R2", "R3"),
        names_prefix = "R",
        names_to = "RESP_ID",
        values_to = "RESPONSE"
    )
}
