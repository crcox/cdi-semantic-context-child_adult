generate_response_profiles <- function(.data) {
    xtabs(~ RESPONSE + CUE, data = .data)
}
