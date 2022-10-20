correlate_response_profiles <- function(x) {
    as.dist(cor(x > 0))
}