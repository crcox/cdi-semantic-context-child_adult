load_to_list <- function(...) {
    files <- rlang::list2(...)
    purrr::map(files, function(filename) {
        e <- new.env()
        nm <- load(filename, envir = e)
        e[[nm]]
    })
}
