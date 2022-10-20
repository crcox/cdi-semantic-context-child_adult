RSM <- repsim %>% map(as.matrix)
str(RSM)

z <- lower.tri(RSM[[1]])

cor_profiles <- RSM %>%
    map(~{
        n <- nrow(.x)
        z <- !diag(n)
        print(str(z))
        matrix(.x[z], nrow = n - 1, ncol = n)
    })

rr <- numeric(598)
for (i in 1:598) {
    rr[i] <- cor(cor_profiles$adult_preproc[, i], cor_profiles$child_preproc[, i])
}

x_df <- tibble(word = names(rr), rr = rr)
x_df %>% slice_max(rr, n = 10)
x_df %>% slice_min(rr, n = 10)


mover_words <- x_df %>% slice_min(rr, n = 30) %>% pull(word)

load("rsm-null-lists.Rdata")
str(rsm_null_lists[1])
r_movers_null <- rsm_null_lists %>%
    map_dbl(function(.x, w) {
        v <- list(
            get_lowertri(as.matrix(.x[[1]])[w, w]),
            get_lowertri(as.matrix(.x[[2]])[w, w])
        )
        return(cor(v[[1]], v[[2]], method = "spearman"))
    }, w = mover_words)
