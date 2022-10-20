library(dplyr)
library(tidyr)
library(purrr)
library(furrr)
library(future)
library(progressr)
library(ggplot2)

source("R/load_to_list.R")
source("R/pivot_responses.R")
source("R/generate_response_profiles.R")
source("R/correlate_response_profiles.R")
source("R/simulate_null_cor.R")

# Establish "future" plan for parallel computing ----
# If Windows --> multisession
# If Linux --> multicore
# If you want serial (non-parallel processing) --> sequential
plan("multisession")

# Load association data ----
load("data/cdi-metadata-preproc.Rdata")
associations <- load_to_list(
    adult_preproc = "./data/associations-adult-preproc.Rdata",
    child_preproc = "./data/associations-child-preproc.Rdata"
)


# Compute true representational similarity matrices ----
repsim <- purrr::map(
    associations,
    function(d) {
        correlate_response_profiles(
            generate_response_profiles(d)
        )
    }
)
repsim_cor <- cor(repsim$adult, repsim$child, method = "spearman")

# Begin simulation ----

# Run in parallel (this will take several minutes)
with_progress({
    p <- progressor(steps = 1000)
    null_repsim_cor <- furrr::future_map_dbl(
        seq_len(1000),
        function(i, data) {
            p()
            simulate_null_cor(data)
        },
        data = purrr::map(associations, pivot_responses_wide),
        .options = furrr_options(seed = TRUE)
    )
})

save(null_repsim_cor, file = "null-repsim-cor-adult-child.Rdata")
save(repsim_cor, file = "repsim-cor-adult-child.Rdata")

get_lowertri <- function(x) {
    x[lower.tri(x)]
}

h <- function(associations) {
    associations %>%
        bind_rows() %>%
        pivot_wider(
            id = tidyselect::everything(),
            names_prefix = "R",
            names_from = "RESP_ID",
            values_from = "RESPONSE"
        ) %>%
        group_by(CUE, COND) %>% # nolint
        mutate(gmixed = sample(gl(2, ceiling(n() / 2)))[seq_len(n())]) %>%
        ungroup() %>%
        pivot_responses_long() %>% # nolint
        group_by(gmixed) %>% # nolint
        group_split() %>%
        map(function(d) {
            correlate_response_profiles( # nolint
                generate_response_profiles(d) # nolint
            )
        })
}
g <- function(rsm_list, ix) {
    rsm_list <- rsm_list %>%
        map(function(x, ix) {
            get_lowertri(as.matrix(x)[ix, ix])
        }, ix = ix)
    return(cor(rsm_list[[1]], rsm_list[[2]], method = "spearman"))
}
f <- function(rsm_real, rsm_null, n, pop = NULL) {
    if (is.null(pop)) {
        pop <- attr(rsm_real[[1]], "Size")
        ix <- sample(pop, size = n)
    } else {
        ix <- sample(pop, size = n)
    }
    return(tibble(
        real = g(rsm_real, ix),
        null = list(map_dbl(rsm_null, g, ix = ix)),
        ix = list(ix)))
}

# Get a set of null RSMs to sample from. ----
#with_progress({
#    iter <- 100
#    p <- progressor(steps = iter)
#    rsm_null_lists <- future_map(1:iter, function(i, x) {
#        p()
#        h(x)
#    },
#    x = associations,
#    .options = furrr_options(seed = TRUE))
#})
#save(rsm_null_lists, file = "rsm-null-lists.Rdata")
load("rsm-null-lists.Rdata")

# Run simulations ----
meta <- tibble(cue = levels(associations$adult_preproc$CUE)) %>%
    left_join(associations$adult_preproc %>%
        select(cue = CUE, num_item_id) %>%
        unique()
    ) %>%
    left_join(cdi_metadata_preproc %>% select(num_item_id, category))

categories <- unique(cdi_metadata_preproc$category)
category_subsets <- list(
    all_categories = categories,
    high_ratio = c(
        "action_words",
        "body_parts",
        "descriptive_words",
        "games_routines",
        "household",
        "outside",
        "people",
        "places",
        "toys",
        "sounds"
    ),
    high_ratio_lm = c(
        "toys",
        "games_routines",
        "people",
        "body_parts",
        "furniture_rooms",
        "outside",
        "descriptive_words",
        "helping_verbs",
        "question_words",
        "sounds"
    ),
    high_glmprob = c(
        "vehicles",
        "animals",
        "pronouns",
        "toys",
        "body_parts",
        "helping_verbs",
        "sounds",
        "descriptive_words",
        "furniture_rooms"
    ),
    high_glmprob_z = c(
        "toys",
        "games_routines",
        "people",
        "body_parts",
        "furniture_rooms",
        "outside",
        "descriptive_words",
        "helping_verbs",
        "question_words"
    )
)

for (i in seq_along(category_subsets)) {
    z <- meta$category %in% category_subsets[[i]]
    w <- meta$cue[z]
    pop <- which(labels(repsim[[1]]) %in% w)
    with_progress({
        iter <- 1e4
        p <- progressor(steps = iter)
        simulation_results <- future_map_dfr(1:iter, function(i, rsm_real, rsm_null, n, pop) {
            p()
            f(rsm_real, rsm_null, n, pop)
        },
        rsm_real = repsim,
        rsm_null = rsm_null_lists,
        n = 30,
        pop = pop,
        .options = furrr_options(seed = TRUE))
    })
    save(
        simulation_results,
        file = sprintf(
            "simulation_results_30_%s.Rdata",
            names(category_subsets)[i]
        )
    )
}

## -------
y$num_item_id <- map(y$ix, function(ix, num_item_ids) {
    num_item_ids[ix]
}, num_item_ids = meta$num_item_id)

q <- function(ix, x) {
    x[ix]
}
y <- y %>%
    mutate(
        num_item_id = map(.$ix, q, x = meta$num_item_id),
        cue = map(.$ix, q, x = meta$cue),
        category = map(
            .$ix,
            q,
            x = factor(meta$category, levels = categories)
        )
    ) %>%
    mutate(
        real = y$real,
        diff = y$real - map_dbl(y$null, mean),
        zval = diff / map_dbl(y$null, sd),
        index = seq_len(n())
    )

y$catcount <- map(y$category, ~{
    c(table(.x))
})

y$catprop <- map(y$catcount, ~{.x / 30})

y %>% arrange(real) %>% slice_min(real, n = 10)

catcount_w <- y$catcount %>%
    bind_rows() %>%
    mutate(
        real = y$real,
        diff = y$real - map_dbl(y$null, mean),
        zval = diff / map_dbl(y$null, sd),
        index = seq_len(n())
    )

catcount <- catcount_w %>%
    pivot_longer(
        cols = -c(real, diff, zval, index),
        names_to = "category",
        values_to = "count"
    )

catcount$p <- catcount$count / 30
catcount$w <- 30
m <- map(categories, ~{
    glm(
        p ~ zval,
        data = catcount,
        family = "binomial",
        weights = w,
        subset = (category == .x)
    )
})
names(m) <- categories
df <- imap_dfr(m, function(model, label) {
    tibble(
        category = factor(label, levels = categories),
        p = fitted(model),
        zval = model.frame(model)$zval
    )
})

ggplot(df, aes(x = zval, y = p, color = category)) +
    geom_line()

aa <- map_dfr(m, coef, .id = "category") %>% arrange(zval)

ggplot(df %>% filter(category %in% aa$category[1:9]), aes(x = zval, y = p, color = category)) +
    geom_line()

ggplot(df %>% filter(!(category %in% aa$category[1:9])), aes(x = zval, y = p, color = category)) +
    geom_line()

