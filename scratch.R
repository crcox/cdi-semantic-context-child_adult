library(dplyr)
library(tidyr)
library(purrr)
library(furrr)
library(future)
library(progressr)
library(ggplot2)

# Helper functions ----
load_to_list <- function(...) {
    files <- rlang::list2(...)
    map(files, function(filename) {
        e <- new.env()
        nm <- load(filename, envir = e)
        e[[nm]]
    })
}
selection_helper <- function(ix, x) {
    x[ix]
}


# Load data ----
load("data/cdi-metadata-preproc.Rdata")
associations <- load_to_list(
    adult_preproc = "./data/associations-adult-preproc.Rdata",
    child_preproc = "./data/associations-child-preproc.Rdata"
)
simulation_results <- load_to_list(
    all_cat = "simulation_results_30_all-categories.Rdata",
    high_ratio = "simulation_results_30_high_ratio.Rdata",
    high_ratio_lm = "simulation_results_30_high_ratio_lm.Rdata",
    high_glmprob = "simulation_results_30_high_glmprob.Rdata",
    high_glmprob_z = "simulation_results_30_high_glmprob.Rdata"
)
meta <- tibble(cue = levels(associations$adult_preproc$CUE)) %>%
    left_join(associations$adult_preproc %>%
        select(cue = CUE, num_item_id) %>%
        unique()
    ) %>%
    left_join(cdi_metadata_preproc %>% select(num_item_id, category))
categories <- unique(cdi_metadata_preproc$category)
category_subsets <- list(
    all_categories = categories,
    high_ratio = c("action_words", "body_parts", "descriptive_words", "games_routines", "household", "outside", "people", "places", "toys", "sounds"),
    high_ratio_lm = c("toys, games_routines, people, body_parts, furniture_rooms, outside, descriptive_words, helping_verbs, question_words, sounds"),
    high_glmprob = c("vehicles", "animals", "pronouns", "toys", "body_parts", "helping_verbs", "sounds", "descriptive_words", "furniture_rooms"),
    high_glmprob_z = c("toys", "games_routines", "people", "body_parts", "furniture_rooms", "outside", "descriptive_words", "helping_verbs", "question_words")
)


# Add category information to simulation results ----
# Also compute differences and z-values relating "real" to null correlations
simulation_results <- simulation_results %>%
    map_dfr(~{
        .x %>%
        mutate(
            num_item_id = map(ix, selection_helper, x = meta$num_item_id),
            cue = map(ix, selection_helper, x = meta$cue),
            category = map(ix, selection_helper, x = factor(meta$category, levels = categories)),
            diff = real - map_dbl(null, mean),
            zval = diff / map_dbl(null, sd),
            index = seq_len(n())
    )

    }, .id = "condition")

p_hist <- ggplot(simulation_results %>% filter(condition %in% c("all_cat", "high_ratio_lm")), aes(x = zval, fill = condition)) +
    geom_histogram(position = position_identity(), alpha = .4) +
    theme_bw(base_size = 18) +
    theme(legend.position = "none")
ggsave("histograms_lm.pdf", p_hist, width = 4, height = 6, units = "in", dpi = 300)

p_box <- ggplot(simulation_results %>% filter(condition %in% c("all_cat", "high_ratio_lm")), aes(x = condition, y = zval, color = condition)) +
    geom_boxplot() +
    theme_bw(base_size = 18) +
    theme(legend.position = "none")
ggsave("boxplots_lm.pdf", p_box, width = 4, height = 6, units = "in", dpi = 300)

# Simulations show that picking restricting to categories with high ratios
# leads to a larger shift between conditions, rather than modeling the
# likelihood of a category appearing in a set (as a function of the difference)
# of words with glm.


# Count categories ----
simulation_results <- simulation_results %>%
    mutate(
        catcount = map(category, ~{
            c(table(.x))
        })
    )


cdi_cat_counts <- meta %>%
    mutate(category = factor(category, levels = categories)) %>%
    group_by(category) %>%
    count() %>%
    ungroup() %>%
    mutate(p = n / sum(n))

catcount_w <- simulation_results$catcount %>%
    bind_rows() %>%
    mutate(
        real = simulation_results$real,
        diff = simulation_results$real - map_dbl(simulation_results$null, mean),
        zval = diff / map_dbl(simulation_results$null, sd),
        index = seq_len(n()),
        condition = simulation_results$condition
    )

catcount <- catcount_w %>%
    pivot_longer(
        cols = -c(real, diff, zval, index, condition),
        names_to = "category",
        values_to = "n"
    ) %>%
    left_join(cdi_cat_counts %>% rename(n_cdi = n, p_cdi = p)) %>%
    group_by(condition, index) %>%
    mutate(
        p = n / sum(n),
        ratio = p / p_cdi,
    ) %>%
    ungroup()


model_tbl <- catcount %>%
    filter(condition == "all_cat") %>%
    group_by(category) %>%
    group_split() %>%
    map_dfr(~{
        m <- lm(ratio ~ zval, data = .x)
        as_tibble(as.list(coef(m)))
    })
model_tbl <- catcount %>%
    group_by(category) %>%
    group_keys() %>%
    bind_cols(model_tbl) %>%
    arrange(zval)

ggplot(catcount %>% filter(condition == "all_cat") %>% mutate(category = factor(category, levels = model_tbl$category)), aes(x = zval, y = ratio, color = category)) +
    geom_smooth() +
    facet_wrap("category") +
    theme(legend.position = "none")

ggsave("./ratio_zval_category_all_cats.pdf", width = 6, height = 6, units = "in", dpi = 300)


# Model the likelihood of category inclusion ----
# in a set as a function of z-scored difference of real and null correlations.
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
