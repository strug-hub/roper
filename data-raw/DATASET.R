## code to prepare `DATASET` dataset goes here

count <- matrix(rnbinom(n = 1e5, mu = 100, size = 1 / 0.5), ncol = 100)



usethis::use_data(DATASET, overwrite = TRUE)
