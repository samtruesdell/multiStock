

# generate a random seed to use

seed <- runif(1, 0, 1e6)

cat(seed, file='ms_seed.txt')
