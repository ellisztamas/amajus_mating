library("tidyverse")

mcmc <- list(
  read_tsv("005.results/003_mcmc_restrict_kurtosis/output/chain1.out") %>% mutate(chain = "1"),
  read_tsv("005.results/003_mcmc_restrict_kurtosis/output/chain2.out") %>% mutate(chain = "2"),
  read_tsv("005.results/003_mcmc_restrict_kurtosis/output/chain3.out") %>% mutate(chain = "3"),
  read_tsv("005.results/003_mcmc_restrict_kurtosis/output/chain4.out") %>% mutate(chain = "4")
) %>%
  do.call(what="rbind")

# This only really converges after 1500 iterations, because starting values for shape and scale are wacky
mcmc %>% 
  ggplot(aes(x = iter, y = log_posterior, colour = chain)) +
  geom_line()
mcmc %>% 
  ggplot(aes(x = iter, y = scale, colour = chain)) +
  geom_line()
mcmc %>% 
  ggplot(aes(x = iter, y = missing, colour = chain)) +
  geom_line()
mcmc %>% 
  ggplot(aes(x = iter, y = shape, colour = chain)) +
  geom_line()
mcmc %>% 
  ggplot(aes(x = iter, y = mixture, colour = chain)) +
  geom_line()

mcmc %>% 
  filter(iter > 1500) %>% 
  ggplot(aes(x=missing)) + 
  geom_histogram(aes(y=..density..)) + 
  stat_function(fun=dbeta,
                color="red",
                args=list(shape1 = 3, 
                          shape2 = 15)
                ) +
  lims(
    x = c(0,0.6)
  )

mcmc %>% 
  filter(iter > 1500) %>% 
  ggplot(aes(x=mixture)) + 
  geom_histogram(aes(y=..density..)) + 
  stat_function(fun=dbeta,
                color="red",
                args=list(shape1 = 1.1, 
                          shape2 = 1.1)
  ) +
  lims(
    x = c(0,1)
  )

mcmc %>% 
  filter(iter > 1500) %>% 
  ggplot(aes(x=shape)) + 
  geom_histogram(aes(y=..density.., colour = chain, fill=chain)) + 
  stat_function(fun=dgamma,
                color="red",
                args=list(shape = 10, 
                          scale = 1/5)
  )
mcmc %>% 
  filter(iter > 1500) %>% 
  ggplot(aes(x=scale)) + 
  geom_histogram(aes(y=..density..)) + 
  stat_function(fun=dgamma,
                color="red",
                args=list(shape = 6, 
                          scale = 50)
  )

# No correlation between missing and mixture.
mcmc %>% 
  filter(iter > 1500) %>% 
  ggplot(aes(x = missing, y = mixture)) + 
  geom_point()
# Shape and scale are strongly correlated. r = 0.92
mcmc %>% 
  filter(iter > 1500) %>% 
  ggplot(aes(x = shape, y = scale, colour = chain)) + 
  geom_point()
