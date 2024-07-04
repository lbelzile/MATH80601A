eightschools <- data.frame(
  id=factor(c("A","B","C","D","E","F","G","H")),
  effect = c(28.39,7.94,-2.75,6.82,-.64,.63,18.01,12.16),
  stderr = c(14.9, 10.2, 16.3, 11.0, 9.4, 11.4, 10.4, 17.6))
data <- list(y = eightschools$effect, 
             se_y = eightschools$stderr,
             N = 8L)
mod1 <- cmdstanr::cmdstan_model("models/eightschool_invgamma1.stan")
mod2 <- cmdstanr::cmdstan_model("models/eightschool_invgamma2.stan")
mod3 <- cmdstanr::cmdstan_model("models/eightschool_student.stan")
samp1 <- mod1$sample(data = data, seed = 202309,
            chains = 2,
            iter_sampling = 2000L)
samp2 <- mod2$sample(data = data, seed = 202309,
                     chains = 2,
                     iter_sampling = 2000L)
samp3 <- mod3$sample(data = data, seed = 202309,
                     chains = 2,
                     iter_sampling = 2000L)
p1 <- bayesplot::mcmc_dens(samp1$draws(variables = "tau"))
p2 <- bayesplot::mcmc_dens(samp2$draws(variables = "tau"))
p3 <- bayesplot::mcmc_dens(samp3$draws(variables = "tau"))

g1 <- ggplot(data = data.frame(tau = c(
  as.numeric(samp1$draws(variables = "tau")),
  as.numeric(samp2$draws(variables = "tau")),
  as.numeric(samp3$draws(variables = "tau"))
  ),
  prior = factor(rep(c("inverse gamma(1,1)",
            "inverse gamma(0.001,0.001)",
            "student-t(3)"), each = 4e3L))),
  mapping = aes(x = tau, group = prior, linetype = prior)) +
  geom_density() +
  scale_y_continuous(limits = c(0, 0.6), expand = c(0,0)) +
  scale_x_continuous(limits = c(0,32), expand = c(0,0)) + 
  labs(x = expression(tau), y = "", subtitle = "posterior density of random effect scale") +
  theme_classic() +
  theme(legend.position = "bottom")

png(filename = "eightschools-stddev-raneff.png", 
    width = 1200, height = 600,
    bg = "transparent",
    res = 200)
print(g1)
dev.off()
