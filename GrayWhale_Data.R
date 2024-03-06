# Gray whale data
# Reads gray whale photogrammetry data from Josh Stewart (2024-03-05). Data were created by 
# Morgan Lynn.
# 

library(tidyverse)
library(readr)
library(ggplot2)
library(loo)
library(bayesplot)


gray.whale.data <- read_csv(file = "data/All GW Photogrammetry TE.csv", 
                            col_types = cols(DATE = col_date(format = "%m/%d/%y"),
                                             ID = col_character(),
                                             Length = col_double(),
                                             Max_Width = col_double(),
                                             L50W = col_double(),
                                             L60W = col_double(),
                                             MaxWL = col_double(),
                                             L50WL = col_double(),
                                             L60WL = col_double(),
                                             Migration = col_character(),
                                             Stage = col_character(),
                                             Notes = col_character())) %>%
  mutate(Year = year(DATE))

# Do some stats on these numbers
South.bound <- gray.whale.data %>% 
  group_by(Year) %>%
  select(Year, Length, Max_Width, Migration) %>%
  filter(Length > 7) %>%
  na.omit() %>%
  #mutate(f.Year = as.factor(Year)) %>%
  filter(Migration == "Southbound")
  
ggplot(South.bound %>% mutate(f.Year = as.factor(Year)), 
       aes(x = Length, y = Max_Width, group = Year)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~f.Year)

# https://agabrioblog.onrender.com/tutorial/ancova-jags/ancova-jags/
X <- model.matrix(~ f.Year + Length , 
                  data = South.bound %>% mutate(f.Year = as.factor(Year)))

# Run ANCOVA in JAGS
model.string <- "model{
for (i in 1:n) { y[i] ~ dnorm(mean[i], tau)
mean[i] <- inprod(beta[], X[i,]) 
log.lkhd[i] <- logdensity.norm(y[i], mean[i], tau)}
for (j in 1:n.groups) { beta[j] ~ dnorm(0, 0.00001)}
sigma ~ dgamma(0.1, 0.01)
tau <- 1/(sigma * sigma)
}"

write_lines(model.string, file = "models/ANCOVA.txt")

jags.data <- list(y = South.bound$Max_Width,
                  X = X,
                  n = nrow(X),
                  n.groups = ncol(X))

jags.params <- c("beta", "sigma", "log.lkhd")

MCMC.params <- list(n.samples = 55000,
                    n.chains = 5,
                    n.burnin = 5000,
                    n.thin = 5)

jm <- jagsUI::jags(jags.data,
                   inits = NULL,
                   parameters.to.save = jags.params,
                   model.file = "models/ANCOVA.txt",
                   n.chains = MCMC.params$n.chains,
                   n.burnin = MCMC.params$n.burnin,
                   n.thin = MCMC.params$n.thin,
                   n.iter = MCMC.params$n.samples,
                   DIC = T,
                   parallel = T)

# This is a handy tool to see posteriors
mcmc_areas(jm$samples, regex_pars = "beta|sigma")

# summary stats for beta parameters
jm$summary %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Parameter") -> summary.df

beta.summary <- summary.df[grep("beta", summary.df$Parameter), ]

# The intercept of 1988 is beta[1]. The means of other years are different from 1988 by beta[2] to beta[9].
# One foot increase in length in 1988 resulted in beta[10] change in maximum width.  

ggplot(South.bound %>% mutate(f.Year = as.factor(Year)), 
       aes(x = Length, 
           y = Max_Width, 
           group = f.Year, 
           color = f.Year)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Length (m)") +
  ylab("Maximum width (m)") +
  labs(color = "Year") +
  theme(legend.position = c(0.2, 0.8),
        legend.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(ncol = 3))
