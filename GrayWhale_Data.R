# Gray whale data
# Reads gray whale photogrammetry data from Josh Stewart (2024-03-05). Data were created by 
# Morgan Lynn.
# 

rm(list = ls())

library(tidyverse)
library(readr)
library(ggplot2)
library(loo)
library(bayesplot)

save.fig <- FALSE

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

South.bound %>% mutate(f.Year = as.factor(Year)) -> South.bound.1 

ggplot(South.bound.1, 
       aes(x = Length, y = Max_Width, group = Year)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~f.Year)

# Create a design matrix
X <- model.matrix(~ f.Year * Length , 
                  data = South.bound.1)

# https://agabrioblog.onrender.com/tutorial/ancova-jags/ancova-jags/
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
jags.params <- c("beta", "sigma", "log.lkhd")

MCMC.params <- list(n.samples = 55000,
                    n.chains = 5,
                    n.burnin = 5000,
                    n.thin = 5)

# South bound analysis
jags.data.SB <- list(y = South.bound.1$Max_Width,
                  X = X,
                  n = nrow(X),
                  n.groups = ncol(X))

SB.out.filename <- "Rdata/jm_out_SB.rds"
if (file.exists(SB.out.filename)){
  jm.south.bound <- read_rds(SB.out.filename)
} else {
  jm.south.bound <- jagsUI::jags(jags.data.SB,
                                 inits = NULL,
                                 parameters.to.save = jags.params,
                                 model.file = "models/ANCOVA.txt",
                                 n.chains = MCMC.params$n.chains,
                                 n.burnin = MCMC.params$n.burnin,
                                 n.thin = MCMC.params$n.thin,
                                 n.iter = MCMC.params$n.samples,
                                 DIC = T,
                                 parallel = T)
  
  write_rds(jm.south.bound,
            file = SB.out.filename)
  
}

# This is a handy tool to see posteriors
mcmc_areas(jm.south.bound$samples, regex_pars = "beta|sigma")

# summary stats for beta parameters
jm.south.bound$summary %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Parameter") -> SB.summary.df

SB.beta.summary <- SB.summary.df[grep("beta", SB.summary.df$Parameter), ]

SB.beta.intercept <- SB.beta.summary[1:9,]
SB.beta.slope <- SB.beta.summary[10:18,]

# I will just use the slopes and forget about the variance for now. It needs 
# either Hessian or create a var-cov matrix from posterior samples to compute
# variance of slopes. 
SB.slopes <- data.frame(Year = unique(South.bound.1$Year),
                        Slope = c(SB.beta.slope$mean[1],
                                  SB.beta.slope$mean[2:9] + SB.beta.slope$mean[1]),
                        Migration = "South")

p.Length_Width_south.bound <- ggplot(South.bound.1, 
                                     aes(x = Length, 
                                         y = Max_Width, 
                                         group = f.Year, 
                                         color = f.Year)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Length (m)") +
  ylab("Maximum width (m)") +
  labs(color = "Year") +
  ylim(c(1, 3)) +
  theme(legend.position.inside = c(0.2, 0.8),
        legend.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(ncol = 3)) +
  ggtitle("Southbound Migration")

if (save.fig)
  ggsave(p.Length_Width_south.bound,
         filename = "figures/length_width_SB.png",
         dpi = 600)

# Do the same with north bound. 
# # Do some stats on these numbers
North.bound <- gray.whale.data %>% 
  group_by(Year) %>%
  select(Year, Length, Max_Width, Migration) %>%
  filter(Length > 7) %>%
  na.omit() %>%
  #mutate(f.Year = as.factor(Year)) %>%
  filter(Migration == "Northbound")

North.bound %>% mutate(f.Year = as.factor(Year)) -> North.bound.1 

ggplot(North.bound.1, 
       aes(x = Length, y = Max_Width, group = Year)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~f.Year)

# Create a design matrix
# https://agabrioblog.onrender.com/tutorial/ancova-jags/ancova-jags/
X <- model.matrix(~ f.Year * Length , 
                  data = North.bound.1)

jags.data.NB <- list(y = North.bound.1$Max_Width,
                     X = X,
                     n = nrow(X),
                     n.groups = ncol(X))

NB.out.filename <- "Rdata/jm_out_NB.rds"
if (file.exists(NB.out.filename)){
  jm.north.bound <- read_rds(NB.out.filename)
} else {
  
  jm.north.bound <- jagsUI::jags(jags.data.NB,
                                 inits = NULL,
                                 parameters.to.save = jags.params,
                                 model.file = "models/ANCOVA.txt",
                                 n.chains = MCMC.params$n.chains,
                                 n.burnin = MCMC.params$n.burnin,
                                 n.thin = MCMC.params$n.thin,
                                 n.iter = MCMC.params$n.samples,
                                 DIC = T,
                                 parallel = T)
  write_rds(jm.north.bound,
            file = NB.out.filename)
  
}

# This is a handy tool to see posteriors
mcmc_areas(jm.north.bound$samples, regex_pars = "beta|sigma")

# summary stats for beta parameters
jm.north.bound$summary %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "Parameter") -> NB.summary.df

NB.beta.summary <- NB.summary.df[grep("beta", NB.summary.df$Parameter), ]

NB.beta.intercept <- NB.beta.summary[1:15,]
NB.beta.slope <- NB.beta.summary[16:30,]

# I will just use the slopes and forget about the variance for now. It needs 
# either Hessian or create a var-cov matrix from posterior samples to compute
# variance of slopes. 
NB.slopes <- data.frame(Year = unique(North.bound.1$Year),
                        Slope = c(NB.beta.slope$mean[1],
                                  NB.beta.slope$mean[2:15] + NB.beta.slope$mean[1]),
                        Mirgraion = "NOrth")

p.Length_Width_north.bound <- ggplot(North.bound.1, 
                                     aes(x = Length, 
                                         y = Max_Width, 
                                         group = f.Year, 
                                         color = f.Year)) +
  geom_point() +
  geom_smooth(method = "lm") +
  xlab("Length (m)") +
  ylab("Maximum width (m)") +
  labs(color = "Year") +
  ylim(c(1, 3)) +
  theme(legend.position.inside = c(0.2, 0.8),
        legend.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(ncol = 3)) +
  ggtitle("Northbound Migration")

if (save.fig)
  ggsave(p.Length_Width_north.bound,
         filename = "figures/length_width_NB.png",
         dpi = 600)

# Bring in the abundance and calf production data
# C:\Users\tomo.eguchi\Documents\R\Granite_Canyon_Counts\Data
Nhats.2024 <- read_csv(file = "~/R/Granite_Canyon_Counts/Data/all_estimates_2024.csv",
                       col_types = cols(ID = col_integer(),
                                        Year = col_integer(),
                                        Season = col_character(),
                                        Nhat = col_integer(),
                                        LCL = col_double(),
                                        UCL = col_double(), 
                                        Method = col_character(),
                                        Nhat.2023 = col_integer())) %>%
  transmute(Year = Year,
            Nhat = Nhat,
            centered.Nhat = (Nhat - mean(Nhat))/sqrt(var(Nhat)))

calf.production.2023 <- read_csv(file = "~/R/Piedras_Blancas_Calf/data/Calf_Estimates_v3_Mv1_2023-06-22.csv",
                                 col_types = cols(Mean = col_double(),
                                                  Median = col_double(),
                                                  LCL = col_double(),
                                                  UCL = col_double(),
                                                  Year = col_integer(),
                                                  Method = col_character())) %>%
  transmute(Year = Year,
            Calf = Median,
            centered.calf = (Calf - mean(Calf))/sqrt(var(Calf)))

calf.production.2023 %>% 
  left_join(Nhats.2024, by = "Year")   -> obsd.N

# ggplot() +
#   geom_path(data = obsd.N, aes(x = Year, y = centered.calf),
#              color = "green", size = 2) +
#   geom_path(data = obsd.N, aes(x = Year, y = centered.Nhat),
#              color = "red", size = 2) +
#   geom_path(data = slopes, aes(x = Year, y = (Slope- mean(Slope))/sqrt(var(Slope))),
#              color = "gold", size = 2) +
#   geom_path(data = north.bound.slopes, aes(x = Year, y = (Slope - mean(Slope))/sqrt(var(Slope))),
#              color = "purple", size = 2) 
# 
# slopes %>%
#   mutate(std.SB.slope = (Slope - mean(Slope))/sqrt(var(Slope))) -> south.bound.slopes
# 
# north.bound.slopes %>%
#   mutate(std.NB.slope = (Slope - mean(Slope))/sqrt(var(Slope))) -> north.bound.slopes
# 
# obsd.N %>% left_join(south.bound.slopes, by = "Year") %>%
#   left_join(north.bound.slopes, by = "Year") -> all.stats
# 
# ggplot(all.stats) +
#   geom_point(aes(x = centered.calf, y = centered.Nhat),
#              size = 2) +
#   geom_point(aes(x = centered.calf, y = std.SB.slope),
#              color = "red", size = 2) +
#   geom_point(aes(x = centered.calf, y = std.NB.slope),
#              color = "gold", size = 2) 
