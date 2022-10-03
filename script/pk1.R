setwd("~/training/accp-qsp-2022")
source("renv/activate.R")

library(mrgsolve)
library(tidyverse)
library(ggplot2)
library(here)

# Compile
mod <- mread(here("model/pk1.mod"))

# Set intervention
evnt <- ev(amt = 100, ii = 24, addl = 9)

# Simulate
out <- mod %>% 
  ev(evnt) %>%
  mrgsim(end = 480, delta = 0.1)

# Output
out

plot(out)

ggplot(data = as_tibble(out), aes(time, CP)) + geom_line() + xlab("Time (hours)") + 
  ylab("Drug X Concentration (ng/mL)") + theme_bw()


# Explore impact of lower CL
param(mod)

params <- list(
  CL = 0.01)
#mod <- mod %>% mrgsolve::param(params)

out2 <- mod %>% 
  ev(evnt) %>%
  param(params) %>%
  mrgsim(end = 480, delta = 0.1)

# Output
out2

plot(out2)


# Explore and compare the impact of multiple CL values
ievnt <- expand.ev(amt = 100, ii = 24, addl = 9,
          CL = c(0.01, 0.02, 0.04))

out3 <- mod %>% 
  ev(ievnt) %>%
  mrgsim(end = 480, delta = 0.1)

# Output
out3

plot(out3)

ggplot(data = as_tibble(out3), aes(time, CP, color = as.factor(CL))) + 
  geom_line() + 
  xlab("Time (hours)") + 
  ylab("Drug X Concentration (ng/mL)") + 
  theme_bw()

# Want to learn about pre-coded model-library...
?modlib
