setwd("~/training/accp-qsp-2022")
source("renv/activate.R")

#' --- 
#' title: Sensitivity Analysis 
#' ---
library(mrgsolve)
library(ggplot2)
library(here)
library(mrgsim.sa)
library(patchwork)

# Setup some initial parameters

# Cell volume: 1760 um^3 
Vc <- 1760*((1e-4)^3)/1000 # um^3 -> cm^3(mL) -> L

# Cell Surface area: 1663 um^2
Sc <- 1663*((1e-4)^2)      # um^2 -> cm^2

# Volume of media in well: 0.5 mL
Vm <- 0.5/1000  # Convert to L

#
# Drug parameters
#

Kd <- 0.327  # nM 

# Kon values, taken from an anti-HER3 mAbs:
# https://journals.sagepub.com/doi/pdf/10.1177/1533034615588422: 5.66e5 (M^-1 s^-1)
# https://pubmed.ncbi.nlm.nih.gov/31911530/: 1.99e5 (M^-1 s^-1)
# Use an average of the two values
Kon <- mean(5.66e5, 1.99e5) * 3600/6.022e23 # (molecules^-1 h^-1)

# Get Koff from Kd and Kon
Koff <- Kd*Kon*6.022e23/1e9 # Convert back to nM


# Lysosomal degradation rate based on endosome half-life
# Reported here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2112524/pdf/jc972508.pdf
Endo_thalf <- 7.5 # min
Kend = (log(2)/Endo_thalf)*60 # Convert to 1/h
Klys = (log(2)/Endo_thalf)*60 # Convert to 1/h

#
# Receptor dynamics
#

# HER3 half-life, from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6205365/
# 4.8 hours
R_thalf <- 4.8

# Convert to degradation rate constant
Kdeg <- log(2)/R_thalf

# in the .mod file, we use SS assumption under initial conditions to 
# calculate receptor synthesis rate (Ksyn), i.e. - from steady-state surface receptor 
# expression assuming recycling and receptor-mediated endocytosis
## // double Ksyn = Nr * Kdeg;  ///(Krec + Kend) * Kend*1;

# Half-life value of internalization
Kint = log(2)/0.5


# Rate of internalized receptors recycling to surface.  
# Manually adjusted as fraction of Kdeg or Kint by fit to in vitro data
Krec    <- Kdeg/1.5 
Krec_AR <- Kint/1.5


Nc0 = 1e5 # Number of cells seeded in wells
Nr = 20e3 # Initial surface HER3 receptor expression

# Membrane permeability, from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4946713/
# PAMPA Peff, 10-6 cm/s: 13 @ pH = 5.0, 12.2 @ pH = 7.4
Peff <- 12.2

# Drug/antibody ratio
DAR <- 8

params <- list(
  Vc = Vc, Sc = Sc, Vm = Vm,
  Kon = Kon, Kd = Kd, # Koff = Koff, 
  Kdeg = Kdeg, Kint = Kint, Klys = Klys, Kend = Kend, 
  Krec = Krec, Krec_AR = Krec_AR,
  DAR = DAR, Peff = Peff, Nc0 = Nc0, Nr=Nr
)

mod <- mread(here("model/QSP_ADCmodel.mod"), project = ".")

# Volume of media in well: 0.5 mL
Vm <- 0.5/1000  # Convert to L
params$Emax_Payload = 0.000038
params$EC50_Payload = 0.8*0.432
mod <- param(mod, params)
conc = 10 # 10 nM dose
e <- ev(amt=conc*Vm/1e9*6.022e23)
sim <- mod %>% mrgsim(events=e, delta=10, end=3600*8.0, out="df")


#' Update model object
mod <- update(mod, delta = 10, end = 3600*6, rtol = 1e-4, atol = 1e-10)

#' Check the parameters ... we can tweak this with 
#' sensitivity analysis functions
param(mod)

#' Look at `Nr`, `Koff`, and `Peff`
#' `?parseq_fct`
#' `.factor` defaults to 2 so this will fill in values between half and double
#' the model parameter value.
sens1 <- 
  mod %>% 
  ev(e) %>% 
  parseq_fct(Nr, Emax_Payload, Kon, Kd,.factor = 2) %>% # for Kd, could also look at that as factor of dose
  sens_each() %>%
  mutate(time = time/3600)

sens_plot(sens1, "Ntot", ncol = 2,
          xlab = "Time (h)",
          ylab = "Total Cell Count")

sens100 <- 
  mod %>% 
  ev(e) %>% 
  parseq_fct(Nr,Kd, Kon,.factor = 100) %>% 
  sens_each() %>%
  mutate(time = time/3600)

sens_plot(sens100, "Ntot", ncol = 2,
          xlab = "Time (h)",
          ylab = "Total Cell Count")


sens_plot(sens1, "C_P_m", grid = FALSE,
          xlab = "Time (h)",
          ylab = "Total Cell Count")

#' Just define the parameter range using a different method (30% CV)
sens2 <- 
  mod %>% 
  ev(e) %>% 
  parseq_cv(Nr, Kd, Emax_Payload) %>% 
  sens_each() %>%
  mutate(time = time/3600)

sens_plot(sens2, "Ntot",
          xlab = "Time (h)",
          ylab = "Total Cell Count")
a <- sens_plot(sens2, "Ntot", grid = FALSE,
               xlab = "Time (h)",
               ylab = "Total Cell Count"); a

#' Local sensitivity analysis
lsa(mod %>% ev(e), par = "Nr,Kd,Emax_Payload", var = "Ntot", eps = 0.1) %>% 
  lsa_plot() / a

#' Focus in on Nr and Kd with larger CV
sens2a <- 
  mod %>% 
  ev(e) %>% 
  parseq_cv(Nr, Kd, .cv = 200) %>% 
  sens_each() %>%
  mutate(time = time/3600)

sens_plot(sens2a, "Ntot",
          xlab = "Time (h)",
          ylab = "Total Cell Count")
a2 <- sens_plot(sens2a, "Ntot", grid = FALSE,
               xlab = "Time (h)",
               ylab = "Total Cell Count"); a2

lsa(mod %>% ev(e), par = "Kd,Nr", var = "Ntot", eps = 400) %>% 
  lsa_plot() / a2

#' There are other techniques to consider, e.g., global sensitivity using Sobol, 
#' with several R packages that can be used. ... that's QSP 201!!

