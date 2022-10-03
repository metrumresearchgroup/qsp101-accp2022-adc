setwd("~/training/accp-qsp-2022")
source("renv/activate.R")

library(mrgsolve)
library(tidyverse)
library(ggplot2)
library(here)

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
#Koff <- Kd*Kon*6.022e23/1e9 # Convert back to nM


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
  Kon = Kon, Kd = Kd, #Koff = Koff, 
  Kdeg = Kdeg, Kint = Kint, Klys = Klys, Kend = Kend, 
  Krec = Krec, Krec_AR = Krec_AR,
  DAR = DAR, Peff = Peff, Nc0 = Nc0, Nr=Nr
)

mod <- mread(here("model/QSP_ADCmodel_Bystander.mod"), project = ".")

# Volume of media in well: 0.5 mL
Vm <- 0.5/1000  # Convert to L
params$Emax_Payload = 0.000038
params$EC50_Payload = 0.8*0.432
mod <- mod %>% mrgsolve::param(params)
conc = 10 # 10 nM dose
e <- ev(amt=conc*Vm/1e9*6.022e23)
sim <- mod %>% mrgsolve::mrgsim(events=e, delta=10, end=3600*8.0, out="df") %>%
  mutate(hours = time/3600)



## Plot results

### Compare Simulated intact ADC concentration in media with observed concentration
obs = read_csv(here("data/observations.csv"),show_col_types=FALSE)

ggplot() + geom_point(data=obs, aes(x=hours,y=A_m/6.022e23*1e9/Vm),inherit.aes=FALSE) + 
  geom_line(data=sim, aes(x=hours,y=A_m/6.022e23*1e9/Vm)) +
  xlab("Time (hours)") + ylab("Free intact ADC in Media (nM)")

ggplot(data=sim, aes(hours, R_s)) + geom_line() + xlab("Time (hours)") + 
  ylab("HER3 Surface Receptor Expression (number of receptors)") + theme_bw()

ggplot(data=sim, aes(hours, AR_s)) + geom_line() + xlab("Time (hours)") + 
  ylab("ADC-HER3 Surface Complex") + theme_bw()

ggplot(data=sim, aes(hours, Ntot)) + geom_line() + xlab("Time (hours)") + 
  ylab("Total Cell (+) Count") + theme_bw()

ggplot(data=sim, aes(hours, AR_e)) + geom_line() + xlab("Time (hours)") + 
  ylab("ADC-HER3 Receptor Complex in Endosome") + theme_bw()

ggplot(data=sim, aes(hours, A_e)) + geom_line() + xlab("Time (hours)") + 
  ylab("Free ADC in Endosome") + theme_bw()

ggplot(data=sim, aes(hours, R_e)) + geom_line() + xlab("Time (hours)") + 
  ylab("Free HER3 in Endosome") + theme_bw()

ggplot(data=sim, aes(hours, AR_l)) + geom_line() + xlab("Time (hours)") + 
  ylab("ADC-HER3 Complex in Lysozome") + theme_bw()

ggplot(data=sim, aes(hours, A_l)) + geom_line() + xlab("Time (hours)") + 
  ylab("Free ADC in Lysozome") + theme_bw()

ggplot(data=sim, aes(hours, R_l)) + geom_line() + xlab("Time (hours)") + 
  ylab("Free HER3 in Lysozome") + theme_bw()

ggplot(data=sim, aes(hours, C_P_c)) + geom_line() + xlab("Time (hours)") + 
  ylab("Deconjugated Payload Concentration in Cytosol (nM)") + theme_bw()

ggplot(data=sim, aes(hours, C_P_m)) + geom_line() + xlab("Time (hours)") + 
  ylab("Deconjugated Payload Concentration in Media (nM)") + theme_bw()

# Target (-) Cells
ggplot(data=sim, aes(hours, C_P_cneg)) + geom_line() + xlab("Time (hours)") + 
  ylab("Deconjugated Payload Concentration in Cytosol (Target -) (nM)") + theme_bw()

ggplot(data=sim, aes(hours, Ntot_Neg)) + geom_line() + xlab("Time (hours)") + 
  ylab("Total Cell (-) Count") + theme_bw()


