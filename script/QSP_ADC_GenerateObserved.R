#' --- 
#' title: Generate pseudo-random "Observed" data
#' ---
library(mrgsolve)
library(ggplot2)
library(here)
library(patchwork)

# Setup some initial parameters

# Cell volume: 1760 um^3 
TVVc <- 1760*((1e-4)^3)/1000 # um^3 -> cm^3(mL) -> L

# Cell Surface area: 1663 um^2
TVSc <- 1663*((1e-4)^2)      # um^2 -> cm^2

# Volume of media in well: 0.5 mL
Vm <- 0.5/1000  # Convert to L

#
# Drug parameters
#

TVKd <- 0.327  # nM 

# Kon values, taken from an anti-HER3 mAbs:
# https://journals.sagepub.com/doi/pdf/10.1177/1533034615588422: 5.66e5 (M^-1 s^-1)
# https://pubmed.ncbi.nlm.nih.gov/31911530/: 1.99e5 (M^-1 s^-1)
# Use an average of the two values
TVKon <- mean(5.66e5, 1.99e5) * 3600/6.022e23 # (molecules^-1 h^-1)

# Get Koff from Kd and Kon
TVKoff <- Kd*Kon*6.022e23/1e9 # Convert back to nM


# Lysosomal degradation rate based on endosome half-life
# Reported here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2112524/pdf/jc972508.pdf
Endo_thalf <- 7.5 # min
TVKend = (log(2)/Endo_thalf)*60 # Convert to 1/h
TVKlys = (log(2)/Endo_thalf)*60 # Convert to 1/h

#
# Receptor dynamics
#

# HER3 half-life, from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6205365/
# 4.8 hours
R_thalf <- 4.8

# Convert to degradation rate constant
TVKdeg <- log(2)/R_thalf

# in the .mod file, we use SS assumption under initial conditions to 
# calculate receptor synthesis rate (Ksyn), i.e. - from steady-state surface receptor 
# expression assuming recycling and receptor-mediated endocytosis
## // double Ksyn = Nr * Kdeg;  ///(Krec + Kend) * Kend*1;

# Half-life value of internalization
TVKint = log(2)/0.5


# Rate of internalized receptors recycling to surface.  
# Manually adjusted as fraction of Kdeg or Kint by fit to in vitro data
TVKrec    <- Kdeg/1.5 
TVKrec_AR <- Kint/1.5


Nc0 = 1e5 # Number of cells seeded in wells
TVNr = 20e3 # Initial surface HER3 receptor expression

# Membrane permeability, from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4946713/
# PAMPA Peff, 10-6 cm/s: 13 @ pH = 5.0, 12.2 @ pH = 7.4
Peff <- 12.2

# Drug/antibody ratio
DAR <- 8

params <- list(
  TVVc = TVVc, TVSc = TVSc, Vm = Vm,
  TVKon = TVKon, TVKd = TVKd, # Koff = Koff, 
  TVKdeg = TVKdeg, TVKint = TVKint, TVKlys = TVKlys, TVKend = TVKend, 
  TVKrec = TVKrec, TVKrec_AR = TVKrec_AR,
  DAR = DAR, Peff = Peff, Nc0 = Nc0, TVNr=TVNr
)

mod <- mread(here("model/QSP_ADCmodel_VAR.mod"), project = ".")

# Volume of media in well: 0.5 mL
Vm <- 0.5/1000  # Convert to L
params$TVEmax_Payload = 0.000038
params$TVEC50_Payload = 0.8*0.432
mod <- param(mod, params)
conc = 10 # 10 nM dose
e <- ev(amt=conc*Vm/1e9*6.022e23)
sim <- mod %>% mrgsim(events=e, delta=10, end=3600*8.0, out="df")


#' Update model object
mod <- update(mod, delta = 10, end = 3600*6, rtol = 1e-4, atol = 1e-10)

e <- expand.ev(amt=conc*Vm/1e9*6.022e23,ID=1:10)

sim <- mod %>% mrgsolve::mrgsim(events=e, delta=3600, end=3600*8.0, out="df") %>%
  mutate(hours = time/3600)

sim <- sim%>% select(ID, time, A_m_out, P_m_out, R_s_out, AR_s_out, AR_e_out, A_e_out, R_e_out, AR_l_out, A_l_out, R_l_out, P_c_out, A_cat_out, Nc_1_out, Nc_2_out, Nc_3_out, Nc_4_out, hours)
sim <- sim%>%rename_at(vars(ends_with("_out")), ~str_replace(., "_out",""))%>%filter(time>0)

sim %>% write_csv(file=here("data/observations.csv"))

