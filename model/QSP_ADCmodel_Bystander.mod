//
// In vitro antibody-drug conjugate disposition model targeting HER3 receptors
//
// Includes:
//   - Receptor turnover
//   - Cell killing dynamics
//   - Payload escape and diffusion
//
// Assumptions:
//   - Deconjugation only happens in endosome/lysosome
//
// Units: Liters, hours, nanomoles
//


[PARAM]

// Compartment volumes and surface areas

Vm   = 5e-4      // Media volume
Vc   = 3.68e-12  // Volume of single cell
Sc   = 1.66e-5   // Surface area of a single cell (cm^2)
Nc0  = 1.5e5     // Initial number of cells in well

// Rate constants

Kon  = 0.0  // ADC/receptor on rate constant
// Koff = 1.0  // ADC/receptor off rate constant
Kd   = 0.3   // ADC/receptor dissociation rate 
Kdeg = 1.0  // Surface-bound receptor degradation rate
Kint = 0.0  // ADC/receptor complex internalization rate constant
Kend = 0.0  // Endosomal sorting rate for internalized ADC/receptor
Klys = 0.0  // Lysosomal degradation rate constant for internalized ADC/receptor
Krec = 0.0  // Rate for receptors recycling to surface
Krec_AR = 0.0 // Rate of AR complex recycling to surface

Kgrow = 0.0 // Cell growth rate
Kkill = 0.0 // Baseline cell killing/death rate
tau   = 1.0 // Cell death transit compartment time constant

Nmax  = 1e6 // Max number of cells in well (confluence number)

// Other parameters

Peff = 12   // Payload cell membrane permeability, 1e-6 cm/s
DAR  = 8    // Drug/antibody ratio
fu_m = 1    // Payload fraction unbound in media

// Killing effect model
EC50_Payload = 1.0  // Killing due to payload
Emax_Payload = 0.0

EC50_ADCC    = 1.0  // Killing due to ADCC
Emax_ADCC    = 0.0
Nr = 20e3 // Surface receptor expression (receptors/cell)



[MAIN]

R_s_0 = Nr; // Initial cell surface receptor expression

// Initial number of cells in well
Nc_1_0  = Nc0; // All cells are healthy
Nc_2_0  = 0.0; 
Nc_3_0  = 0.0; 
Nc_4_0  = 0.0; 

NcNeg_1_0 = Nc0/2; // Assume half as many target (-) cells as healthy cells

// Calculate receptor synthesis rate from steady-state surface receptor expression assuming recycling and receptor-mediated endocytosis
double Ksyn = Nr * Kdeg;  
  
// Calculate initial number of free receptor in endosomes
R_e_0 = Ksyn/Kend;

// Get Koff from Kd and Kon
double Koff = Kd*Kon*6.022e23/1e9 ; // Convert back to nM
  

[CMT] 

//-------------------------------------------------------------------
// Species 
//-------------------------------------------------------------------

// Media
A_m   // ADC in media
P_m   // Payload in media

// Cell surface
R_s   // Cell surface receptors
AR_s  // Surface-bound ADC/receptor complex

// Endosome
AR_e  // Endosomal ADC/receptor complex
A_e   // Endosomal ADC free ADC
R_e   // Endosomal free receptor

// Lysosome
AR_l  // Lysosomal ADC/receptor complex
A_l   // Lysosomal free ADC
R_l   // Lysosomal free receptor

// Cytoplasm  
P_c   // Payload in cytoplasm
P_cneg   // Payload in cytoplasm of target negative cells

// Catabolized
A_cat // Amount of ADC (antibody) catabolized

// Cell counts
// Target Positive Cells
Nc_1   // Number of (antigen-positive) proliferating n cells 
Nc_2   // Transit (non-growing) compartment
Nc_3   // Transit (non-growing) compartment
Nc_4   // Transit (non-growing) compartment

// Target Negative Cells
NcNeg_1   // Number of (antigen-negative) proliferating n cells 
NcNeg_2   // Transit (non-growing) compartment
NcNeg_3   // Transit (non-growing) compartment
NcNeg_4   // Transit (non-growing) compartment



  
[ODE]

// Total number of cells in well
double Ntot = Nc_1 + Nc_2 + Nc_3 + Nc_4; // Target (+)
double Ntot_Neg = NcNeg_1 + NcNeg_2 + NcNeg_3 + NcNeg_4; // Target(-)

double C_P_c = P_c/(Vc);
double C_P_m = P_m/Vm;
double C_P_cneg = P_cneg/Vc; // Payload concentration in cytoplasm of target negative cells

//-------------------------------------------------------------------
// Fluxes
//-------------------------------------------------------------------

// ADC in media binding to surface receptor
double flux_A_R_s_binding = A_m/Vm*R_s*Kon;

// Surface ADC/receptor unbinding
double flux_AR_s_unbinding = AR_s*Koff;

// Unbound surface receptor degrading
double flux_R_s_degrade = R_s*Kdeg;

// ADC/receptor complex internalizing to endosome
double flux_AR_s_int = AR_s*Kint;

// Surface receptor synthesis and feedback
double flux_R_s_syn = Ksyn;

// Endosomal ADC/receptor unbinding
double flux_AR_e_unbinding = AR_e*Koff;

// Endosomal receptor recycles to surface
double flux_R_e_recycle = R_e*Krec;

// Endosomal receptor transport to lysosome
double flux_R_e_to_R_l = R_e*Kend;

// Endosomal unbound ADC transport to lysosome
double flux_A_e_to_A_l = A_e*Kend;

// Endosomal ADC/receptor complex transport to lysosome
double flux_AR_e_to_AR_l = AR_e*Kend;

// Endosomal ADC/receptor recycles to surface
double flux_AR_e_recycle = AR_e*Krec_AR;

// Lysosomal ADC/receptor complex catabolized
double flux_AR_l_cat = AR_l*Klys;

// Lysosomal antibody catabolized
double flux_A_l_cat = A_l*Klys;

// Lysosomal receptor catabolized
double flux_R_l_cat = R_l*Klys;

// Payload escape from lysosome to cytoplasm
double flux_P_l_to_P_c = DAR*A_l*Klys;

// Payload diffusion from cytoplasm to media (can be negative)
double flux_P_c_to_P_m = (3600*Peff*1e-6)*(C_P_c - fu_m*C_P_m)*Sc/1000;

// Payload diffusion from media to cytoplasm of target negative cells
double flux_P_m_to_P_cneg = (3600*Peff*1e-6)*(fu_m*C_P_m - C_P_cneg)*Sc/1000;



//-------------------------------------------------------------------
// Derivatives
//-------------------------------------------------------------------

// Receptor occupancies (using baseline receptor abundance)
double RO = AR_s/R_s_0;

// Cell dynamics

// Effective kill rate = baseline death rate + payload killing rate + ADCC killing rate
// ADCC is assumed to be negligible, here.
double Kkill_eff = Kkill + Emax_Payload*(P_c/Vc/6.022e23*1e9)/(EC50_Payload + (P_c/Vc/6.022e23*1e9))
  + Emax_ADCC*RO/(EC50_ADCC + RO);
  
// Killing of target negative cells by payload
double Kkill_effNeg = Kkill + Emax_Payload*(P_cneg/Vc/6.022e23*1e9)/(EC50_Payload + (P_cneg/Vc/6.022e23*1e9));

// Effective growth rate = rate from cell multiplication limited to confluence number
double Kgrow_eff = Kgrow*(1 - Ntot/Nmax);

// Overall cell growth/death
dxdt_Nc_1 = Kgrow_eff*Nc_1 - Kkill_eff*Nc_1;
dxdt_NcNeg_1 = Kgrow_eff*NcNeg_1 - Kkill_effNeg*NcNeg_1;

// Transit compartments (non-growing) for cells in process of being killed
dxdt_Nc_2 = Kkill_eff*Nc_1 - Nc_2/tau;
dxdt_NcNeg_2 = Kkill_eff*NcNeg_1 - NcNeg_2/tau;

dxdt_Nc_3 = (Nc_2 - Nc_3)/tau;
dxdt_NcNeg_3 = (NcNeg_2 - NcNeg_3)/tau;

dxdt_Nc_4 = (Nc_3 - Nc_4)/tau;
dxdt_NcNeg_4 = (NcNeg_3 - NcNeg_4)/tau;

// ADCs in media
// Flux = unbinding - binding
dxdt_A_m = flux_AR_s_unbinding*Ntot - flux_A_R_s_binding*Ntot;

// ADCs bound to cell surface
// Flux = binding - unbinding - internalization
dxdt_AR_s = flux_A_R_s_binding - flux_AR_s_unbinding - flux_AR_s_int
  + flux_AR_e_recycle;

// ADC/receptor complex in endosome
// Flux = internalization - unbinding - transport to lysosome - recycling
dxdt_AR_e = flux_AR_s_int - flux_AR_e_unbinding - flux_AR_e_to_AR_l
  - flux_AR_e_recycle;

// Unbound receptor in endosome
// Flux = unbinding - recycling - transport to lysosome
dxdt_R_e = flux_AR_e_unbinding + flux_R_s_degrade - flux_R_e_recycle - flux_R_e_to_R_l;

// Unbound ADC in endosome
// Flux = unbinding - transport to lysosome
dxdt_A_e = flux_AR_e_unbinding - flux_A_e_to_A_l;

// Unbound receptor in lysosome
// Flux = transport from endosome - catabolism
dxdt_R_l = flux_R_e_to_R_l - flux_R_l_cat;

// Unbound ADC in lysosome
// Flux = transport from endosome - catabolism
dxdt_A_l = flux_A_e_to_A_l - flux_A_l_cat;

// ADC/receptor complex in lysosome
// Flux = transport from endosome - catabolism
dxdt_AR_l = flux_AR_e_to_AR_l - flux_AR_l_cat;

// Unconjugated payload in media
// Flux = diffusion from cytoplasm to media
dxdt_P_m = flux_P_c_to_P_m*Ntot + DAR*Kkill_eff*Nc_1*(A_e + A_l + AR_e + AR_s)/Nc_1 - flux_P_m_to_P_cneg*Ntot_Neg;

// Unconjugated payload in cytoplasm
// Flux = payload escape from lysosome - diffusion from cytoplasm to media
dxdt_P_c = flux_P_l_to_P_c - flux_P_c_to_P_m;

// Unconjugated payload in cytoplasm of target negative cells
// Flux = diffusion from media into cytoplasm
dxdt_P_cneg = flux_P_m_to_P_cneg;


// Free receptors on surface
// Flux = synthesis + unbinding + recycling - binding - degradation
dxdt_R_s = flux_R_s_syn + flux_AR_s_unbinding + flux_R_e_recycle- flux_A_R_s_binding - flux_R_s_degrade;
  
  
  // Total catabolized antibodies
dxdt_A_cat = AR_l*Klys + A_l*Klys;
  
[CAPTURE]
C_A_m = A_m/Vm;
C_R_s = R_s/Vm;
C_P_m = P_m/6.022e23*1e9/Vm;
C_P_c = P_c/6.022e23*1e9/Vc;
C_P_cneg = P_cneg/6.022e23*1e9/Vc;
C_AR_s = AR_s/Vm;
Ntot;
Ntot_Neg;
flux_AR_e_to_AR_l;
flux_AR_l_cat;

  
  