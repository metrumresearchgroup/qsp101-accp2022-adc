// Explain your model up here if you would like to do so
// Simple 1 cmt PK model written as ODEs

[PARAM]
 CL = 0.02
 VC = 0.5
 KA = 0.9

[CMT]
 GUT
 CENTRAL

[ODE]
dxdt_GUT = -KA*GUT;
dxdt_CENTRAL = KA*GUT - (CL/VC)*CENTRAL;

[TABLE]
 capture CP = CENTRAL/VC; 

[CAPTURE]
 CL