import matplotlib.pyplot as plt
import numpy as np
#------------------------------------------------------------------------
# function to calculate T, BETA, and specifies LIS
def MOMEN(P, E0, T, AZ, BETA, FNP, FLIS):
    BETA = P/np.sqrt(P*P+E0*E0*AZ*AZ)
    FNP = BETA*P
    T = P/BETA/AZ - E0
    FLIS = 21.1*T**(-2.8)/(1. + 5.85*T**(-1.22) + 1.18*T**(-2.54))
    return T, BETA, FNP, FLIS

#------------------------------------------------------------------------
# Program units:
# P = rigidity in GV
# T = kinetic energy in GeV/nuc
# R = radial distance, AU
# DK = (radial part) of diffusion coefficient, 
# V = solar wind speed, units if 400km/s
# speed of light = 750 program units

N = 91 # Number of grid points in the radial direction
D = 1 # The grid size in the r-direction, delta_r
PMIN = .02 # The minimum rigidity where the model will terminate
P = 20. # The maximum rigidity where the model will start, PMAX
DLNP = .02 # The grid size in the ln(P)-direction
CK1 = 73. # The value of kappa_0, the reference value at Earth
# The program units for kappa is 6 \times 10^{20} cm^2/s
CK2 = 0. # The exponent that determines the radial dependence of kappa
E0 = .938 # The rest mass energy in GeV
AZ = 1. # The species value of the cosmic ray population
BETA = 0. # Initiate the v/c ratio
FNP = 0. # Shorthand variable BETA*P = v/c*P
FLIS = 0. # The LIS value in units of differential intensity
T = 0. # Kinetic energy per nucleus, in GeV/nuc

DK = np.zeros(N + 1) # Stores radial dependence of diffusion coefficient
DKDR = np.zeros(N + 1) # Radial derivative of diffusion coefficient
V = np.zeros(N + 1) # Radially dependend solar wind speed
DV = np.zeros(N + 1) # Divergence of solar wind speed
R = np.zeros(N + 1) # The radial coordinate
F = np.zeros(N + 1) # The computed distribution function as a function of r
X_POS = np.zeros(N + 1) # Placeholder variable for calculations
X_NEG = np.zeros(N + 1) # Placeholder variable for calculations
Y = np.zeros(N + 1) # Placeholder variable for calculations
     
      
f = open('output.txt',"w")
      
T, BETA, FNP, SPECTRUM = MOMEN(P, E0, T, AZ, BETA, FNP, FLIS)
# SPECTRUM = LIS at PMAX

R[0]=-D+0.000001
# First radial coordinate is ghost point at r = -D

for i in range(1, N + 1):
    F[i] = SPECTRUM/P/P/AZ
    # The LIS value converted to a distribution funciton
    R[i] = R[i-1] + D
    # Setting uo the radial coordinate
    V[i] = 1.-np.exp(-13.862*R[i])
    # Specify the solar wind speed
    DV[i] = 2.*V[i]/R[i]+(1.-V[i])*13.862
    # Divergence of the solar wind speed
    DK[i] = CK1*R[i]**CK2 
    # Radial dependence of diffusion coefficient
    DKDR[i] = DK[i]*CK2/R[i] 
    # Radial derivate of radial dependence of diffusion coefficient

while(P > PMIN):
# Model will step in log(rigidity) from PMAX to PMIN
    F[0]=F[2]
    # Reflective boundary condition at radial boundary
    
    for i in range(0,N + 1):
        # Placeholder variables fro calculations
        X_NEG[i] = 0.
        X_NEG[1] = -1.
        
    T, BETA, FNP, SPECTRUM = MOMEN(P/np.exp(DLNP/2.), E0, T, AZ, BETA, FNP, FLIS)
    # Calculate LIS at half rigidity step
    
    for i in range(1, N):
        # Solve for r-equation
        A = FNP*(DK[i]/D/D/2.-(2.*DK[i]/R[i]+DKDR[i])/D/4.)+V[i]/D/4.
        B = FNP*DK[i]/D/D - A
        C = -A -B -A*X_POS[i-1] - DV[i]/DLNP/3.
        X_POS[i]=(B - A*X_NEG[i])/C # Should this move up one line??
        Y[i]=(-A*F[i-1] - (DV[i]/DLNP/3. -A -B)*F[i]-B*F[i + 1]- A*Y[i-1])/C
    
    F[N - 1] = Y[N - 1] - X_POS[N - 1]*SPECTRUM/((P/np.exp(DLNP/2.))**2/AZ)
    T, BETA, FNP, SPECTRUM = MOMEN(P/np.exp(DLNP),E0,T,AZ,BETA,FNP,SPECTRUM)
    F[N] = SPECTRUM/(P/np.exp(DLNP))**2/AZ
    
    for i in range(N-1,0,-1):
        # Solve for P-equation
        F[i]=Y[i]-X_POS[i]*F[i+1]
        DP2=(P/np.exp(DLNP))**2
        #DP2 = new rigidity value (squared)

    write_text = str(T*1000) + ',' + str(DP2*F[2]) + ',' + str(DP2*F[81]) + ',' + str(DP2*F[N])   + '\n'
    f.write(write_text)
    
    P=P/np.exp(DLNP)
    # Step down in rigidity
      
f.close()      
#-------------------------------------------------------------
