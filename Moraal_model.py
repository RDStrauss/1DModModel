#------------------------------------------------------------------------
# Python3 version of the 1D galactic cosmic ray modulation model of 
# Caballero-Lopez and Moraal (2004). The model is called in the main.py file.
# RD Strauss, February 2023, dutoit.strauss@nwu.ac.za

#------------------------------------------------------------------------
import matplotlib.pyplot as plt
import numpy as np
#------------------------------------------------------------------------
# function to calculate T, BETA, and specifies LIS
def MOMENTUM(P, E0, T, AZ, BETA):
    BETA = P/np.sqrt(P*P+E0*E0*AZ*AZ)
    T = P/BETA/AZ - E0
    return T, BETA

#------------------------------------------------------------------------
# function to calculate the LIS 
# FLIS is the LIS in differential intensity. Units will determine units throughout
# The current LIS'se are taken from Bischoff et al. (2019): https://arxiv.org/pdf/1902.10438.pdf
def LIS(P, T, SPECIES, BETA):
    if (SPECIES == 'proton'):
        FLIS = 2620./BETA/BETA*T**(1.1)*((T**(0.98) + 0.7**(0.98))/(1. + 0.7**(0.98)))**(-4.) + 30.*T*T*((T + 8.)/9.)**(-12.)
    if (SPECIES == 'electron'):
        FLIS = 255./BETA/BETA/T*((T + 0.63)/(1.63))**(-2.43) + 6.4*T*T*((T + 15.)/16.)**(-26.)
    if (SPECIES == 'helium'):
        FLIS = 163.4/BETA/BETA*T**(1.1)*((T**(0.97) + 0.58**(0.97))/(1. + 0.58**(0.97)))**(-4.)
    if (SPECIES == 'carbon'):
        FLIS = 3.3/BETA/BETA*T**(1.22)*((T**(0.9) + 0.63**(0.9))/(1. + 0.63**(0.9)))**(-4.43)
    if (SPECIES == 'oxygen'):
        FLIS = 3.3/BETA/BETA*T**(1.23)*((T**(0.86) + 0.62**(0.86))/(1. + 0.62**(0.86)))**(-4.63)
    if (SPECIES == 'boron'):
        FLIS = 1./BETA/BETA*T**(1.7)*((T + 0.685)/(1. + 0.685))**(-4.8) + 3e-4*T**(3)*((T + 0.204)/(1. + 0.204))**(-11)
    if (SPECIES == 'positron'):
        FLIS = 25./BETA/BETA/T**(0.1)*((T**(1.1) + 0.2**(1.1))/(1. + 0.2**(1.1)))**(-3.31) + 23.*T**(0.5)*((T + 2.2)/3.2)**(-9.5)
    return FLIS

#------------------------------------------------------------------------
# function to calculate solar wind speed
def SOLAR_WIND(V, DV, R, N):
    for i in range(0,N + 1):
        V[i] = 1.-np.exp(-13.862*R[i])
        # Specify the solar wind speed
        DV[i] = 2.*V[i]/R[i]+(1.-V[i])*13.862
        # Divergence of the solar wind speed
    return V, DV

#------------------------------------------------------------------------
def ONEDMODMODEL(LAMBDA, CK3, SPECIES):
    # Variables that might be changed
    N = 121 # Number of grid points in the radial direction
    # Using D = 1 and N = 121, puts the heliopause at 120 AU and R[2] at Earth
    D = 1 # The grid size in the r-direction, delta_r
    PMIN = .05 # The minimum rigidity where the model will terminate
    P = 20. # The maximum rigidity where the model will start, PMAX
    DLNP = .02 # The grid size in the ln(P)-direction
    # The program units for kappa is 6 \times 10^{20} cm^2/s
    CK2 = .0 # The exponent that determines the radial dependence of kappa
    E0 = .938 # The rest mass energy of ions in GeV
    #SPECIES = 'protons' # The species value specifies the GCR species under consideration
    AZ = 2. # The species value of the cosmic ray population
    if (SPECIES == 'electron') or (SPECIES == 'proton') or (SPECIES == 'positron'):
        AZ = 1
    # For electrons and/or protons: ZA = 1, for all other (fully ionized) GCR species, AZ = 2
    if (SPECIES == 'electron') or (SPECIES == 'positron'):
        E0 = 0.5110/1000. # Electron trest mass energy in GeV
    print('Species: ', SPECIES, AZ, E0, ' GeV')
    #------------------------------------------------------------------------
    # Program units:
    # P = rigidity in GV
    # T = kinetic energy in GeV/nuc
    # R = radial distance, AU
    # DK = (radial part) of diffusion coefficient, 
    # V = solar wind speed, units if 400km/s
    # speed of light = 750 program units
    
    #------------------------------------------------------------------------
    # Program placeholder variables
    BETA = 0. # Initiate the v/c ratio
    FLIS = 0. # The LIS value in units of differential intensity
    T = 0. # Kinetic energy per nucleus, in GeV/nuc

    #------------------------------------------------------------------------
    # Program variables
    DK = np.zeros(N + 1) # Stores radial dependence of diffusion coefficient
    DKDR = np.zeros(N + 1) # Radial derivative of diffusion coefficient
    V = np.zeros(N + 1) # Radially dependend solar wind speed
    DV = np.zeros(N + 1) # Divergence of solar wind speed
    R = np.zeros(N + 1) # The radial coordinate
    F = np.zeros(N + 1) # The computed distribution function as a function of r
    X_POS = np.zeros(N + 1) # Placeholder variable for calculations
    X_NEG = np.zeros(N + 1) # Placeholder variable for calculations
    Y = np.zeros(N + 1) # Placeholder variable for calculations
    J_1AU = 0. # Differential intensity at Earth
    J_LIS = 0. # Differential intensity at LIS
    T_PRINT = 0. # Energy spectra where J_1AU is calculated
    #------------------------------------------------------------------------
    # Set up radial coordinate
    R[0]=-D+0.000001
    # First radial coordinate is ghost point at r = -D
    for i in range(1, N + 1):
        R[i] = R[i-1] + D

    # Calculate the solar wind speed vs R
    V, DV = SOLAR_WIND(V, DV, R, N)

    # Calculate T, BETA, and BETA*P at given P
    T, BETA = MOMENTUM(P, E0, T, AZ, BETA)
    # Calculate the LIS value at P
    SPECTRUM = LIS(P, T, SPECIES, BETA)

    CK1 = LAMBDA*BETA*750./3.   # The value of kappa_0, the reference value at Earth

    for i in range(1, N + 1):
        F[i] = SPECTRUM/P/P/AZ
        # The LIS value converted to a distribution funciton
        # At P = PMAX, j = J_LIS for all R positions
        DK[i] = CK1*R[i]**CK2 
        # Radial dependence of diffusion coefficient
        DKDR[i] = DK[i]*CK2/R[i] 
        # Radial derivate of radial dependence of diffusion coefficient

    DP2 = P*P
    J_1AU = DP2*F[2]
    J_LIS = DP2*F[N]
    T_PRINT = T
    #------------------------------------------------------------------------
    while(P > PMIN):
    # Model will step in log(rigidity) from PMAX to PMIN
        F[0]=F[2]
        # Reflective boundary condition at radial boundary
    
        X_NEG[1] = -1.
        
        T, BETA = MOMENTUM(P/np.exp(DLNP/2.), E0, T, AZ, BETA)
        SPECTRUM = LIS(P, T, SPECIES, BETA)
        # Calculate LIS at half rigidity step
    
        for i in range(1, N):
            # Solve for r-equation
            A = (P**CK3)*(DK[i]/D/D/2.-(2.*DK[i]/R[i]+DKDR[i])/D/4.)+V[i]/D/4.
            B = (P**CK3)*DK[i]/D/D - A
            C = -A -B -A*X_POS[i-1] - DV[i]/DLNP/3.
            X_POS[i]=(B - A*X_NEG[i])/C # Should this move up one line??
            Y[i]=(-A*F[i-1] - (DV[i]/DLNP/3. -A -B)*F[i]-B*F[i + 1]- A*Y[i-1])/C
    
        F[N - 1] = Y[N - 1] - X_POS[N - 1]*SPECTRUM/((P/np.exp(DLNP/2.))**2/AZ)
    
        T, BETA = MOMENTUM(P/np.exp(DLNP),E0,T,AZ,BETA)
        SPECTRUM = LIS(P, T, SPECIES, BETA)
        F[N] = SPECTRUM/(P/np.exp(DLNP))**2/AZ
    
        for i in range(N-1,0,-1):
            # Solve for P-equation
            F[i]=Y[i]-X_POS[i]*F[i+1]
            DP2=(P/np.exp(DLNP))**2
            #DP2 = new rigidity value (squared)

        J_1AU = np.append(J_1AU, DP2*F[2])
        J_LIS = np.append(J_LIS, DP2*F[N])
        T_PRINT = np.append(T_PRINT, T)

        P=P/np.exp(DLNP)
        # Step down in rigidity
    return T_PRINT, J_1AU, J_LIS 

#------------------------------------------------------------------------

