import matplotlib.pyplot as plt
import numpy as np
#-------------------------------------------------------------
# The 1D modulation model is contained in the function ONEDMODMODEL, in the 
# file Moraal_model.py
from Moraal_model import ONEDMODMODEL
#-------------------------------------------------------------
# Parameters
# KE = Energies where J_1AU and J_LIS are calculated in GeV/nuc
# J_1AU = Differential intensity at 1 AU
# J_LIS = Differential intensity of LIS (120 AU)
# LAMBDA = The effective radial mean-free-path in units of AU at 1 AU
# CK3 = The exponent that determines the rigidity dependence of kappa

LAMBDA = 0.25 # The effective radial mean-free-path in units of AU at 1 AU
CK3 = 1.5 # The exponent that determines the rigidity dependence of kappa
#------------------------------------------------------------------------

KE, J_1AU, J_LIS = ONEDMODMODEL(LAMBDA, CK3)

#------------------------------------------------------------------------
# file output to file
f = open('output.txt',"w")

for i in range(0, len(KE) - 1):
    write_text = str(KE[i]*1000) + ',' + str(J_1AU[i]) + ',' + str(J_LIS[i])   + '\n'
    f.write(write_text)
f.close()      
#-------------------------------------------------------------

fig = plt.figure(figsize = (15,10))

subplot1 = fig.add_subplot(111)
subplot1.set_xlim(1,1e4)
subplot1.set_ylim(1e-2, 1e2)
#subplot1.set_title('Parker Solar Probe', fontsize = 14)
subplot1.set_xlabel('Kinetic energy (MeV)', fontsize = 16)
subplot1.set_ylabel('Differential intensity', fontsize = 16)
subplot1.set_yscale('log')
subplot1.set_xscale('log')
subplot1.tick_params(labelsize = 14)

subplot1.plot(KE*1000.,J_1AU, color = 'red', label = 'Intensity @ Earth (r = 1 AU)')
subplot1.plot(KE*1000.,J_LIS, color = 'black', label = 'LIS (r = 120 AU)', linewidth = 2)

subplot1.legend(ncol = 2, fontsize = 14)
      
fig.savefig('output.png', dpi=300, bbox_inches='tight')

plt.show()
