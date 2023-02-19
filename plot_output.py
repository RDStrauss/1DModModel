import matplotlib.pyplot as plt
import numpy as np
#-------------------------------------------------------------
data = np.loadtxt("output.txt", skiprows = 0, delimiter = ',')

T = data[0:,0]
F = data[0:,1]
FLIS = data[0:,2]

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

subplot1.plot(T,F, color = 'red', label = 'Intensity @ Earth (r = 120 AU)')
subplot1.plot(T,FLIS, color = 'black', label = 'LIS (r = 120 AU)', linewidth = 2)

subplot1.legend(ncol = 2, fontsize = 14)
      
fig.savefig('output.png', dpi=300, bbox_inches='tight')

plt.show()
