# 1DModModel
1D cosmic ray transport model from Moraal and co-workers, ca 1980 translated into Python3.

## Model details

The FORTRAN version of the 1D cosmic ray modulation model of Caballero-Lopez & Moraal (2004) is translated into Python3. The model can be run from the command line:

fskrdts@fskrdts:/Test_3 - Moraal comparison$ python Moraal_model.py 

Followed by the plotting function:

fskrdts@fskrdts:~/Test_5$ python plot_output.py 

## Changes and updates

Model is changed so that the effective radial mean-free-path (LAMBDA) is specified. 

The rigidity dependence of the diffusion coefficient is then approximated as ~P^CK3, where CK3 can also be treated as a free paramater. 

The LIS is now specified at 120 AU and by default the model only calculates the intensty at Earth (1 AU).

## References

[Caballero-Lopez, R.A. & Moraal, H., 2004, Limitations of the force field equation to describe cosmic ray modulation, Journal of Geophysical Research: Space Physics, Volume 109, Issue A1, CiteID A01101](https://ui.adsabs.harvard.edu/abs/2004JGRA..109.1101C/abstract)
