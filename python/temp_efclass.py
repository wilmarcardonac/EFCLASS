from classy import Class
from pylab import *
import numpy as np

# Define your cosmology (what is not specified will be set to CLASS default parameters)
params = {
    'lensing' : 'no',
    'k_pivot' : 0.05,
    'A_s' : 2.3e-9, 
    'n_s': 1., 
    'h' : 0.7,
    'Omega_b' : 0.05,
    'Omega_cdm' : 0.25,
    'z_reio' : 10.,
    'log10b_pi' : -0.797825028,
    'Omega_k' : 0.,
    'Omega_fld' : 0.,
    'output' : 'tCl,pCl',
#    'l_min' : 2,
    'k_output_values' : 1.5e-2,
    'gauge' : 'newtonian',
    'perturb_sampling_stepsize' : 0.01,
    'tol_perturb_integration' : 1.e-9}

# Create an instance of the CLASS wrapper
cosmo = Class()

# Set the parameters to the cosmological code
cosmo.set(params)

# Run the whole code. Depending on your output, it will call the
# CLASS modules more or less fast. For instance, without any
# output asked, CLASS will only compute background quantities,
# thus running almost instantaneously.
# This is equivalent to the beginning of the `main` routine of CLASS,
# with all the struct_init() methods called.
cosmo.compute()

# Access the lensed cl until l=2000
cls = cosmo.density_cl(2000)

# Print on screen to see the output
#print cls

# plot something with matplotlib...

ell=arange(0,len(cls[0]),1)
plot(ell,ell*(ell+1.)/2/pi*cls[0])
show()

#with open("../output/test_cl.dat") as f:
#    data = f.read()
#plot(data)
#show()


# Clean CLASS (the equivalent of the struct_free() in the `main`
# of CLASS. This step is primordial when running in a loop over different
# cosmologies, as you will saturate your memory very fast if you ommit
# it.
cosmo.struct_cleanup()

# If you want to change completely the cosmology, you should also
# clean the arguments, otherwise, if you are simply running on a loop
# of different values for the same parameters, this step is not needed
cosmo.empty()
