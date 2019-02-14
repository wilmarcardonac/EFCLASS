from classy import Class
from pylab import *
import numpy as np

# Define your cosmology (what is not specified will be set to CLASS default parameters)
params = {
    'output': 'nCl',
    'number count contributions': 'density, rsd, lensing',
    'non linear': 'halofit',
    'selection': 'gaussian',
    'dNdz_selection': 'analytic',
    'selection_mean':  '0.2805,  0.6755,  0.9055,  1.174,   2.1635',
    'selection_width':  '0.2805,  0.1145,  0.1155,  0.153,   0.8365',
    'selection_bias': '1.1316,   1.2944,  1.3804,  1.4744,  1.7786',
    'non_diagonal': 4,
    'l_max_lss': 2000,
    'A_s': 2.0909e-9,
    'n_s': 0.96276, 
    'h': 0.73094,
    'omega_b': 0.02975,
    'omega_cdm': 0.133346,
    'l_switch_limber_for_cl_density_over_z': 10000.,
    'selection_sampling_bessel': .3,
    'q_linstep': 1000.,
    'k_max_tau0_over_l_max': 2.}

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
#cls = cosmo.density_cl(2000)

# Print on screen to see the output
#print cls

# plot something with matplotlib...

#ell=arange(0,len(cls[0]),1)
#plot(ell,ell*(ell+1.)/2/pi*cls[0])
#show()

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
