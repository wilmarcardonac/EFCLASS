import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as py

# Loading data files into arrays 

#CMB

# DES

Cl1 = np.loadtxt('./output/Cl_designer_lcdm_cl.dat',unpack=True)

# FIGURE 

Tfactor = (2.726e6)**2

py.loglog(Cl1[0],Cl1[1]*Tfactor,label=r'$w=-1$ DES')
py.xlabel(r'$\ell$')
py.ylabel(r'$\ell(\ell+1)C_\ell/2\pi$')
py.xlim(2,3000)
#py.ylim(620,1200)
#py.yscale('linear')
py.legend(loc=0)
py.savefig('./output/Cl_designer_lcdm_cl.pdf')
py.close()

exit()
