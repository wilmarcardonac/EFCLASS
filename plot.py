import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as py

# Loading data files into arrays 

#CMB

# DES

Cl1 = np.loadtxt('./output/Cl_designer_lcdm_cl.dat',unpack=True)
Cl1a = np.loadtxt('./output/Cl_designer_lcdm_cl_ap.dat',unpack=True)
# Hu-Sawicky

#Cl2 = np.loadtxt('./output/Cl_husawicki_cl.dat',unpack=True)

# FIGURE 

Tfactor = (2.726e6)**2

py.loglog(Cl1[0],Cl1[1]*Tfactor,label=r'$w=-1$ DES')
py.loglog(Cl1a[0],Cl1a[1]*Tfactor,label=r'$w=-1$ DES Approx')
#py.loglog(abs(Cl1a[0]-Cl1[0])/Cl1[0],Cl1[1]*Tfactor,label=r'$w=-1$ DES')
#py.loglog(Cl2[0],Cl2[1]*Tfactor,label=r'HS')
py.xlabel(r'$\ell$')
py.ylabel(r'$\ell(\ell+1)C_\ell/2\pi$')
#py.xlim(2,3000)
#py.xlim(2,40)
#py.ylim(620,1200)
#py.yscale('linear')
py.legend(loc=0)
py.savefig('./output/Cl_designer_lcdm_cl.pdf')
#py.savefig('./output/Cl_husawicki_cl.pdf')
py.close()

exit()
