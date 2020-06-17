import matplotlib as mpl
mpl.use('Agg')
import numpy as np
import matplotlib.pyplot as py

# Loading data files into arrays 

#CMB

# DES

#Cl1 = np.loadtxt('./output/Cl_designer_lcdm_cl.dat',unpack=True)
#Cl1p = np.loadtxt('./output/Cl_designer_lcdm_pk.dat',unpack=True)

# Hu-Sawicky

#Cl2 = np.loadtxt('./output/Cl_husawicki_cl.dat',unpack=True)
#Cl2p = np.loadtxt('./output/Cl_husawicki_pk.dat',unpack=True)

# Granda

Cl3 = np.loadtxt('./output/Cl_granda_cl.dat',unpack=True)
Cl3p = np.loadtxt('./output/Cl_granda_pk.dat',unpack=True)

Cl = np.loadtxt('./output/Cl_lcdm_cl.dat',unpack=True)
Clp = np.loadtxt('./output/Cl_lcdm_pk.dat',unpack=True)

# FIGURE 

Tfactor = (2.726e6)**2

#TT

#py.loglog(Cl1[0],Cl1[1]*Tfactor,label=r'$w=-1$ DES')
#py.loglog(Cl1a[0],Cl1a[1]*Tfactor,label=r'$w=-1$ DES Approx')
#py.loglog(abs(Cl1a[0]-Cl1[0])/Cl1[0],Cl1[1]*Tfactor,label=r'$w=-1$ DES')
#py.loglog(Cl2[0],Cl2[1]*Tfactor,label=r'HS')
py.loglog(Cl3[0],Cl3[1]*Tfactor,label=r'Granda')
py.loglog(Cl[0],Cl[1]*Tfactor,label=r'$\Lambda CDM$')
py.xlabel(r'$\ell$')
py.ylabel(r'$\ell(\ell+1)C_\ell/2\pi$')
#py.xlim(2,3000)
#py.xlim(2,40)
#py.ylim(620,1200)
#py.yscale('linear')
py.legend(loc=0)
#py.savefig('./output/Cl_designer_lcdm_cl.pdf')
#py.savefig('./output/Cl_husawicki_cl.pdf')
py.savefig('./output/Cl_granda_cl.pdf')
py.close()

# EE

py.loglog(Cl3[0],Cl3[2],label=r'Granda')
py.loglog(Cl[0],Cl[2],label=r'$\Lambda CDM$')
py.xlabel(r'$\ell$')
#py.ylabel(r'$\ell(\ell+1)C_\ell/2\pi$')
#py.xlim(2,3000)
#py.xlim(2,40)
#py.ylim(620,1200)
#py.yscale('linear')
py.legend(loc=0)
#py.savefig('./output/Cl_designer_lcdm_cl.pdf')
#py.savefig('./output/Cl_husawicki_cl.pdf')
py.savefig('./output/Cl_granda_cl_EE.pdf')
py.close()

# TE

py.loglog(Cl3[0],abs(Cl3[3]),label=r'Granda')
py.loglog(Cl[0],abs(Cl[3]),label=r'$\Lambda CDM$')
py.xlabel(r'$\ell$')
#py.ylabel(r'$\ell(\ell+1)C_\ell/2\pi$')
#py.xlim(2,3000)
#py.xlim(2,40)
#py.ylim(620,1200)
#py.yscale('linear')
py.legend(loc=0)
#py.savefig('./output/Cl_designer_lcdm_cl.pdf')
#py.savefig('./output/Cl_husawicki_cl.pdf')
py.savefig('./output/Cl_granda_cl_TE.pdf')
py.close()

#py.loglog(Cl2p[0],Cl2p[1],label=r'HS')
#py.loglog(Cl1p[0],Cl1p[1],label=r'DES')
py.loglog(Cl3p[0],Cl3p[1],label=r'Granda')
py.loglog(Clp[0],Clp[1],label=r'$\Lambda CDM$')
#py.xlabel(r'$\ell$')
#py.ylabel(r'$\ell(\ell+1)C_\ell/2\pi$')
#py.xlim(2,3000)
#py.xlim(2,40)
#py.ylim(620,1200)
#py.yscale('linear')
py.legend(loc=0)
#py.savefig('./output/Cl_designer_lcdm_pk.pdf')
#py.savefig('./output/Cl_husawicki_pk.pdf')
py.savefig('./output/Cl_granda_pk.pdf')

exit()
