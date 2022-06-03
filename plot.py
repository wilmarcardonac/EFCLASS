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
#back2 = np.loadtxt('./output/Cl_husawicki_background.dat',unpack=True)

# Granda

#Cl3 = np.loadtxt('./output/Cl_granda_cl.dat',unpack=True)
#Cl3p = np.loadtxt('./output/Cl_granda_pk.dat',unpack=True)
#back3 = np.loadtxt('./output/Cl_granda_background.dat',unpack=True)

# DESIGNER HORNDESKI

#Cl4 = np.loadtxt('./output/Cl_designer_horndeski_lcdm_cl.dat',unpack=True)
#Cl4p = np.loadtxt('./output/Cl_designer_horndeski_lcdm_pk.dat',unpack=True)
#back4 = np.loadtxt('./output/Cl_designer_horndeski_lcdm_background.dat',unpack=True)

# DESIGNER SVT 

Cl5 = np.loadtxt('./output/Cl_designer_svt_lcdm_00_cl.dat',unpack=True)
Cl5p = np.loadtxt('./output/Cl_designer_svt_lcdm_00_pk.dat',unpack=True)
back5 = np.loadtxt('./output/Cl_designer_svt_lcdm_00_background.dat',unpack=True)

Cl6 = np.loadtxt('./output/Cl_designer_svt_lcdm_01_cl.dat',unpack=True)
Cl6p = np.loadtxt('./output/Cl_designer_svt_lcdm_01_pk.dat',unpack=True)

Cl7 = np.loadtxt('./output/Cl_designer_svt_lcdm_02_cl.dat',unpack=True)
Cl7p = np.loadtxt('./output/Cl_designer_svt_lcdm_02_pk.dat',unpack=True)

#Cl8 = np.loadtxt('./output/Cl_designer_svt_lcdm_03_cl.dat',unpack=True)
#Cl8p = np.loadtxt('./output/Cl_designer_svt_lcdm_03_pk.dat',unpack=True)

# LCDM

Cl = np.loadtxt('./output/Cl_lcdm_nu_00_cl.dat',unpack=True)
Clp = np.loadtxt('./output/Cl_lcdm_nu_00_pk.dat',unpack=True)
back = np.loadtxt('./output/Cl_lcdm_nu_00_background.dat',unpack=True)

# FIGURES

#py.plot(-np.log(1.+back2[0]),back2[3],label=r'HS')
py.plot(1+back[0],back[3],label=r'$\Lambda CDM$',color='red',ls='dotted')
#py.plot(-np.log(1.+back[0]),abs(1.-back2[3]/back[3])*1.e2,label=r'$\% diff. HS$')
#py.plot(-np.log(1.+back[0]),abs(1.-back3[3]/back[3])*1.e2,label=r'$\% diff. Granda$')
#py.plot(-np.log(1.+back4[0]),back4[3],label=r'$HORNDESKI DES$',color='blue',ls='dashed')
py.plot(1+back5[0],back5[3],label=r'SVTDES',color='blue',ls='dashed')
py.xlabel(r'$1+z$')
py.ylabel(r'$H(z)$')
py.xscale('log')
py.yscale('log')
py.legend(loc=0)
#py.savefig('./output/Cl_husawicki_background.pdf')
#py.savefig('./output/Cl_granda_background.pdf')
#py.savefig('./output/Cl_designer_horndeski_background.pdf')
py.savefig('./output/Cl_designer_svt_background.pdf')
py.close()

Tfactor = (2.726e6)**2

#TT

#py.loglog(Cl1[0],Cl1[1]*Tfactor,label=r'$w=-1$ DES')
#py.loglog(Cl1a[0],Cl1a[1]*Tfactor,label=r'$w=-1$ DES Approx')
#py.loglog(abs(Cl1a[0]-Cl1[0])/Cl1[0],Cl1[1]*Tfactor,label=r'$w=-1$ DES')
#py.loglog(Cl2[0],Cl2[1]*Tfactor,label=r'HS')
#py.loglog(Cl3[0],Cl3[1]*Tfactor,label=r'Granda')
#py.loglog(Cl4[0],Cl4[1]*Tfactor,label=r'$w=-1$ HORNDESKI DES')
py.loglog(Cl5[0],Cl5[1]*Tfactor,label=r'$ \tilde{J} = 10^{-4}$',color='blue',ls='dashed',lw=3)
py.loglog(Cl6[0],Cl6[1]*Tfactor,label=r'$ \tilde{J} = 5\times 10^{-4}$',color='green',ls='dashed',lw=3)
py.loglog(Cl7[0],Cl7[1]*Tfactor,label=r'$ \tilde{J} = 10^{-3}$',color='magenta',ls='dashed',lw=3)
#py.loglog(Cl8[0],Cl8[1]*Tfactor,label=r'SVTDES',color='cyan',ls='dashed',lw=3)
py.loglog(Cl[0],Cl[1]*Tfactor,label=r'$\Lambda$CDM',ls='dotted',color='red',lw=3)
py.xlabel(r'Multipole $\ell$',fontsize=16)
py.ylabel(r'$\ell(\ell+1)C^{TT}_\ell/2\pi \, [\mu \mathrm{K}^2]$',fontsize=16)
py.xticks(fontsize=16)
py.yticks(fontsize=16)
#py.xlim(2,3000)
#py.xlim(2,40)
#py.ylim(750,1050)
#py.yscale('linear')
py.legend(loc='lower left',fontsize='x-large')
#py.savefig('./output/Cl_designer_lcdm_cl.pdf')
#py.savefig('./output/Cl_husawicki_cl.pdf')
#py.savefig('./output/Cl_granda_cl.pdf')
#py.savefig('./output/Cl_designer_horndeski_lcdm_cl.pdf')
py.savefig('./output/Cl_designer_svt_lcdm_cl.pdf')
py.close()

# EE

#py.loglog(Cl4[0],Cl4[2],label=r'$w=-1$ HORNDESKI DES')
py.loglog(Cl5[0],Cl5[2],label=r'SVTDES',color='blue',ls='dashed')
#py.loglog(Cl3[0],Cl3[2],label=r'Granda')
#py.loglog(Cl2[0],Cl2[2],label=r'HS')
py.loglog(Cl[0],Cl[2],label=r'$\Lambda$CDM',color='red',ls='dotted')
py.xlabel(r'$\ell$')
#py.ylabel(r'$\ell(\ell+1)C_\ell/2\pi$')
#py.xlim(2,3000)
#py.xlim(2,40)
#py.ylim(620,1200)
#py.yscale('linear')
py.legend(loc=0)
#py.savefig('./output/Cl_designer_lcdm_cl.pdf')
#py.savefig('./output/Cl_husawicki_EE.pdf')
#py.savefig('./output/Cl_granda_cl_EE.pdf')
#py.savefig('./output/Cl_designer_horndeski_lcdm_EE.pdf')
py.savefig('./output/Cl_designer_svt_lcdm_EE.pdf')
py.close()

# TE

#py.loglog(Cl3[0],abs(Cl3[3]),label=r'Granda')
#py.loglog(Cl2[0],abs(Cl2[3]),label=r'HS')
#py.loglog(Cl4[0],abs(Cl4[3]),label=r'$w=-1$ HORNDESKI DES')
py.loglog(Cl5[0],abs(Cl5[3]),label=r'SVTDES',color='blue',ls='dashed')
py.loglog(Cl[0],abs(Cl[3]),label=r'$\Lambda$CDM',color='red',ls='dotted')
py.xlabel(r'$\ell$')
#py.ylabel(r'$\ell(\ell+1)C_\ell/2\pi$')
#py.xlim(2,3000)
#py.xlim(2,40)
#py.ylim(620,1200)
#py.yscale('linear')
py.legend(loc=0)
#py.savefig('./output/Cl_designer_lcdm_cl.pdf')
#py.savefig('./output/Cl_husawicki_cl_TE.pdf')
#py.savefig('./output/Cl_granda_cl_TE.pdf')
#py.savefig('./output/Cl_designer_horndeski_lcdm_cl_TE.pdf')
py.savefig('./output/Cl_designer_svt_lcdm_cl_TE.pdf')
py.close()

#  P(k)

#py.loglog(Cl4p[0],Cl4p[1],label=r'$w=-1$ HORNDESKI DES')
py.loglog(Cl5p[0],Cl5p[1],label=r'$\tilde{J}=10^{-4}$',color='blue',ls='dashed',lw=3)
py.loglog(Cl6p[0],Cl6p[1],label=r'$\tilde{J}=5\times 10^{-4}$',color='green',ls='dashed',lw=3)
py.loglog(Cl7p[0],Cl7p[1],label=r'$\tilde{J}=10^{-3}$',color='magenta',ls='dashed',lw=3)
#py.loglog(Cl8p[0],Cl8p[1],label=r'SVTDES',color='cyan',ls='dashed',lw=3)
#py.loglog(Cl2p[0],Cl2p[1],label=r'HS')
#py.loglog(Cl1p[0],Cl1p[1],label=r'DES')
#py.loglog(Cl3p[0],Cl3p[1],label=r'Granda')
py.loglog(Clp[0],Clp[1],label=r'$\Lambda$CDM',color='red',ls='dotted',lw=3)
py.xlabel(r'Wavenumber $k\, [h\, \mathrm{Mpc}^{-1}]$',fontsize=16)
py.ylabel(r'$P_m(k) [(h^{-1}\, \mathrm{Mpc})^3]$',fontsize=16)
py.xticks(fontsize=16)
py.yticks(fontsize=16)
#py.xlim(2,3000)
py.xlim(1.e-4,2.)
py.ylim(10.,5.e4)
#py.yscale('linear')
py.legend(loc='lower center',fontsize='x-large')
#py.savefig('./output/Cl_designer_lcdm_pk.pdf')
#py.savefig('./output/Cl_husawicki_pk.pdf')
#py.savefig('./output/Cl_designer_horndeski_lcdm_pk.pdf')
py.savefig('./output/Cl_designer_svt_lcdm_pk.pdf')
#py.savefig('./output/Cl_granda_pk.pdf')

exit()
