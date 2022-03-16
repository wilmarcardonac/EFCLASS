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

# LCDM

#Cl = np.loadtxt('./output/Cl_lcdm_cl.dat',unpack=True)
#Clp = np.loadtxt('./output/Cl_lcdm_pk.dat',unpack=True)
#back = np.loadtxt('./output/Cl_lcdm_background.dat',unpack=True)

# LCDM_NU

Cl = np.loadtxt('./output/Cl_lcdm_00_cl.dat',unpack=True)
Clp = np.loadtxt('./output/Cl_lcdm_00_pk.dat',unpack=True)
back = np.loadtxt('./output/Cl_lcdm_00_background.dat',unpack=True)

# ewCDM

Cl5 = np.loadtxt('./output/Cl_ewcdm_00_cl.dat',unpack=True)
Cl5p = np.loadtxt('./output/Cl_ewcdm_00_pk.dat',unpack=True)
back5 = np.loadtxt('./output/Cl_ewcdm_00_background.dat',unpack=True)

# FIGURES

#py.plot(-np.log(1.+back2[0]),back2[3],label=r'HS')
#py.plot(-np.log(1.+back[0]),back[3],label=r'$\Lambda CDM$',color='red',ls='dotted')
#py.plot(-np.log(1.+back[0]),abs(1.-back2[3]/back[3])*1.e2,label=r'$\% diff. HS$')
#py.plot(-np.log(1.+back[0]),abs(1.-back3[3]/back[3])*1.e2,label=r'$\% diff. Granda$')
#py.plot(-np.log(1.+back4[0]),back4[3],label=r'$HORNDESKI DES$',color='blue',ls='dashed')
#py.plot(-np.log(1.+back5[0]),back5[3],label=r'$ewcdm$',color='blue',ls='dashed')
#py.xlabel(r'$x$')
#py.ylabel(r'$H(x)$')
#py.xscale('linear')
#py.legend(loc=0)
#py.savefig('./output/Cl_husawicki_background.pdf')
#py.savefig('./output/Cl_granda_background.pdf')
#py.savefig('./output/Cl_designer_horndeski_background.pdf')
#py.savefig('./output/Cl_ewcdm_background.pdf')
#py.close()

Tfactor = (2.726e6)**2

#TT

#py.loglog(Cl1[0],Cl1[1]*Tfactor,label=r'$w=-1$ DES')
#py.loglog(Cl1a[0],Cl1a[1]*Tfactor,label=r'$w=-1$ DES Approx')
#py.loglog(abs(Cl1a[0]-Cl1[0])/Cl1[0],Cl1[1]*Tfactor,label=r'$w=-1$ DES')
#py.loglog(Cl2[0],Cl2[1]*Tfactor,label=r'HS')
#py.loglog(Cl3[0],Cl3[1]*Tfactor,label=r'Granda')
#py.loglog(Cl4[0],Cl4[1]*Tfactor,label=r'$w=-1$ HORNDESKI DES')
py.loglog(Cl5[0],Cl5[1]*Tfactor,label=r'$e_\pi wCDM$')
py.loglog(Cl[0],Cl[1]*Tfactor,label=r'$\Lambda CDM$')
py.xlabel(r'$\ell$')
py.ylabel(r'$\ell(\ell+1)C_\ell/2\pi$')
#py.xlim(2,3000)
#py.xlim(2,40)
#py.ylim(750,1050)
#py.yscale('linear')
py.legend(loc=0)
#py.savefig('./output/Cl_designer_lcdm_cl.pdf')
#py.savefig('./output/Cl_husawicki_cl.pdf')
#py.savefig('./output/Cl_granda_cl.pdf')
#py.savefig('./output/Cl_designer_horndeski_lcdm_cl.pdf')
py.savefig('./output/Cl_ewcdm_cl.pdf')
py.close()

# EE

py.loglog(Cl5[0],Cl5[2],label=r'$e_\pi wcdm$')
#py.loglog(Cl4[0],Cl4[2],label=r'$w=-1$ HORNDESKI DES')
#py.loglog(Cl3[0],Cl3[2],label=r'Granda')
#py.loglog(Cl2[0],Cl2[2],label=r'HS')
py.loglog(Cl[0],Cl[2],label=r'$\Lambda CDM$')
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
py.savefig('./output/Cl_ewcdm_EE.pdf')
py.close()

# TE

#py.loglog(Cl3[0],abs(Cl3[3]),label=r'Granda')
#py.loglog(Cl2[0],abs(Cl2[3]),label=r'HS')
#py.loglog(Cl4[0],abs(Cl4[3]),label=r'$w=-1$ HORNDESKI DES')
py.loglog(Cl5[0],abs(Cl5[3]),label=r'$e_\pi wCDM$')
py.loglog(Cl[0],abs(Cl[3]),label=r'$\Lambda CDM$')
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
py.savefig('./output/Cl_ewcdm_cl_TE.pdf')
py.close()

py.loglog(Cl5p[0],Cl5p[1],label=r'$e_\pi wCDM$')
#py.loglog(Cl4p[0],Cl4p[1],label=r'$w=-1$ HORNDESKI DES')
#py.loglog(Cl2p[0],Cl2p[1],label=r'HS')
#py.loglog(Cl1p[0],Cl1p[1],label=r'DES')
#py.loglog(Cl3p[0],Cl3p[1],label=r'Granda')
py.loglog(Clp[0],Clp[1],label=r'$\Lambda CDM$')
#py.xlabel(r'$\ell$')
#py.ylabel(r'$\ell(\ell+1)C_\ell/2\pi$')
#py.xlim(2,3000)
py.xlim(1.e-4,1)
py.ylim(10,3.5e4)
#py.yscale('linear')
py.legend(loc=0)
#py.savefig('./output/Cl_designer_lcdm_pk.pdf')
#py.savefig('./output/Cl_husawicki_pk.pdf')
#py.savefig('./output/Cl_designer_horndeski_lcdm_pk.pdf')
py.savefig('./output/Cl_ewcdm_pk.pdf')
#py.savefig('./output/Cl_granda_pk.pdf')

exit()
